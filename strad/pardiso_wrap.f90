module pardiso_wrap_mod
  implicit none
  save

  include 'mkl_pardiso.fi'

  integer, parameter :: NMAX = 110000
  integer, parameter :: NZMAX = 3100000

  type(MKL_PARDISO_HANDLE) :: pt(64)
  integer :: iparm(64)
  integer :: initialized = 0
  integer :: prev_n = 0
  integer :: prev_nzmod = 0
  integer :: prev_phase_done = 0

  integer :: ia(NMAX+1)
  integer :: ja(NZMAX)
  double precision :: a_csr(NZMAX)

  integer :: rowcount(NMAX)
  integer :: perm(NMAX)
  integer :: solve_count = 0
  integer :: prev_ia_len = 0
  integer :: prev_ja_len = 0
  integer, allocatable :: prev_ia(:), prev_ja(:)

  ! Last-factorized matrix (for identical-matrix factorization reuse).
  integer :: prev_nzmod_out = 0
  integer, allocatable :: fac_ia(:), fac_ja(:)
  double precision, allocatable :: fac_a(:)
  integer :: n_fac = 0, n_reuse = 0   ! factorization vs reuse counters

  ! Diagnostic state (cached once at first call; inert unless env vars set).
  integer :: dbg = 0            ! PARDISO_DEBUG=1 enables per-solve checks
  integer :: dbg_dump_at = -1   ! PARDISO_DUMP_AT=<solve index>
  integer :: diag_read = 0
  double precision :: dbg_dump_res = -1.0d0  ! PARDISO_DUMP_RES=<threshold>

contains

  subroutine pardiso_init_diag()
    ! Read diagnostic env vars once. Keeps them out of the per-solve path.
    character(len=64) :: envv
    integer :: ios
    call get_environment_variable('PARDISO_DEBUG', envv)
    if (len_trim(envv) > 0) read(envv,*,iostat=ios) dbg
    call get_environment_variable('PARDISO_DUMP_AT', envv)
    dbg_dump_at = -1
    if (len_trim(envv) > 0) read(envv,*,iostat=ios) dbg_dump_at
    call get_environment_variable('PARDISO_DUMP_RES', envv)
    dbg_dump_res = -1.0d0
    if (len_trim(envv) > 0) read(envv,*,iostat=ios) dbg_dump_res
    diag_read = 1
  end subroutine pardiso_init_diag

  subroutine pardiso_solve(n, a_coo, irn_coo, icn_coo, nzmod, b, x, ifail)
    integer, intent(in) :: n, nzmod
    integer, intent(in) :: irn_coo(nzmod), icn_coo(nzmod)
    double precision, intent(in) :: a_coo(nzmod), b(n)
    double precision :: b_work(n)
    double precision, intent(out) :: x(n)
    integer, intent(out) :: ifail

    integer :: maxfct, mnum, mtype, phase, nrhs, msglvl, error
    integer :: i, irow, ipos, nzmod_out
    double precision :: resid_coo, resid_csr

    if (diag_read == 0) call pardiso_init_diag()

    ! Sanity checks
    if (n > NMAX) then
      write(*,*) 'PARDISO_WRAP: n=',n,' exceeds NMAX=',NMAX
      ifail = 14
      return
    endif
    if (nzmod > NZMAX) then
      write(*,*) 'PARDISO_WRAP: nzmod=',nzmod,' exceeds NZMAX=',NZMAX
      ifail = 14
      return
    endif

    ! COO -> CSR conversion
    do i = 1, n
      rowcount(i) = 0
    end do
    do i = 1, nzmod
      irow = irn_coo(i)
      if (irow >= 1 .and. irow <= n) then
        rowcount(irow) = rowcount(irow) + 1
      end if
    end do

    ia(1) = 1
    do i = 1, n
      ia(i+1) = ia(i) + rowcount(i)
      rowcount(i) = ia(i)
    end do

    do i = 1, nzmod
      irow = irn_coo(i)
      if (irow >= 1 .and. irow <= n) then
        ipos = rowcount(irow)
        ja(ipos) = icn_coo(i)
        a_csr(ipos) = a_coo(i)
        rowcount(irow) = rowcount(irow) + 1
      end if
    end do

    ! Sort column indices within each row AND merge duplicates
    ! (Y12M sums duplicate (row,col) entries during assembly;
    !  PARDISO requires unique sorted columns per row)
    block
      integer :: j1, j2, itmp_i, wpos
      double precision :: dtmp
      ! Step 1: bubble-sort each row's (ja, a_csr) by column index
      do irow = 1, n
        do j1 = ia(irow), ia(irow+1)-2
          do j2 = j1+1, ia(irow+1)-1
            if (ja(j1) > ja(j2)) then
              itmp_i = ja(j1)
              ja(j1) = ja(j2)
              ja(j2) = itmp_i
              dtmp = a_csr(j1)
              a_csr(j1) = a_csr(j2)
              a_csr(j2) = dtmp
            end if
          end do
        end do
      end do
      ! Step 2: compact each row in-place (merge duplicates, write new ia)
      ! Since each row is sorted, duplicates are adjacent.
      ! We write compacted data into the same arrays starting from position 1.
      ! ia(irow) gives the old start; we track wpos as the write position.
      wpos = 1
      do irow = 1, n
        rowcount(irow) = wpos   ! new start of this row (save for ia rebuild)
        j1 = ia(irow)
        do while (j1 < ia(irow+1))
          ja(wpos) = ja(j1)
          a_csr(wpos) = a_csr(j1)
          ! Merge any following entries with the same column index
          j2 = j1 + 1
          do while (j2 < ia(irow+1) .and. ja(j2) == ja(j1))
            a_csr(wpos) = a_csr(wpos) + a_csr(j2)
            j2 = j2 + 1
          end do
          wpos = wpos + 1
          j1 = j2
        end do
      end do
      ! Rebuild ia from the compacted row starts
      do irow = 1, n
        ia(irow) = rowcount(irow)
      end do
      ia(n+1) = wpos
    end block

    nzmod_out = ia(n+1) - 1   ! actual nnz after dedup
    if (dbg /= 0 .or. solve_count == 0) then
      write(*,'(a,i0,a,i0,a,i0)') 'PARDISO_WRAP: n=',n,' nzmod=',nzmod_out, &
           ' ia(n+1)=',ia(n+1)
    end if

    ! Diagnostic: verify every COO (row,col) exists in the CSR (debug only).
    if (dbg /= 0) then
    block
      integer :: k, irow2, jlo, jhi, jmid, ncol, nmiss
      nmiss = 0
      do k = 1, nzmod
        if (irn_coo(k) < 1 .or. irn_coo(k) > n) cycle
        irow2 = irn_coo(k)
        jlo = ia(irow2); jhi = ia(irow2+1)-1
        ncol = icn_coo(k)
        ! binary search for ncol in ja(jlo..jhi)
        do while (jlo < jhi)
          jmid = (jlo+jhi)/2
          if (ja(jmid) < ncol) then; jlo = jmid+1; else; jhi = jmid; end if
        end do
        if (jlo > ia(irow2+1)-1 .or. ja(jlo) /= ncol) nmiss = nmiss + 1
      end do
      if (nmiss > 0) then
        write(*,'(a,i0)') 'PARDISO_COO_CSR_MISMATCH: missing entries=',nmiss
      end if
    end block
    end if

    ! Initialize PARDISO on first call
    if (initialized == 0) then
      do i = 1, 64
        pt(i)%DUMMY = 0
      end do
      do i = 1, 64
        iparm(i) = 0
      end do
      iparm(1) = 1     ! explicit settings (MKL defaults fail on this Jacobian)
      iparm(2) = 2     ! serial nested dissection ordering
      iparm(4) = 0
      iparm(5) = 0
      iparm(8) = 10    ! max iterative refinement steps
      iparm(10) = 8    ! pivot perturbation 1e-8
      iparm(11) = 1    ! scaling vectors (required for unsymmetric)
      iparm(13) = 1    ! matching / weighted matching (required for unsymmetric)
      iparm(18) = -1
      iparm(19) = -1
      iparm(21) = 1    ! 1x1+2x2 pivoting (robust for indefinite/unsymmetric)
      iparm(24) = 0
      iparm(27) = 0    ! no matrix check
      initialized = 1
      prev_n = 0
      prev_nzmod = 0
      prev_phase_done = 0
    end if

    b_work = b
    maxfct = 1
    mnum = 1
    mtype = 11   ! real unsymmetric
    nrhs = 1
    msglvl = 0
    error = 0

    ! Decide whether we can reuse the existing factorization (phase 33 only).
    ! Within one timestep the corrector solves (modified Newton) all use the
    ! SAME matrix AJAC=1-EL(1)*H*J (EVALJA only rebuilds it when the Jacobian
    ! is refreshed), so refactorizing every solve wastes the dominant cost.
    ! We reuse ONLY when the CSR matrix is bit-identical to the one we
    ! factorized (same ia, ja AND a). A changed matrix -> fresh phase 12.
    ! NOTE: we deliberately do NOT reuse phase 22 for merely-same-pattern
    ! matrices — a stale symbolic permutation on this ill-conditioned
    ! Jacobian produced garbage solves (KFLAG=-2) in earlier testing.
    block
      logical :: same
      integer :: kk
      same = (prev_phase_done /= 0 .and. nzmod_out == prev_nzmod_out)
      if (same) then
        do kk = 1, n+1
          if (ia(kk) /= fac_ia(kk)) then; same = .false.; exit; end if
        end do
      end if
      if (same) then
        do kk = 1, nzmod_out
          if (ja(kk) /= fac_ja(kk) .or. a_csr(kk) /= fac_a(kk)) then
            same = .false.; exit
          end if
        end do
      end if

      if (.not. same) then
        ! Fresh factorization (phase 12): release old, analyze+factorize, store.
        if (prev_phase_done /= 0) then
          phase = -1
          call pardiso(pt, maxfct, mnum, mtype, phase, n, a_csr, ia, ja, &
                       perm, nrhs, iparm, msglvl, b_work, x, error)
          do i = 1, 64
            pt(i)%DUMMY = 0
          end do
        end if
        phase = 12
        call pardiso(pt, maxfct, mnum, mtype, phase, n, a_csr, ia, ja, &
                     perm, nrhs, iparm, msglvl, b_work, x, error)
        if (error /= 0) then
          write(*,*) 'PARDISO phase 12 error=',error
          if (error == -4 .or. error == -7 .or. error == -8 .or. &
              error == -9) then
            ifail = 14
          else
            ifail = 33
          end if
          prev_phase_done = 0
          return
        end if
        ! Store the factorized matrix for future identity comparison.
        if (.not. allocated(fac_ia)) allocate(fac_ia(NMAX+1))
        if (.not. allocated(fac_ja)) then
          allocate(fac_ja(NZMAX)); allocate(fac_a(NZMAX))
        end if
        do kk = 1, n+1
          fac_ia(kk) = ia(kk)
        end do
        do kk = 1, nzmod_out
          fac_ja(kk) = ja(kk)
          fac_a(kk) = a_csr(kk)
        end do
        prev_nzmod_out = nzmod_out
        prev_phase_done = 1
        n_fac = n_fac + 1
      else
        n_reuse = n_reuse + 1
        if (dbg /= 0 .and. mod(n_reuse,100) == 0) then
          write(*,'(a,i0,a,i0)') 'PARDISO_REUSE: reuse=',n_reuse,' fac=',n_fac
        end if
      end if
    end block
    prev_n = n
    prev_nzmod = nzmod
    prev_phase_done = 1

    ! Phase 33: solve
    phase = 33
    call pardiso(pt, maxfct, mnum, mtype, phase, n, a_csr, ia, ja, &
                 perm, nrhs, iparm, msglvl, b_work, x, error)
    if (error /= 0) then
      write(*,*) 'PARDISO phase 33 error=',error
      ifail = 33
      return
    end if

    ! Iterative refinement in quadruple precision + residual-based retry.
    ! Y12M refines with real(16) residuals; if the residual cannot be
    ! reduced it signals failure (ifail=33) so the integrator retries with
    ! a smaller step (as H->0, AJAC->I and the matrix becomes well-posed).
    ! We mirror that: solve A*x=b, compute r=b-A*x in real(16); if the
    ! residual is too large, solve A*dx=r and update x; repeat a bounded
    ! number of times. If still unconverged -> ifail=33 (clean retry)
    ! rather than letting a garbage x pollute the Newton state.
    block
      integer :: k, irow3, jj, iref
      integer, parameter :: MAXREF = 8
      double precision :: bnorm, rmax_coo, rmax_csr
      double precision :: resid_prev, dxnorm, xnorm
      double precision, allocatable :: rcoo(:), rcsr(:), dx(:)
      real(16), allocatable :: r16(:)
      allocate(rcoo(n)); allocate(rcsr(n)); allocate(dx(n)); allocate(r16(n))
      bnorm = 0.0d0
      do k = 1, n
        bnorm = max(bnorm, abs(b(k)))
        rcoo(k) = b(k); rcsr(k) = b(k)
      end do
      do k = 1, nzmod
        if (irn_coo(k) >= 1 .and. irn_coo(k) <= n .and. &
            icn_coo(k) >= 1 .and. icn_coo(k) <= n) then
          rcoo(irn_coo(k)) = rcoo(irn_coo(k)) - a_coo(k)*x(icn_coo(k))
        end if
      end do
      do irow3 = 1, n
        do jj = ia(irow3), ia(irow3+1)-1
          if (ja(jj) >= 1 .and. ja(jj) <= n) then
            rcsr(irow3) = rcsr(irow3) - a_csr(jj)*x(ja(jj))
          end if
        end do
      end do
      rmax_coo = 0.0d0; rmax_csr = 0.0d0
      do k = 1, n
        rmax_coo = max(rmax_coo, abs(rcoo(k)))
        rmax_csr = max(rmax_csr, abs(rcsr(k)))
      end do
      resid_coo = rmax_coo / max(bnorm,1.0d-300)
      resid_csr = rmax_csr / max(bnorm,1.0d-300)

      ! Early dump (debug): capture a bad matrix before refinement/ifail return.
      if (dbg_dump_res > 0.0d0 .and. resid_coo > dbg_dump_res) then
      block
        character(len=256) :: fn2
        write(fn2,'(a,i0,a)') 'pardiso_badmat_', solve_count+1, '.dat'
        open(unit=92, file=trim(fn2), status='replace', form='unformatted')
        write(92) n, nzmod
        write(92) irn_coo(1:nzmod); write(92) icn_coo(1:nzmod)
        write(92) a_coo(1:nzmod); write(92) b(1:n); write(92) x(1:n)
        close(92)
        write(*,'(a,1p,e10.2)') 'PARDISO_DUMP_BADMAT resid=', resid_coo
      end block
      end if

      ! Refine only if residual is above a modest threshold
      resid_prev = resid_coo
      if (resid_coo > 1.0d-10) then
        do iref = 1, MAXREF
          ! r16 = b - A*x in quad precision (COO is the true Jacobian)
          do k = 1, n
            r16(k) = real(b(k),16)
          end do
          do k = 1, nzmod
            if (irn_coo(k) >= 1 .and. irn_coo(k) <= n .and. &
                icn_coo(k) >= 1 .and. icn_coo(k) <= n) then
              r16(irn_coo(k)) = r16(irn_coo(k)) - &
                real(a_coo(k),16)*real(x(icn_coo(k)),16)
            end if
          end do
          do k = 1, n
            b_work(k) = real(r16(k),8)
          end do
          ! solve A*dx = r  (reuses existing factorization, phase 33)
          phase = 33
          call pardiso(pt, maxfct, mnum, mtype, phase, n, a_csr, ia, ja, &
                       perm, nrhs, iparm, msglvl, b_work, dx, error)
          if (error /= 0) exit
          dxnorm = 0.0d0; xnorm = 0.0d0
          do k = 1, n
            dxnorm = max(dxnorm, abs(dx(k)))
            x(k) = x(k) + dx(k)
            xnorm = max(xnorm, abs(x(k)))
          end do
          ! recompute residual (double) to test convergence
          do k = 1, n
            rcoo(k) = b(k)
          end do
          do k = 1, nzmod
            if (irn_coo(k) >= 1 .and. irn_coo(k) <= n .and. &
                icn_coo(k) >= 1 .and. icn_coo(k) <= n) then
              rcoo(irn_coo(k)) = rcoo(irn_coo(k)) - a_coo(k)*x(icn_coo(k))
            end if
          end do
          rmax_coo = 0.0d0
          do k = 1, n
            rmax_coo = max(rmax_coo, abs(rcoo(k)))
          end do
          resid_coo = rmax_coo / max(bnorm,1.0d-300)
          ! converged, or correction no longer helping / diverging
          if (resid_coo < 1.0d-12) exit
          if (dxnorm / max(xnorm,1.0d-300) < 1.0d-14) exit
          if (resid_coo > resid_prev * 10.0d0) exit  ! diverging; bail
          resid_prev = resid_coo
        end do
      end if

      if (resid_coo > 1.0d-8) then
        write(*,'(a,1p,e10.2,a,i0)') &
          'PARDISO_R: resid=', resid_coo, ' solve#', solve_count+1
      end if
      deallocate(rcoo); deallocate(rcsr); deallocate(dx); deallocate(r16)

      ! If residual still unacceptable -> tell integrator to retry (smaller H)
      if (resid_coo > 1.0d-7 .and. resid_coo == resid_coo) then
        ifail = 33
        return
      end if
      if (resid_coo /= resid_coo) then   ! NaN residual
        ifail = 33
        return
      end if
    end block

    ! Dump matrix+RHS+solution on a chosen solve (debug only).
    ! PARDISO_DUMP_AT=<1-based solve index>; PARDISO_DUMP_RES=<resid threshold>.
    solve_count = solve_count + 1
    if (dbg_dump_at > 0 .or. dbg_dump_res > 0.0d0) then
      if (solve_count == dbg_dump_at .or. &
          (dbg_dump_res > 0.0d0 .and. resid_coo > dbg_dump_res)) then
      block
        character(len=256) :: fname
        write(fname,'(a,i0,a)') 'pardiso_dump_', solve_count, '.dat'
        open(unit=91, file=trim(fname), status='replace', form='unformatted')
        write(91) n, nzmod
        write(91) irn_coo(1:nzmod)
        write(91) icn_coo(1:nzmod)
        write(91) a_coo(1:nzmod)
        write(91) b(1:n)
        write(91) x(1:n)
        close(91)
        write(*,'(a,i0,a,1p,e10.2)') 'PARDISO_DUMP: wrote solve ', &
              solve_count, ' resid=', resid_coo
      end block
      end if
    end if

    ifail = 0
  end subroutine pardiso_solve

  subroutine pardiso_cleanup_wrap()
    integer :: maxfct, mnum, mtype, phase, n, nrhs, msglvl, error
    integer :: iparm_loc(64), perm_loc(1)
    integer :: ia_loc(2), ja_loc(1)
    double precision :: a_loc(1), b_loc(1), x_loc(1)
    integer :: i

    if (initialized == 0) return

    phase = -1
    maxfct = 1
    mnum = 1
    mtype = 11
    n = 1
    nrhs = 1
    msglvl = 0
    error = 0
    do i = 1, 64
      iparm_loc(i) = 0
    end do
    ia_loc(1) = 1
    ia_loc(2) = 1

    call pardiso(pt, maxfct, mnum, mtype, phase, n, a_loc, ia_loc, ja_loc, &
                 perm, nrhs, iparm_loc, msglvl, b_loc, x_loc, error)
    do i = 1, 64
      pt(i)%DUMMY = 0
    end do
    initialized = 0
  end subroutine pardiso_cleanup_wrap

end module pardiso_wrap_mod
