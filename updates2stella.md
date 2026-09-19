# Updates to STELLA: PARDISO solver version

Notes on what changed in this fork and why. The fork replaces the Y12M sparse
direct solve with Intel MKL PARDISO inside the existing stiff BDF integrator.
Everything else (physics, hydro, opacity, radiation, time integrator) is
unchanged.

## Files changed

| File | Change |
|------|--------|
| `strad/pardiso_wrap.f90` | New. COO to CSR assembly, MKL PARDISO driver, quadruple-precision iterative refinement, recoverable `ifail=33`, factorization reuse on identical matrices, diagnostics gated by env vars. |
| `strad/stiffbghY12m.f` | One call site: `y12mff(...)` replaced by `pardiso_solve(N,ajac,irn,icn,nzmod,RHS,XSAVE,ifail)`, plus `USE pardiso_wrap_mod`. Nothing else differs. |
| `obj/Stellaf90GF.mak` | Adds `pardiso_wrap.f` to the object list, `-I$(MKLROOT)/include`, links `-lmkl_rt -lpthread -lm -ldl`. |
| `src/zone.inc` | `Mzon`/`NZ` templated by the build script (`sed Mzon=$NZON`). |
| `scripts/stella_job_template.slurm` | Exports `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `MKL_THREADING_LAYER=intel` so PARDISO threading follows the SLURM `cpus_per_task` allocation. |
| `build_parallel.sh` | New. Self-contained build script in the repo root. |

## Solver interface

`stiffbghY12m.f` calls, once per corrector iteration:

```fortran
call pardiso_solve(N, ajac, irn, icn, nzmod, FSAVE(NYDIM+1), XSAVE, ifail)
```

This solves `AJAC * X = RHS` where `AJAC = 1 - EL(1)*H*J` is the Newton/corrector
matrix of the BDF method. `ifail` follows Y12M semantics: `0` success, `33`
recoverable (integrator retries with a smaller step), `14` fatal/structural.

## COO to CSR conversion

STELLA passes the Jacobian in COO form (`irn`, `icn`, `ajac`). PARDISO wants CSR
with sorted, unique column indices per row. The wrapper:

1. counts entries per row to build `ia`;
2. scatters `(icn, ajac)` into `ja`/`a_csr`;
3. sorts each row by column;
4. merges duplicate `(row, col)` entries by summation. Y12M does this internally
   and PARDISO does not, so this is needed for correctness.

`nzmod_out = ia(n+1) - 1` is the deduplicated nonzero count.

## PARDISO settings

MKL defaults give bad LU factors on STELLA's ill-conditioned unsymmetric
Jacobian. The wrapper sets:

| iparm | value | reason |
|-------|-------|--------|
| 1 | 1 | use explicit settings (defaults failed) |
| 2 | 2 | nested-dissection fill-reducing ordering |
| 8 | 10 | max internal iterative-refinement steps |
| 10 | 8 | pivot perturbation 1e-8 |
| 11 | 1 | scaling vectors, needed for unsymmetric |
| 13 | 1 | weighted matching, needed for unsymmetric |
| 21 | 1 | 1x1 + 2x2 Bunch-Kaufman pivoting |
| 27 | 0 | no internal matrix check |

On the same dumped matrices these settings took residuals from about 865 down to
about 1e-14.

## Correctness checks

- Quadruple-precision refinement (`real(16)`): computes `r = b - A*x` in quad
  precision and applies a bounded correction loop, like Y12M does. This is what
  makes it survive this system, where double-precision refinement alone
  diverged.
- Residual validation: every solve checks `||b - A*x||/||b||` against the
  original COO matrix, and against the CSR as a conversion check.
- `ifail = 33` retry: if the residual cannot be reduced, or is NaN, the wrapper
  returns `33` so STELLA retries with a smaller step instead of feeding a bad
  `x` into the Newton state. Same semantics as Y12M's recoverable failure.

## Factorization strategy

Two bugs were found by dumping matrices mid-run and solving them offline.

**Bug 1: weak default pivoting.** Default `iparm` gave residual ~865 on a matrix
that scipy solved to ~1e-14. Fixed by the settings above.

**Bug 2: stale factorization reuse.** Reusing PARDISO phase-22 (numeric
refactorization) whenever `(n, nzmod)` matched recycled a stale symbolic
permutation. The Jacobian's nonzero pattern and conditioning change between
Newton steps even at constant `nzmod`, which gave residuals up to ~1e29 and
`KFLAG=-2`. Fixed by doing a fresh phase-12 factorization whenever the matrix
changes.

Current strategy: each solve compares the freshly built CSR against the
last-factorized one (`fac_ia`, `fac_ja`, `fac_a`).

- identical (`ia`, `ja`, `a` all equal): phase-33 only, triangular solve on the
  stored factors. This is the common corrector-iteration case, where `AJAC` is
  held fixed across modified-Newton iterations within one step.
- any difference: release and do a fresh phase-12, then store the new matrix.

Phase-22 is deliberately not used for same-pattern/different-values matrices;
that is exactly the case that produced the bad solves. Measured in a test run:
`reuse=200, fac=98`, so about two thirds of solves skip refactorization, with
zero KFLAG failures.

## Diagnostics

Read once at the first call into cached module variables, so there are no
per-solve env lookups.

| Env var | Effect |
|---------|--------|
| `PARDISO_DEBUG=1` | per-solve `PARDISO_WRAP:` size line, COO/CSR consistency check, `PARDISO_REUSE` counter every 100 reuses |
| `PARDISO_DUMP_AT=<k>` | dump matrix, RHS and solution of solve `k` to `pardiso_dump_<k>.dat` |
| `PARDISO_DUMP_RES=<tol>` | dump the first solve whose residual exceeds `tol` (`pardiso_badmat_*.dat`) |

Production runs should leave all of these unset.

## Build

```bash
bash build_parallel.sh <NZON>
```

`NZON` sets `Mzon`, the array bound. It must match the model's zone count.
Passing no argument defaults to 200 and prints a warning, which would silently
under-dimension a larger model. Build the code with the actual number of zones
(for the 500-zone benchmark: `bash build_parallel.sh 500`).

Output: `bin/xstella.exe` (ifort + MKL `mkl_rt` + OpenMP).

## Benchmark

500-zone model, `small` partition. All four runs reached about +8 d observer
time with zero KFLAG failures.

| Run | CPUs | Elapsed | Speedup |
|-----|------|---------|---------|
| Y12M baseline | 1 | 164.5 min | 1.0x |
| PARDISO | 1 | 78.3 min | 2.10x |
| PARDISO | 2 | 71.7 min | 2.29x |
| PARDISO | 3 | 69.4 min | 2.37x |

Physical equivalence (matched-time, 1254 common epochs): median per-zone
density/velocity difference about 0; max ~2-4% in isolated shock zones.
Bolometric lightcurve overlays from -20 to +8 d; final luminosity about
5.07e41 erg/s on all runs.

The reuse optimization landed after these timings, so current builds should be
at least this fast.

The main gain is a faster solver rather than strong parallelism: 2.1x already at
1 CPU, and threading only adds ~10-13% going from 1 to 3 cores.

## Caveats

- The stall at strong shocks is a time-integration problem, not a solver
  problem. This work lowers the cost per stiff micro-step; it does not reduce
  how many steps the integrator takes. Real improvements there need integrator
  changes (IMEX splitting, Rosenbrock-type methods, a looser corrector
  iteration limit, better step-size control).
- Keep `PARDISO_DEBUG` and `PARDISO_DUMP_*` unset for production.
