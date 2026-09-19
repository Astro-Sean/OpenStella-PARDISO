      program ttdrv
      character*200 runname
      CHARACTER*80 Model,Sumprf,Sumcur,Depfile,Flxfile
      COMMON/Files/Model,Sumprf,Sumcur,Depfile,Flxfile
      integer l
      call getarg(1, runname)
      l = len_trim(runname)
      Model='../eve/'//runname(1:l)//'.mod'
      Sumprf=runname(1:l)//'.prf'
      Sumcur=runname(1:l)//'.crv'
      Depfile=runname(1:l)//'.dep'
      Flxfile=runname(1:l)//'.flx'
      write(*,*) 'ttdrv: runname=', runname(1:l)
      call tt4strad(runname(1:l))
      stop
      end
