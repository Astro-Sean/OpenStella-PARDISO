      SubroutineStradIO(mode,Lunit,NFILE)
      IMPLICITREAL*8(A-H,O-Z)
      character*80NFILE,mode*2,Status*3
      integerLunit,IOST
      LOGICALLEXIST,DUMMY
      PARAMETER(NVARS=3)
      PARAMETER(NFREQ=96+4)
      PARAMETER(Mzon=500)
      PARAMETER(NYDIM=(NVARS+2*NFREQ)*Mzon,MAXDER=4)
      Parameter(Is=5)
      PARAMETER(NZ=3000000)
      Parameter(Nstage=28,Natom=15)
      PARAMETER(KOMAX=80)
      LogicalLSYSTEM
      Parameter(LSystem=.FALSE.)
      Parameter(Pi=3.1415926535897932d+00,hPlanc=1.0545716280D-27,Cs=2.9
     *979245800D+10,Boltzk=1.3806504000D-16,Avogar=6.0221417900D+23,AMbr
     *un=1.6605387832D-24,AMelec=9.1093821500D-28,echarg=4.8032042700D-1
     *0,CG=6.6742800000D-08,CMS=1.9884000000D+33,RSol=6.9551000000D+10,U
     *LGR=1.4000000000D+01,UPURS=1.0000000000D+07,ULGPU=7.0000000000D+00
     *,ULGEU=1.3000000000D+01,UPC=3.0856776000D+18,UTP=1.0000000000D+05,
     *URHO=1.0000000000D-06,CARAD=7.5657680191D-15,CSIGM=5.6704004778D-0
     *5,ERGEV=1.6021764864D-12,GRADeV=1.1604505285D+04,RADC=7.5657680191
     *D-02,CTOMP=4.0062048575D-01,CCAPS=2.6901213726D+01,CCAPZ=9.8964034
     *725D+00)
      IntegerZn(Natom),ZnCo(Natom+1)
      DimensionAZ(Natom)
      Common/AZZn/AZ,Zn,ZnCo
      Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap
      LOGICALFRST
      Parameter(Mfreq=130)
      Common/Kmzon/km,kmhap,Jac,FRST
      COMMON/STCOM1/t,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART
      COMMON/YMAX/YMAX(NYDIM)
      COMMON/YSTIF/Y(NYDIM,MAXDER+1)
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/HNT/HNT(7)
      PARAMETER(DELTA=1.d-05)
      PARAMETER(LICN=4*NZ,LIRN=2*NZ)
      LogicalNEEDBR
      COMMON/STJAC/THRMAT,HL,AJAC(LICN),IRN(LIRN),ICN(LICN),WJAC(NYDIM),
     *FSAVE(NYDIM*2),IKEEP(5*NYDIM),IW(8*NYDIM),IDISP(11),NZMOD,NEEDBR
      LOGICALCONV,CHNCND,SCAT,SEP
      COMMON/CUTOFF/FLOOR(NVARS+1),Wacc(NVARS+1),FitTau,TauTol,Rvis,CONV
     *,CHNCND,SCAT,SEP
      LogicalLTHICK
      COMMON/THICK/LTHICK(Nfreq,Mzon)
      COMMON/CONVEC/UC(Mzon),YAINV(Mzon)
      COMMON/RAD/EDDJ(Mzon,Nfreq),EDDH(Mzon),HEDD(Nfreq),HEDRAD,CLIGHT,C
     *KRAD,UFREQ,CFLUX,CCL,CLUM,CLUMF,CIMP,NTHICK(NFREQ),NTHNEW(NFREQ),N
     *CND,KRAD,NFRUS
      LOGICALEDTM
      COMMON/RADOLD/HEDOLD,HINEDD,EDTM
      Common/newedd/EddN(Mzon,Nfreq),HEdN(Nfreq),tfeau
      Common/oldedd/EddO(Mzon,Nfreq),HEdo(Nfreq),trlx
      Common/cnlast/Cnlast
      Common/Dhap/DHaphR(Mzon,Nfreq)
      COMMON/BAND/FREQ(NFREQ+1),FREQMN(NFREQ),WEIGHT(130),HAPPAL(NFREQ),
     *HAPABSRON(NFREQ),HAPABS(NFREQ),DLOGNU(NFREQ)
      PARAMETER(NFRMIN=Nfreq/2)
      IntegerdLfrMax
      Common/observer/wH(Mfreq),cH(Mfreq),zerfr,Hcom(Mfreq),Hobs(Mfreq),
     *freqob(Mfreq),dLfrMax
      Parameter(NP=15+15-1)
      Common/famu/fstatic(0:NP+1,Nfreq),fobs_corr(0:NP+1,Mfreq),fcom(0:N
     *P+1,Mfreq),amustatic(0:NP+1)
      Common/rays/Pray(0:Np+1),fout(0:NP+1,Mfreq),abermu(0:NP+1),NmuNzon
      COMMON/LIM/Uplim,Haplim
      COMMON/AMM/DMIN,DM(Mzon),DMOUT,AMINI,AM(Mzon),AMOUT
      COMMON/Centr/RCE,Nzon
      Common/InEn/AMHT,EBurst,tBurst,tbeght
      COMMON/RADPUM/AMNI,XMNi,XNi,KmNick
      COMMON/RADGAM/FJgam(Mzon,2),toldg,tnewg
      COMMON/RADGlg/FJglog(Mzon,2)
      COMMON/CHEM/CHEM0(Mzon),RTphi(0:Mzon+1),EpsUq
      COMMON/REGIME/NREG(Mzon)
      doubleprecisionNRT
      COMMON/AQ/AQ,BQ,DRT,NRT
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/QNRGYE/QNUC,RGASA,YELECT
      COMMON/CKN1/CK1,CK2,CFR,CRAP,CRAOLD
      LOGICALEVALJA,OLDJAC,BADSTE
      COMMON/EVAL/EVALJA,BADSTE,OLDJAC
      LogicalRadP
      COMMON/RadP/RadP
      COMMON/ARG/TP,PL,CHEM,LST,KENTR,JURS
      COMMON/RESULT/P,Egas,Sgas,ENG,CAPPA,PT,ET,ST,ENGT,CAPT,NZR
      COMMON/ABUND/XYZA,Yat(Natom)
      COMMON/AZ/AS,ZS,SCN
      COMMON/STR/PPL,EPL,SPL,ENGPL,CAPPL,CP,GAM,DA,DPE,DSE,DSP,BETgas
      COMMON/XELECT/XE,XET,XEPL,PE,Ycomp
      COMMON/URScap/Tpsqrt,Psicap,Scap,ScapT,ScapPl,ZMean,YZMean,ZMT,ZMP
     *l,YZMT,YZMPl
      COMMON/BURNCC/CC,CCTP,CCPL,YDOT
      COMMON/ABarr/YABUN(Natom,Mzon)
      COMMON/UNSTL/UL,UPRESS,UE,UEPS,UCAP,UTIME,UV,UFLUX,UP
      COMMON/TAIL/KTAIL
      COMMON/UNINV/UPI,UEI
      COMMON/UNBSTL/UR,UM,UEPRI,ULGP,ULGE,ULGV,ULGTM,ULGEST,ULGFL,ULGCAP
     *,ULGEPS
      COMMON/CONUR/EIT,DST,BBRCNR(5)
      COMMON/BAL/EL(MAXDER+1),YENTOT(MAXDER+1),ETOT0,ELVOL,ELSURF,ELTOT,
     *TPSURF,HOLDBL,ELOST,EKO,RADBEG
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,NOUT
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
      REAL*8TPMAX(MAXDER+1),TQ(4)
      COMMON/TAU/TAU(Mzon+1),FLUX(Mzon)
      common/tauubvri/tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(M
     *zon)
      COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH
      PARAMETER(NFUNC=6)
      REAL*4WORK(Mzon+2,NFREQ),WRK(Mzon,4)
      REAL*8WRKX(Mzon),WORKX(Mzon+2)
      COMMON/STEPD/WRKX,WORKX,TPHOT,TEFF,WORK,WRK,NPHOT,NZM
      PARAMETER(TMCRIT=1.D-6,TPNSE=5.D0,EPGROW=0.02D0)
      Common/RUTP/Ry(Mzon),Uy(Mzon),Ty(Mzon),Press(Mzon),Rho(Mzon)
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      Parameter(Lcurdm=200000)
      RealTcurv
      IntegerNFRUSED
      REAL*8Flsave
      Common/Curve/tcurv(8,Lcurdm),Depos(Lcurdm),Flsave(MFREQ+1,Lcurdm),
     *NFRUSED(Lcurdm),Lsaved
      LOGICALBEGRUN
      Common/BEGR/BEGRUN
      CHARACTER*80Model,Sumprf,Sumcur,Depfile,Flxfile
      COMMON/Files/Model,Sumprf,Sumcur,Depfile,Flxfile
      CHARACTER*1app
      LogicalGivdtl
      Common/ABGrap/NSTA,NSTB,TcurA,TcurB,Givdtl
      REAL*8MBOL,MU,MB,MV,MR,MI,MBOL1
      COMMON/COLOR/MBOL,MU,MB,MV,MR,MI,UMB,BMV,MBOL1,LubvU,LubvB,LubvV,L
     *ubvR,LubvI,Lyman
      COMMON/DETAIL/QRTarr(Mzon),UUarr(Mzon),ArrLum(Mzon),Acc(Mzon)
      Common/XYZ/XA,YA,URM
      COMMON/STSAVE/RMAX,TREND,OLDL0,RC,HOLD,EDN,E,EUP,BND,EPSOLD,TOLD,M
     *EO,NOLD,NQ,LNQ,IDOUB
      IF(.NOT.(mode=='rm'.OR.mode=='RM'))GOTO09999
      write(*,*)' ------','-->Entering Node %_readmodel:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      WRITE(*,*)' StradIO opened file ',NFILE
09988 CONTINUE
      READ(Lunit,IOSTAT=IOST,END=10)t,NSTEP,KTAIL,NZON,NQUSED,KRAD,NCND,
     *NFRUS,NREG,XA,YA,XYZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,ULGP,ULGV
     *,ULGE,UEPRI,UTIME,UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,UM,UV,URM,DM,AM
     *,DMIN,AMINI,DMOUT,AMOUT,YABUN,(Y(ISAVE,1),ISAVE=1,NZON*NVARS+2*KRA
     *D)
      WRITE(*,*)' StradIO READ STEP=',NSTEP
      IF(.NOT.(IOST/=0.or.(IABS(NSTB)>=IABS(NSTA).AND.Nstep>=IABS(NSTA))
     *))GOTO09988
10    CONTINUE
      if(Nzon==0)then
      WRITE(*,*)' StradIO READ STEP=',NSTEP
      write(*,*)'  Error: in stradio.trf %_readmodel: Nzon=',Nzon
      stop111
      endif
      write(*,*)' ------','<--Leaving  Node %_readmodel:'
      GOTO09989
09999 CONTINUE
      IF(.NOT.(mode=='cm'.OR.mode=='CM'))GOTO09998
      write(*,*)' ------','-->Entering Node %_contmodel:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      WRITE(*,*)' StradIO opened file ',NFILE
09985 CONTINUE
      READ(Lunit,IOSTAT=IOST,END=11,ERR=111)t,NSTEP,KTAIL,NZON,NQUSED,KR
     *AD,NCND,NFRUS,NREG,XA,YA,XYZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,U
     *LGP,ULGV,ULGE,UEPRI,UTIME,UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,UM,UV,U
     *RM,DM,AM,DMIN,AMINI,DMOUT,AMOUT,YABUN,(Y(ISAVE,1),ISAVE=1,NZON*NVA
     *RS+2*KRAD),JSTART,METH,KFLAG,EVALJA,OLDJAC,BADSTE,Nclast,NFUN,NJAC
     *,NITER,NFAIL,MEO,NOLD,LNQ,IDOUB,NTHICK,LTHICK,HUSED,RMAX,TREND,OLD
     *L0,RC,HOLD,EDN,E,EUP,BND,EPSOLD,TOLD,TAUOLD,AMHT,EBurst,TBurst,AMN
     *I,XMNI,YENTOT,H,HOLDBL,ETOT0,FREQ,FREQMN,DLOGNU,WEIGHT,HNT
      WRITE(*,*)' StradIO READ STEP=',NSTEP
      IF(.NOT.(IOST/=0.or.(IABS(NSTB)>=IABS(NSTA).AND.Nstep>=IABS(NSTA))
     *))GOTO09985
      goto11
111   write(*,*)'ReadError ''CONTREAD''(first read),IOSTAT=',IOST
11    CONTINUE
      write(*,*)' ------','<--Leaving  Node %_contmodel:'
      GOTO09989
09998 CONTINUE
      IF(.NOT.(mode=='rc'.OR.mode=='RC'))GOTO09997
      GOTO09989
09997 CONTINUE
      IF(.NOT.(mode=='wm'.OR.mode=='WM'))GOTO09996
      write(*,*)' ------','-->Entering Node %_writemodel:'
      write(*,*)' ent wm'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      CRAOLD=CRAP
      WRITE(Lunit)t,NSTEP,KTAIL,NZON,NQUSED,KRAD,NCND,NFRUS,NREG,XA,YA,X
     *YZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,ULGP,ULGV,ULGE,UEPRI,UTIME,
     *UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,UM,UV,URM,DM,AM,DMIN,AMINI,DMOUT,
     *AMOUT,YABUN,(Y(ISAVE,1),ISAVE=1,NZON*NVARS+2*KRAD),JSTART,METH,KFL
     *AG,EVALJA,OLDJAC,BADSTE,Nclast,NFUN,NJAC,NITER,NFAIL,MEO,NOLD,LNQ,
     *IDOUB,NTHICK,LTHICK,HUSED,RMAX,TREND,OLDL0,RC,HOLD,EDN,E,EUP,BND,E
     *PSOLD,TOLD,TAUOLD,AMHT,EBurst,TBurst,AMNI,XMNI,YENTOT,H,HOLDBL,ETO
     *T0,FREQ,FREQMN,DLOGNU,WEIGHT,HNT
      Begrun=.False.
      write(*,*)' leav wm'
      write(*,*)' ------','<--Leaving  Node %_writemodel:'
      GOTO09989
09996 CONTINUE
      IF(.NOT.(mode=='wc'.OR.mode=='WC'))GOTO09995
      write(*,*)' ------','-->Entering Node %_writecurve:'
      write(*,*)' ent wc'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      write(*,*)' opened wc, Lsaved=',Lsaved
      WRITE(Lunit)Lsaved,((Tcurv(ISAVE,JSAVE),ISAVE=1,8),JSAVE=1,Lsaved)
      NFILE=Depfile
      LUNIT=22
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      Write(Lunit)Lsaved,(Depos(JSAVE),JSAVE=1,Lsaved)
      write(*,*)' stradio Lsaved=',Lsaved
      NFILE=Flxfile
      LUNIT=23
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      Write(Lunit)Lsaved,(Tcurv(1,JSAVE),NFRUSED(JSAVE),(Flsave(l,JSAVE)
     *,l=1,Mfreq),JSAVE=1,Lsaved)
      Begrun=.False.
      Lsaved=0
      write(*,*)' leav wc'
      write(*,*)' ------','<--Leaving  Node %_writecurve:'
      GOTO09989
09995 CONTINUE
      IF(.NOT.(mode=='we'.OR.mode=='WE'))GOTO09994
      write(*,*)' ------','-->Entering Node %_writeeve:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      write(*,*)'%_writeeve: Nzon=',Nzon
      WRITE(Lunit)t,NSTEP,KTAIL,NZON,NQUSED,KRAD,NCND,NFRUS,NREG,XA,YA,X
     *YZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,ULGP,ULGV,ULGE,UEPRI,UTIME,
     *UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,UM,UV,URM,DM,AM,DMIN,AMINI,DMOUT,
     *AMOUT,YABUN,(Y(ISAVE,1),ISAVE=1,NZON*NVARS+2*KRAD)
      write(*,*)' ------','<--Leaving  Node %_writeeve:'
      GOTO09989
09994 CONTINUE
      IF(.NOT.(mode=='sm'.OR.mode=='SM'))GOTO09993
      write(*,*)' ------','-->Entering Node %_stinfomodel:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      READ(Lunit,IOSTAT=IOST,END=30)t,NSTEP,KTAIL,NZON,NQUSED,KRAD,NCND,
     *NFRUS,NREG,XA,YA,XYZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,ULGP,ULGV
     *,ULGE,UEPRI,UTIME,UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,UM,UV,URM,DM,AM
     *,DMIN,AMINI,DMOUT,AMOUT,YABUN,(Y(ISAVE,1),ISAVE=1,NZON*NVARS+2*KRA
     *D),JSTART,METH,KFLAG,EVALJA,OLDJAC,BADSTE,Nclast,NFUN,NJAC,NITER,N
     *FAIL,MEO,NOLD,LNQ,IDOUB,NTHICK,LTHICK,HUSED,RMAX,TREND,OLDL0,RC,HO
     *LD,EDN,E,EUP,BND,EPSOLD,TOLD,TAUOLD,AMHT,EBurst,TBurst,AMNI,XMNI,Y
     *ENTOT,H,HOLDBL,ETOT0,FREQ,FREQMN,DLOGNU,WEIGHT,HNT
      WRITE(*,*)' StradIO READ STEP=',NSTEP
      Return
30    LUNIT=-1
      write(*,*)' ------','<--Leaving  Node %_stinfomodel:'
      GOTO09989
09993 CONTINUE
      IF(.NOT.(mode=='sc'.OR.mode=='SC'))GOTO09992
      write(*,*)' ------','-->Entering Node %_stinfocurve:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      READ(Lunit,IOSTAT=IOST,END=40)Lsaved,((Tcurv(ISAVE,JSAVE),ISAVE=1,
     *8),JSAVE=1,Lsaved)
      write(*,*)' Curve read',lsaved
      return
40    write(*,*)' EOF reached in file crv',lsaved
      LUNIT=-1
      write(*,*)' ------','<--Leaving  Node %_stinfocurve:'
      GOTO09989
09992 CONTINUE
      IF(.NOT.(mode=='sd'.OR.mode=='SD'))GOTO09991
      write(*,*)' ------','-->Entering Node %_stinfodepos:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      READ(Lunit,IOSTAT=IOST,END=50)Lsaved,(Depos(JSAVE),JSAVE=1,Lsaved)
      write(*,*)' Depos read ',Lsaved
      Return
50    write(*,*)' EOF reached in file depos',lsaved
      LUNIT=-1
      write(*,*)' ------','<--Leaving  Node %_stinfodepos:'
      GOTO09989
09991 CONTINUE
      IF(.NOT.(mode=='sf'.OR.mode=='SF'))GOTO09990
      write(*,*)' ------','-->Entering Node %_stinfoflux:'
      INQUIRE(FILE=NFILE,EXIST=LEXIST)
      IF(LEXIST)THEN
      STATUS='OLD'
      ELSE
      STATUS='NEW'
      ENDIF
      OPEN(Lunit,IOSTAT=IOST,FILE=NFILE,STATUS='unknown',FORM='UNFORMATT
     *ED')
      READ(Lunit,IOSTAT=IOST,END=60)Lsaved,(Tcurv(1,JSAVE),NFRUSED(JSAVE
     *),(Flsave(l,JSAVE),l=1,Mfreq),JSAVE=1,Lsaved)
      write(*,*)' Flx read ',Lsaved
      Return
60    write(*,*)' EOF reached in file flx',lsaved
      LUNIT=-1
      write(*,*)' ------','<--Leaving  Node %_stinfoflux:'
      GOTO09989
09990 CONTINUE
      write(*,'('' Unknown mode:'',A)')mode
      stop16
09989 CONTINUE
      Return
      End
