      SUBROUTINEBEGIN
      IMPLICITREAL*8(A-H,O-Z)
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
      Character*80Opafile
      real*4hptabab,hptababron,hptabsc,hpsavsc
      Common/Opsave/hpsavsc(Nfreq,14,14,Mzon/(Mzon/2)+1,6)
      Common/OpBand/TpTab(14),RhoTab(14),STab(6),Wavel(Nfreq),YATab(Nato
     *m),hptabab(Nfreq,14,14,Mzon/(Mzon/2)+1,2),hptababron(Nfreq,14,14,M
     *zon/(Mzon/2)+1,2),hptabsc(Nfreq,14,14,Mzon/(Mzon/2)+1,2),Msta,Nrho
     *,Ntp,im
      real*4hpbanab1(Nfreq,14,14,Mzon/(Mzon/2)+1),hpbanabron1(Nfreq,14,1
     *4,Mzon/(Mzon/2)+1),hpbansc1(Nfreq,14,14,Mzon/(Mzon/2)+1),hpbanab2(
     *Nfreq,14,14,Mzon/(Mzon/2)+1),hpbanabron2(Nfreq,14,14,Mzon/(Mzon/2)
     *+1),hpbansc2(Nfreq,14,14,Mzon/(Mzon/2)+1)
      Equivalence(hptabab(1,1,1,1,1),hpbanab1(1,1,1,1)),(hptababron(1,1,
     *1,1,1),hpbanabron1(1,1,1,1)),(hptabsc(1,1,1,1,1),hpbansc1(1,1,1,1)
     *),(hptabab(1,1,1,1,2),hpbanab2(1,1,1,1)),(hptababron(1,1,1,1,2),hp
     *banabron2(1,1,1,1)),(hptabsc(1,1,1,1,2),hpbansc2(1,1,1,1))
      Common/tintrp/stmlog(6),tdlog,thaplog1,thaplog2,istold,Opafile
      Common/dumfreq/dumFreq(Nfreq+1),dumFreqmn(Nfreq),dumwavel(Nfreq)
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
      CHARACTER*80STRING
      Dimensionindfr(6),indop(6)
      doubleprecisionAMi1
      CHARACTER*3status,NFILE*80
      Logicalpetread
      Datapetread/.false./
      Dataindfr/1,2,2,3,3,3/,indop/1,1,2,3,3,3/
      DATAWLMAX/5.0D+04/,WLMIN/1.D+00/
      logical,save::EachStepOutput=.false.
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %R:'
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RD:'
      Nc=8
      READ(Nc,'(A)')
      Read(Nc,*)EPS,HMIN,HMAX
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)METHN,JSTART,MAXORD,KNadap
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)NSTA,NSTB,TcurA,TcurB
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)NSTMAX,NDebug,NOUT,IOUT,MBATCH
      write(*,*)' 1'
      Mbatch=Min(Mbatch,Lcurdm)
      READ(Nc,'(A)')
      Read(Nc,*)AMNI,XMNI
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)AMHT,EBurst,tBurst,tbeght
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)EKOB,AI1,AI2,AI3,US
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)THRMAT,CRAP,CONV,EDTM,CHNCND,Givdtl
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)FLOOR(1),FLOOR(2),FLOOR(3),FLOOR(4)
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)Wacc(1),Wacc(2),Wacc(3),Wacc(4)
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)FitTau,TauTol,Rvis,AQ,BQ,DRT,NRT,SCAT
      write(*,*)' 1'
      READ(Nc,'(A)')
      Read(Nc,*)NnTO
      if(NnTo<0)then
      EachStepOutput=.true.
      NnTo=abs(NnTo)
      endif
      doito=1,NnTO
      Read(Nc,*)TO(ito)
      enddo
      If(LSystem)then
      READ(5,*)IRC
      If(IRC==0)then
      BEGRUN=.FALSE.
      else
      BEGRUN=.TRUE.
      endif
      endif
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RD:'
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM:'
      IF(BEGRUN.AND.IABS(NSTB)==1)THEN
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_Start:'
      Lunit=9
      NFILE=Model
      write(*,*)'>>',Nfile,'<<'
      callStradIO('rm',Lunit,NFILE)
      write(*,*)' Eko, Ekob: ',Eko,Ekob
      EKO=EKO+EKOB
      If(EKO>0.)then
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_Start_Vel:
     *'
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_Start_Vel_
     *tri:'
      i=0
09997 CONTINUE
      i=i+1
      IF(.NOT.((am(i)-amini)/(amout-amini)>=ai1.or.i==nzon))GOTO09997
      i1=i
      write(*,'(a,i5,1p,e10.3)')' i1, ai1:',i1,ai1
09994 IF(.NOT.(i<nzon))GOTO09993
      i=i+1
      IF(.NOT.((am(i)-amini)/(amout-amini)>=ai2.or.i==nzon))GOTO09994
09993 CONTINUE
      i2=i
      write(*,'(a,i5,1p,e10.3)')' i2, ai2:',i2,ai2
09991 IF(.NOT.(i<nzon))GOTO09990
      i=i+1
      IF(.NOT.((am(i)-amini)/(amout-amini)>=ai3.or.i==nzon))GOTO09991
09990 CONTINUE
      i3=i
      write(*,'(a,i5,1p,e10.3)')' i3, ai3:',i3,ai3
      write(*,'(a,1p,3e15.8)')'amout-amini =',(amout-amini)*UM
      write(*,'(a,1p,3e15.8)')'       am(i):',am(i1)*UM,am(i2)*UM,am(i3)
     **UM
      write(*,'(a,1p,3e15.8)')' am(i)-amini:',(am(i1)-amini)*UM,(am(i2)-
     *amini)*UM,(am(i3)-amini)*UM
      write(*,'(a,1p,3e12.5)')'  dam(i)/dAM:',(am(i1)-amini)/(amout-amin
     *i),(am(i2)-amini)/(amout-amini),(am(i3)-amini)/(amout-amini)
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_Start_Vel_
     *tri:'
      if(I1<=1)then
      AMi1=ai1*(amout-amini)+AMini
      else
      AMi1=AM(I1)
      endif
      Z1=US/(AM(I2)-AMi1)
      IF(I2/=I3)Z2=US/(AM(I2)-AM(I3))
      doK=I1,I3
      If(K<=I2)then
      Uy(K)=Z1*(AM(K)-AMi1)
      else
      Uy(K)=Z2*(AM(K)-AM(I3))
      endif
      enddo
      Z1=0.D0
      K1=I1+1
      doK=K1,I3
      if(K<Nzon)then
      Z1=Z1+Uy(K)**2*(DM(K)+DM(K+1))
      else
      Z1=Z1+Uy(K)**2*(DM(K)+DMout)
      endif
      enddo
      doK=I1,I3
      Y(Nzon+K,1)=Uy(K)*2.D0*SQRT((EKO/UEPRI)/Z1)
      enddo
      write(*,'(a,1p,3e12.5)')'  uy(i):',Y(Nzon+I1,1),Y(Nzon+I2,1),Y(Nzo
     *n+I3,1)
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_Start_Vel:
     *'
      endif
      NCND=NZON
      NFRUS=NFREQ
      KRAD=(NZON-NCND)*NFRUS
      TAUOLD=0.D0
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_Start:'
      ELSE
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_Conr:'
      if(NSTA<0)then
      Lunit=10
      NFILE=Sumprf
      else
      Lunit=12
      NFILE='test.prf'
      endif
      callStradIO('cm',Lunit,NFILE)
      WRITE(*,*)' Begrad  READ STEP=',NSTEP
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_Conr:'
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_Curv:'
      if(NSTA<0)then
      Lunit=11
      NFILE=Sumcur
      else
      Lunit=13
      NFILE='test.crv'
      endif
      callStradIO('rc',Lunit,NFILE)
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_Curv:'
      ENDIF
      Z1=0.d0
      Z2=0.d0
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      If(abs(KNadap)==5)then
      If(K<NZON)Then
      DM2=DM(K+1)
      else
      DM2=DMOUT
      endif
      Z1=Z1+Uy(K)**2*(DM(K)+DM2)
      Z2=Z2+Ry(K)**2*(DM(K)+DM2)
      endif
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      If(abs(KNadap)==5)then
      URout=sqrt(Z1/Z2)
      doK=1,Nzon
      Uy(K)=URout*Ry(K)
      Y(Nzon+K,1)=Uy(K)
      enddo
      endif
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_nickel:'
      km=1
      print*,'AMINI=',Amini,UM
      dowhile((AM(km)-AMini)*UM<=AMNi*1.0000000001d0.AND.km<nzon)
      km=km+1
      enddo
      If(km>=nzon)then
      write(*,*)'         km, (AM(km)-AMini)*UM,       AMNi'
      print*,km,(AM(km)-AMini)*UM,AMNi*1.0000000001d0
      write(*,*)' Begrad: AMNI in error! km=',km,'>=nzon =',nzon
      stop
      Endif
      kmnick=km-1
      AMNi=(AM(max(kmnick,1))-AMini)*UM
      XNI=XMNI/AMNi
      If(XNI>=1.d0.AND.KNadap>0)then
      write(*,*)' Begrad: AMNI too small ! XNI=',XNI
      stop
      Endif
      If(KNadap>0)then
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_nickel_Fe:
     *'
      iferrum=1
09988 IF(.NOT.(iferrum<=Natom.AND.Zn(Iferrum)/=26))GOTO09987
      iferrum=iferrum+1
      GOTO09988
09987 CONTINUE
      If(iferrum>Natom)Then
      write(*,*)' Begrad: Ferrum not found !'
      stop
      Endif
      DO09985km=1,kmnick
      YABUN(iferrum,km)=XNI/AZ(iferrum)
      sum=0.d0
      DO09982j=1,Natom
      if(j/=iferrum)sum=sum+YABUN(j,km)*AZ(j)
09982 CONTINUE
      If(sum<1.d-5)Then
      write(*,*)' Begrad: Sum too small=',sum
      stop
      Endif
      DO09979j=1,Natom
      if(j/=iferrum)YABUN(j,km)=YABUN(j,km)*(1.d0-XNI)/sum
09979 CONTINUE
09985 CONTINUE
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_nickel_Fe:
     *'
      else
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RM_nickel_Nif
     *or:'
      write(*,*)'Read .XNI'
      read(28)(XNIfor(i),i=1,Nzon)
      Xmni=0.d0
      dokm=1,Nzon
      Xmni=Xmni+Xnifor(km)*dm(km)*UM
      enddo
      write(*,*)' Ni mass=',Xmni
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_nickel_Nif
     *or:'
      endif
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM_nickel:'
      IF(CONV)THEN
      doI=1,NZON
      Y(I+(NVARS-2)*NZON,1)=UC(I)
      Y(I+(NVARS-1)*NZON,1)=YAINV(I)
      enddo
      ENDIF
      NSTMAX=NSTEP+NSTMAX
      IF(MOD(NSTMAX,MBATCH)/=0)THEN
      NSTMAX=(NSTMAX/MBATCH+1)*MBATCH
      WRITE(4,'( '' NSTMAX CHANGED:'',I6)')NSTMAX
      ENDIF
      JSTART=0
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RM:'
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RU:'
      UR=10.D0**ULGR
      CLIGHT=CS*UTIME/UR
      IF(CRAOLD/=CRAP)THEN
      T=T*(CRAOLD/CRAP)
      H=H*(CRAOLD/CRAP)
      CK1=CK1*(CRAOLD/CRAP)
      CK2=CK2*(CRAOLD/CRAP)
      ENDIF
      UFREQ=BOLTZK*UTP/(2.d0*PI*HPLANC)
      CKRAD=6.D+01/PI**4*CSIGM*UTP**4*UTIME**3/(URHO*UR**3)
      CCL=CS*1.D+08/UFREQ
      CFLUX=60.D0*CSIGM*(UTP/PI)**4
      CLUM=32.D0*PI/3.D0*(CSIGM*UR*UTP**4/URHO)
      CLUMF=4.D0*PI*UR**2*CFLUX
      CIMP=CFLUX*UTIME**2/(CS*UR**2*URHO)
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RU:'
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %R_haptab:'
      If(ABS(Knadap)<=3)then
      Open(unit=2,file=Opafile,form='unformatted')
      Else
      donf=1,6
      nu=30+nf
      app=char(ichar('0')+nf)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=nu,file=Opafile,form='unformatted')
      write(*,*)' opened unit=',nu,' file=',Opafile
      enddo
      Open(unit=29,file=Opafile(1:Length(Opafile)-1)//'ab',form='unforma
     *tted')
      endif
      GOTO(09976,09975,09974),indop(abs(KNadap))
      GOTO09973
09976 CONTINUE
      GOTO09973
09975 CONTINUE
      If(.NOT.petread)then
      read(2)Nfreq0,Msta,Nrho,Ntp,YATab,(Wavel(iif),iif=1,96),TpTab,RhoT
     *ab,STab,EpsBand,EppBand
      petread=.true.
      endif
      GOTO09973
09974 CONTINUE
      DO09972ihp=1,6
      read(30+ihp)nw,Stime,Nfreq0,Msta,Nrho,NTp,dumWavel,dumFreq,dumFreq
     *mn,TpTab,RhoTab,hpbansc2
      DO09969im=1,Nzon/(Mzon/2)
      DO09966iro=1,Nrho
      DO09963itp=1,NTp
      DO09960L=1,Nfreq0
      hpsavsc(L,itp,iro,im,ihp)=hpbansc2(L,itp,iro,im)
09960 CONTINUE
09963 CONTINUE
09966 CONTINUE
09969 CONTINUE
      write(*,*)' read unit=',30+ihp
      close(30+ihp)
      If(ihp==nw.or.ABS(knadap)==6.or.ihp==1)then
      stmlog(ihp)=log(Stime)
      write(*,*)' ihp stmlog=',stmlog(ihp),ihp
      else
      write(*,*)' error reading haptab, ihp=',ihp,' nw=',nw
      stop324
      endif
09972 CONTINUE
      read(29)hpbanab2
      close(29)
09973 CONTINUE
      write(*,*)' Nfreq: ',Nfreq0,Nfreq,' Msta,Nrho,Ntp: ',Msta,Nrho,Ntp
      tdlog=log(max(t*Utime/8.6d+04,hmin))
      IF(.NOT.(tdlog<=stmlog(1)-(stmlog(2)-stmlog(1))/2.d0.or.tdlog>stml
     *og(6).or.abs(Knadap)==6))GOTO09957
      istim=0
      GOTO09955
09957 CONTINUE
      IF(.NOT.(tdlog>stmlog(1)-(stmlog(2)-stmlog(1))/2.d0.AND.tdlog<=stm
     *log(1)))GOTO09956
      istim=1
      GOTO09955
09956 CONTINUE
      istim=2
09954 IF(.NOT.(stmlog(istim)<tdlog))GOTO09953
      istim=istim+1
      GOTO09954
09953 CONTINUE
09955 CONTINUE
      write(*,*)' istim=',istim
      istold=istim
      GOTO(09951,09950),istim+1
      if(istim>1.and.istim<=6)then
      thaplog1=stmlog(istim-1)
      thaplog2=stmlog(istim)
      else
      write(*,*)' in begrad istim=',istim
      stop' wrong istim in begrad'
      endif
      GOTO09949
09951 CONTINUE
      DO09948im=1,Nzon/(Mzon/2)
      DO09945iro=1,14
      DO09942itp=1,14
      DO09939L=1,Nfreq
      hpbanab1(L,itp,iro,im)=hpbanab2(L,itp,iro,im)
      hpbansc1(L,itp,iro,im)=hpbansc2(L,itp,iro,im)
09939 CONTINUE
09942 CONTINUE
09945 CONTINUE
09948 CONTINUE
      thaplog1=stmlog(1)
      thaplog2=stmlog(6)
      GOTO09949
09950 CONTINUE
      DO09936im=1,Nzon/(Mzon/2)
      DO09933iro=1,14
      DO09930itp=1,14
      DO09927L=1,Nfreq
      hpbanab1(L,itp,iro,im)=hpbanab2(L,itp,iro,im)
      hpbansc1(L,itp,iro,im)=hpbansc2(L,itp,iro,im)
09927 CONTINUE
09930 CONTINUE
09933 CONTINUE
09936 CONTINUE
      thaplog1=stmlog(1)-(stmlog(2)-stmlog(1))/2.d0
      thaplog2=stmlog(1)
09949 CONTINUE
      If(istim/=0)then
      app=char(ichar('0')+istim)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=30,file=Opafile,form='unformatted')
      read(30)nw,Stime,Nfreq0,Msta,Nrho,NTp,dumWavel,dumFreq,dumFreqmn,T
     *pTab,RhoTab,hpbansc2
      close(30)
      endif
      If(istim/=0.AND.istim/=1)then
      app=char(ichar('0')+istim-1)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=30,file=Opafile,form='unformatted')
      read(30)nw,Stime,Nfreq0,Msta,Nrho,NTp,dumWavel,dumFreq,dumFreqmn,T
     *pTab,RhoTab,hpbansc1
      close(30)
      close(29)
      Open(unit=29,file=Opafile(1:Length(Opafile)-1)//'ab',form='unforma
     *tted')
      read(29)hpbanab1
      close(29)
      endif
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %R_haptab_chec
     *k:'
      write(*,*)'check nw,  Stime,  Nfreq0,   Msta,    Nrho,    NTp:'
      write(*,*)nw,Stime,Nfreq0,Msta,Nrho,NTp
      doim=1,Nzon/(Mzon/2)
      doiro=1,14
      doitp=1,14
      doL=1,Nfreq
      hbab=hpbanab2(L,itp,iro,im)
      hbal=hpbansc2(L,itp,iro,im)
      if(hbal>1.d50.or.hbab>1.d50)then
      write(*,'(4(a,i4),1p,2(a,e12.4))')' L=',L,' itp=',itp,' iro=',iro,
     *' im=',im,'  hbal=',hbal,'  hbab=',hbab
      endif
      enddo
      enddo
      enddo
      enddo
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %R_haptab_chec
     *k:'
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %R_haptab:'
      AS=ACARB
      ZS=ZCARB
      N=NZON*NVARS+2*KRAD
      HMIN=HMIN/UTIME
      HMAX=HMAX/UTIME
      METH=METHN
      Haplim=1.D0/(3.D0*FitTau)
      Uplim=1.D0+Haplim
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RH:'
      WRITE(4,'(''%RUN:'')')
      WRITE(4,'(//30X,''<===== HYDRODYNAMIC RUN OF MODEL '',A,          
     *          ''=====>'',/)')Sumprf
      WRITE(4,'(30X,''MASS(SOLAR)='',F8.3,7X,''RADIUS(SOLAR)='',F9.3/,  
     *            30X,''EXPLOSION ENERGY(10**50 ERG)='',1P,E12.5,/)')AMO
     *UT*UM,RADBEG,EKO+Eburst
      WRITE(4,'(30X,''<====='',33X,''=====>'',//)')
      WRITE(4,'(''  INPUT PARAMETERS     '')')
      WRITE(4,'('' EPS   = '',F15.5,9X,'' Rce   = '',1P,E15.5,A)')EPS,RC
     *E*UR/RSOL,' SOL.Rad.'
      WRITE(4,'('' HMIN  = '',1P,E15.5,A,7X,'' AMht  = '',1P,E15.5,A)')H
     *MIN*UTIME,' S',AMht,' SOL.MASS'
      WRITE(4,'('' HMAX  = '',1P,E15.5,A,7X,'' Tburst= '',1P,E15.5,A)')H
     *MAX*UTIME,' S',TBurst,' S'
      WRITE(4,'('' THRMAT= '',1P,E15.5,9X,'' Ebstht= '',1P,E15.5,A)')THR
     *MAT,EBurst,' 1e50 ergs'
      WRITE(4,'('' METH  = '',I15,9X,'' CONV  = '',L15)')METH,CONV
      WRITE(4,'('' JSTART= '',I15,9X,'' EDTM  = '',L15)')JSTART,EDTM
      WRITE(4,'('' MAXORD= '',I15,9X,'' CHNCND= '',L15)')MAXORD,CHNCND
      WRITE(4,'('' NSTA  = '',I15,9X,'' FitTau= '',1P,E15.5)')NSTA,FitTa
     *u
      WRITE(4,'('' NSTB  = '',I15,9X,'' TauTol= '',1P,E15.5)')NSTB,TauTo
     *l
      WRITE(4,'('' NOUT  = '',I15,9X,'' IOUT  = '',I15)')NOUT,IOUT
      WRITE(4,'('' TcurA = '',F15.5,9X,'' Rvis   ='',F15.5)')TcurA,Rvis
      WRITE(4,'('' TcurB = '',F15.5,9X,'' BQ    = '',1P,E15.5)')TcurB,BQ
      WRITE(4,'('' NSTMAX= '',I15,9X,'' DRT   = '',1P,E15.5)')NSTMAX,DRT
      WRITE(4,'('' XMNI  = '',1P,E15.5,A,'' NRT   = '',I15)')XMNI,' SOL.
     *MASS',NRT
      If(KNadap>0)then
      WRITE(4,'('' XNI   = '',1P,E15.5)')XNI
      WRITE(4,'('' CONTM.= '',1P,E15.5,A,'' SCAT  = '',L15)')AMNI,' SOL.
     *MASS',SCAT
      else
      WRITE(4,'('' XNifor= '',1P,E15.5)')XNifor(1)
      WRITE(4,'('' MNicor= '',1P,E15.5,A,'' SCAT  = '',L15)')AMNI,' SOL.
     *MASS',SCAT
      endif
      PRINT09999,UTP,UTIME,URHO,UFREQ
      PRINT09998,CK1,CK2,CFR,CKRAD,CFLUX,CLUM,CLUMF
09998 FORMAT(' CK1  =',1P,E12.5,'  CK2=',E12.5,'   CFR=',E12.5,'  CKRAD=
     *',E12.5/' CFLUX=',E12.5,'  CLUM=',E12.5,'  CLUMF=',E12.5)
09999 FORMAT(' UTP=',1P,E12.5,' UTIME=',E12.5,' URHO=',E12.5,' UFREQ=',E
     *12.5)
      PRINT'('' *****CK1, CK2 INCREASED'',G8.1,'' TIMES*****'')',CRAP
      WRITE(4,'('' FLOOR :''/1X,1P,10E10.2)')FLOOR
      PRINT*,' N DIFF. EQS=',N,'   N RAD. EQS=',2*KRAD
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RH:'
      IF(NSTEP==0)THEN
      H=1.D+05*HMIN
      HOLDBL=H
      YENTOT(1)=ELOST
      YENTOT(2)=0.d0
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %RF:'
      GOTO(09924,09923,09922),indfr(abs(KNadap))
      GOTO09921
09924 CONTINUE
      FREQ(1)=CS*1.D+08/(WLMAX*UFREQ)
      FREQ(NFREQ+1)=CS*1.D+08/(WLMIN*UFREQ)
      Basis=(WLMAX/WLMIN)**(1.D0/(DBLE(NFREQ)))
      doi=2,NFREQ
      FREQ(i)=FREQ(i-1)*BASIS
      enddo
      doL=1,NFREQ
      FREQMN(L)=SQRT(FREQ(L)*FREQ(L+1))
      enddo
      GOTO09921
09923 CONTINUE
      If(.NOT.petread)then
      read(2)Nfreq0,Msta,Nrho,Ntp,YATab,(Wavel(iif),iif=1,96),TpTab,RhoT
     *ab,STab,EpsBand,EppBand
      petread=.true.
      endif
      If(Nfreq>96)then
      doL=96+1,Nfreq
      Wavel(L)=Wavel(L-1)/2.d0
      enddo
      endif
      FREQMN(1:NFREQ)=CS*1.D+08/(WaveL(1:NFREQ)*UFREQ)
      FREQ(1)=0.5d0*Freqmn(1)
      FREQ(Nfreq+1)=2.d0*Freqmn(Nfreq)
      DO09920L=2,NFREQ
      FREQ(L)=0.5d0*(Freqmn(L-1)+Freqmn(L))
09920 CONTINUE
      GOTO09921
09922 CONTINUE
      DO09917i=1,NFREQ+1
      FREQ(i)=dumFREQ(i)/Ufreq
09917 CONTINUE
      DO09914L=1,NFREQ
      FREQMN(L)=dumFREQMN(L)/Ufreq
09914 CONTINUE
09921 CONTINUE
      write(*,'(2(a,i5))')' KNadap=',KNadap,'  indfr(abs(KNadap))=',indf
     *r(abs(KNadap))
      write(*,'(a,1pe12.4)')' Lam max,AA=',CS*1.D+08/(FREQ(1)*UFREQ)
      write(*,'(a,1pe12.4)')' Lam min,AA=',CS*1.D+08/(FREQ(NFREQ+1)*UFRE
     *Q)
      doL=1,NFREQ
      WEIGHT(L)=(FREQ(L+1)-FREQ(L))*FREQMN(L)**3
      IF(L<NFREQ)THEN
      DLOGNU(L)=1.D0/LOG(FREQMN(L)/FREQMN(L+1))
      ELSE
      DLOGNU(L)=1.D0/LOG(FREQMN(L)/FREQ(L+1))
      ENDIF
      enddo
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %RF:'
      NFRUS=NFREQ
      NTHICK(1:NFRUS)=ncnd
      doKM=NCND+1,NZON
      LTHICK(1:NFRUS,KM)=.FALSE.
      enddo
      ENDIF
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %R_freqob:'
      doL=1,Nfreq
      freqob(L)=freqmn(L)
      enddo
      BASIS=FREQ(2)/FREQ(1)
      write(*,'(a,1p2g12.4)')' basis,exp(-1.d0/dlognu(1)):',basis,exp(-1
     *.d0/dlognu(1))
      write(*,'(a)')' WEIGHT(L),freqob(L),testFreqob:'
      If(Mfreq>Nfreq)then
      FREQprev=FREQ(Nfreq+1)
      doL=Nfreq+1,Mfreq
      freqob(L)=freqmn(Nfreq)*exp(-dble(L-Nfreq)/dlognu(1))
      FREQnext=FREQprev*BASIS
      WEIGHT(L)=(FREQnext-FREQprev)*freqob(L)**3
      testFreqob=sqrt(FREQnext*FREQprev)
      write(*,'(1p,4g12.4)')WEIGHT(L),freqob(L),testFreqob
      FREQprev=FREQnext
      enddo
      endif
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %R_freqob:'
      toldg=-1.d0
      WRITE(4,'(''   FREQ:'',1P,10E12.5)')FREQ
      WRITE(4,'(''   FREQMN:'',1P,10E12.5)')FREQMN
      WRITE(4,'(''  WAVE BOUNDS:'',1P,11E10.3)')(CCL/FREQ(LW),LW=1,NFREQ
     *+1)
      WRITE(4,'(''  WAVES:'',1P,10E12.5)')(CCL/FREQMN(LW),LW=1,NFREQ)
      WRITE(4,'('' WEIGHT:'',1P,10E12.5)')WEIGHT
      write(*,'(3a)')'Opafile=',Opafile,'-->Entering Node %R_BANDS:'
      LubvU=NFRUS
09911 IF(.NOT.(LubvU>1.AND.FREQMN(LubvU)>(CCL/3.65D+03)))GOTO09910
      LubvU=LubvU-1
      GOTO09911
09910 CONTINUE
      LubvB=NFRUS
09908 IF(.NOT.(LubvB>1.AND.FREQMN(LubvB)>(CCL/4.4D+03)))GOTO09907
      LubvB=LubvB-1
      GOTO09908
09907 CONTINUE
      LubvV=NFRUS
09905 IF(.NOT.(LubvV>1.AND.FREQMN(LubvV)>(CCL/5.5D+03)))GOTO09904
      LubvV=LubvV-1
      GOTO09905
09904 CONTINUE
      LubvR=NFRUS
09902 IF(.NOT.(LubvR>1.AND.FREQMN(LubvR)>(CCL/7.D+03)))GOTO09901
      LubvR=LubvR-1
      GOTO09902
09901 CONTINUE
      LubvI=NFRUS
09899 IF(.NOT.(LubvI>1.AND.FREQMN(LubvI)>(CCL/9.D+03)))GOTO09898
      LubvI=LubvI-1
      GOTO09899
09898 CONTINUE
      Lyman=NFRUS
09896 IF(.NOT.(Lyman>1.AND.FREQMN(Lyman)>(CCL/912.D0)))GOTO09895
      Lyman=Lyman-1
      GOTO09896
09895 CONTINUE
      WRITE(4,*)' L UBVRI Lyman:',LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %R_BANDS:'
      LST=2
      KENTR=0
      write(*,'(2(a,1p,e12.4))')'Debug: Ry(1)=',Ry(1),' Rce=',Rce
      RHO(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,NZON
      RHO(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      IF(NSTEP==0)THEN
      CALLLOSSEN
      HOLDBL=H
      YENTOT(1)=ELOST
      YENTOT(2)=H*ELTOT
      ENDIF
      chem=0.5d0
      write(*,'(3a)')'Opafile=',Opafile,'<--Leaving  Node %R:'
      RETURN
9     WRITE(4,*)' BEGIN: error in read file ',Lunit
      stop
      ENDSUBROUTINEBEGIN
      SUBROUTINESAVEM
      IMPLICITREAL*8(A-H,O-Z)
      INTEGERLUNIT
      Character*80NFILE*80
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
      Common/FRAD/FRADJ(Mzon,Nfreq),FRADH(Mzon,Nfreq)
      if(NSTB<0)then
      Lunit=10
      NFILE=Sumprf
      else
      Lunit=12
      NFILE='test.prf'
      endif
      callStradIO('wm',Lunit,NFILE)
      if(NSTB<0)then
      Lunit=11
      NFILE=Sumcur
      else
      Lunit=13
      NFILE='test.crv'
      endif
      callStradIO('wc',Lunit,NFILE)
      if(NSTB<0)then
      Lunit=14
      NFILE=Sumprf(1:index(Sumprf,' ')-1)//'.edd'
      else
      Lunit=15
      NFILE='test.edd'
      endif
      OPEN(Lunit,FILE=NFILE,FORM='UNFORMATTED',STATUS='unknown')
      WRITE(Lunit)t,NSTEP,NZON,NCND,NFRUS,KRAD
      WRITE(Lunit)UFREQ,UTP,URM,UV,URHO,UM,UTIME
      WRITE(Lunit)(FREQMN(L),L=1,NFRUS)
      WRITE(Lunit)(WEIGHT(L),L=1,NFRUS)
      WRITE(Lunit)((EddN(K,L),K=NCND+1,NZON),L=1,NFRUS)
      WRITE(Lunit)((FRADJ(K,L),K=NCND+1,NZON),L=1,NFRUS)
      WRITE(Lunit)((FRADH(K,L),K=NCND+1,NZON),L=1,NFRUS)
      WRITE(Lunit)(HEdN(L),L=1,NFRUS)
      WRITE(Lunit)(Ry(K),K=1,NZON)
      WRITE(Lunit)(Uy(K),K=1,NZON)
      WRITE(Lunit)(Ty(K),K=1,NZON)
      CLOSE(Lunit)
      write(*,*)' SAVEM: wrote .edd file with EddN, FRADJ, FRADH'
      RETURN
      ENDSUBROUTINESAVEM
      BLOCKDATASTDATA
      IMPLICITREAL*8(A-H,O-Z)
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
      COMMON/QNRGYE/QNUC,RGASA,YELECT
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      DATATO/komax*1.D+19/,STO/komax*1.D-5/,NTO/komax*0/
      DATAQNUC,YELECT/8.8861D+6,0.5D0/
      DATAACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI/12.D0,6.D0,28.D0,14.D0,
     *56.D0,28.D0,88.861D0,1.884D0/
      DATANFUN,NJAC,NITER,NFAIL/4*0/
      END
