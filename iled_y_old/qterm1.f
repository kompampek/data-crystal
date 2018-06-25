C     TRIPLES GENERATOR
      PARAMETER (IE=24000)
      CHARACTER fil*6
      BYTE JH(IE,24),JK(IE,24),JL(IE,24),NR(IE),NRR(IE),NRT(IE,24)
      INTEGER*4 IH,IK,IL,NS,NE,NP,HH,HK,HL,LH,LK,LL,NA,I,JSER,KSER,NZ
      INTEGER*4 MAXH,MAXK,MAXL,LAUE,IR(24,9),ITR(24,3),LMP,LAT
      DIMENSION SCT(16),ATMS(16),E(IE),LOC(500000)
      dimension p(-ie:ie)
      DATA NTRIP,MAXH,MAXK,MAXL,MINH,MINK,MINL/4*0,3*999/
C
      WRITE(*,37)
   37 FORMAT(' FILENAME : ',$)
      READ(*,38)fil
   38 FORMAT(A)     
      WRITE(*,35)
   35 FORMAT('  NREF, AMIN: ',$)
      READ(*,*)NOE,AMIN
   36 FORMAT(2I6,F12.5,A4)
      OPEN(UNIT=1,FILE='sym.dat',STATUS='OLD',FORM='FORMATTED')
      READ(1,1) NEQ,LAUE,ICENT,LAT
    1 FORMAT(5I5,F5.2,f8.3)
      READ(1,2) SCT
      READ(1,2) ATMS
    2 FORMAT(16F5.0)
      DO 3 I = 1,16
      Z2 = Z2 + (2-ICENT)*NEQ*ATMS(I)*SCT(I)**2
    3 Z4 = Z4 + (2-ICENT)*NEQ*ATMS(I)*SCT(I)**4
      SCALE = 2.0*Z4/(Z2**2)
C
      DO 4 I = 1,NEQ
      READ(1,5)T1,HH,HK,HL,T2,KH,KK,KL,T3,LH,LK,LL
    5 FORMAT(3(F10.5,3I5),i5)
      IX=IABS(HH)+IABS(KK)+IABS(LL)
      IY=IABS(HK)+IABS(HL)+IABS(KH)+IABS(KL)+IABS(LH)+IABS(LK)
      IF(IY.EQ.0.AND.IX.NE.3) GOTO 6
      IDET = HH*KK*LL+HK*KL*LH+KH*LK*HL-LH*KK*HL-KH*HK*LL-HH*LK*KL
      IF(IABS(IDET).NE.1) GOTO 6
      IR(I,1) = (KK*LL-LK*KL)/IDET
      IR(I,2) = -(KH*LL-LH*KL)/IDET
      IR(I,3) = (KH*LK-KK*LH)/IDET
      IR(I,4) = -(HK*LL-HL*LK)/IDET
      IR(I,5) = (HH*LL-LH*HL)/IDET
      IR(I,6) = -(HH*LK-LH*HK)/IDET
      IR(I,7) = (HK*KL-KK*HL)/IDET
      IR(I,8) = -(HH*KL-HL*KH)/IDET
      IR(I,9) = (HH*KK-HK*KH)/IDET
      ITR(I,1)=NINT(12*T1)
      ITR(I,2)=NINT(12*T2)
    4 ITR(I,3)=NINT(12*T3)
      CLOSE(1)
C
      fil(5:5)='.'
      fil(6:6)='E'
      OPEN(1,FILE=fil,STATUS='OLD',FORM='FORMATTED')
      DO 20 I = 1,14000
      READ(1,*,END=9)IH,IK,IL,NS,NE,a,b,c,NP
C     lattice types are P A B C I F R (obverse hexagonal setting)
      IF(LAT.EQ.1.AND.MOD((IK+IL),2).NE.0) GOTO 21
      IF(LAT.EQ.2.AND.MOD((IH+IL),2).NE.0) GOTO 21
      IF(LAT.EQ.3.AND.MOD((IH+IK),2).NE.0) GOTO 21
      IF(LAT.EQ.4.AND.MOD((IH+IK+IL),2).NE.0) GOTO 21
      IX = MOD((IH+IK),2)+MOD((IH+IL),2)
      IF(LAT.EQ.5.AND.IX.NE.0) GOTO 21
      IF(LAT.EQ.6.AND.MOD((-IH+IK+IL),3).NE.0) GOTO 21
      E(I)=0.001*NE
      CALL STD(IH,IK,IL,LAUE)
      MAXH=MAX(IH,MAXH)
      MAXK=MAX(IK,MAXK)
      MAXL=MAX(IL,MAXL)
      MINH=MIN(IH,MINH)
      MINK=MIN(IK,MINK)
      MINL=MIN(IL,MINL)
      p(i)=0.001*np
      p(-i)=-p(i)
      NR(I)=13
      JH(I,1)=IH
      JK(I,1)=IK
   20 JL(I,1)=IL
    9 CLOSE(1)
      NREF=I-1
      IHR=MAXH+1
      IKR=MAXK-MINK+1
      fil(6:6)='Q'
      write(*,*)maxh,maxk,maxl,mink,minl,NREF,scale
      OPEN(1,FILE=fil,STATUS='UNKNOWN',FORM='FORMATTED')
      DO 8 I=1,NREF
      IH=JH(I,1)
      IK=JK(I,1)
      IL=JL(I,1)
      LOCK=1+IH+IHR*(IK-MINK + IKR*(IL-MINL))
      LOC(LOCK)=I
      NCT=1
      BR=-999999.
      AR=-999999.
      DO 10 J=2,NEQ
      KH=IH*IR(J,1) +IK*IR(J,2) + IL*IR(J,3)
      KK=IH*IR(J,4) +IK*IR(J,5) + IL*IR(J,6)
      KL=IH*IR(J,7) +IK*IR(J,8) + IL*IR(J,9)
      IT=0
      JT=0
      NX=0
      DO 11 K=1,J-1
      IX=IABS(JH(I,K)-KH) + IABS(JK(I,K)-KK) + IABS(JL(I,K)-KL)
      JX=IABS(JH(I,K)+KH) + IABS(JK(I,K)+KK) + IABS(JL(I,K)+KL)
      IF(IX.EQ.0.OR.JX.EQ.0) NX=1
      IF(JX.EQ.0) JT=K
   11 IF(IX.EQ.0) IT=K
      IF(IT.NE.0.AND.AR.EQ.-999999.) AR = (IH*(ITR(J,1)-ITR(IT,1)) 
     * + IK*(ITR(J,2)-ITR(IT,2)) + IL*(ITR(J,3)-ITR(IT,3)))
      IF(JT.NE.0.AND.BR.EQ.-999999.) BR = (IH*(ITR(J,1)+ITR(JT,1))
     * + IK*(ITR(J,2)+ITR(JT,2)) + IL*(ITR(J,3)+ITR(JT,3)))
      IF(NX.EQ.1) GOTO 10
      NCT=NCT+1
      JH(I,NCT)=KH
      JK(I,NCT)=KK
      JL(I,NCT)=KL
      NRT(I,NCT)=MOD((KH*ITR(J,1) +KK*ITR(J,2) +KL*ITR(J,3)),12)
   10 CONTINUE
      IF(BR.NE.-999999.) NR(I)=NINT(AMOD(BR,12.0))
      IF(AR.NE.-999999.0.AND.AMOD(AR,12.0).NE.0.0) GOTO 13
      IF(ICENT.EQ.0) NR(I)=0
      NRR(I)=NCT
    8 CONTINUE
C
      CALL CPU_TIME(T0)
      DO 15 I=1,NOE-3
      I1=JH(I,1)
      I2=JK(I,1)
      I3=JL(I,1)
      DO 16 J=I+1,NOE-2
      DO 33 LP=-1,1,2
      DO 17 JJ=1,NRR(J)
      JT = LP*NRT(J,JJ)
      J1=LP*JH(J,JJ)
      J2=LP*JK(J,JJ)
      J3=LP*JL(J,JJ)
      i4=i1+j1
      i5=i2+j2
      i6=i3+j3
      CALL STD(i4,i5,i6,Laue)
      IF(i4.GT.MAXH.OR.i5.GT.MAXK.OR.i6.GT.MAXL) GOTO 45
      IF(i4.LT.MINH.OR.i5.LT.MINK.OR.i6.LT.MINL) GOTO 45
      LOCK = 1 +i4 +IHR*(i5-MINK + IKR*(i6-MINL))
      M=LOC(LOCK)
      IF(M.EQ.0.OR.M.GT.NOE) GOTO 17
      EHK=E(M)
      NUMB=1
   45 CONTINUE   
      DO 30 K=J+1,NOE-1
      DO 40 MP=-1,1,2
      DO 41 KK=1,NRR(K)
      KT = MP*NRT(K,KK)
      K1=MP*JH(K,KK)
      K2=MP*JK(K,KK)
      K3=MP*JL(K,KK)
      HH = -I1 -J1 -K1
      HK = -I2 -J2 -K2
      HL = -I3 -J3 -K3
      j4=hh
      j5=hk
      j6=hl
      CALL STD(HH,HK,HL,LAUE)
      IF(HH.GT.MAXH.OR.HK.GT.MAXK.OR.HL.GT.MAXL) GOTO 41
      IF(HH.LT.MINH.OR.HK.LT.MINK.OR.HL.LT.MINL) GOTO 41
      LOCK = 1 +HH +IHR*(HK-MINK + IKR*(HL-MINL))
      if(lock.lt.0) write(*,1)hh,hk,hl,lock
      L = LOC(LOCK)
      IF(L.LE.K.OR.L.GT.NOE) GOTO 41
      A = SCALE*E(I)*E(J)*E(K)*E(L)
      DO 18 LL=1,NRR(L)
      DO 34 NP=-1,1,2
      IX=IABS(NP*JH(L,LL)-j4)+IABS(NP*JK(L,LL)-j5)+IABS(NP*JL(L,LL)-j6)
      IF(IX.NE.0) GOTO 34
      LMP = MOD((JT + KT + NP*NRT(L,LL) +20376),12)
      JSER=J*LP
      KSER=K*MP
      LSER=L*NP
      GOTO 42
   34 CONTINUE
   18 CONTINUE
      WRITE(*,19)K,i1,i2,i3,jh(k,1),jk(k,1),jl(k,1)
   19 FORMAT(' SYMMETRY INDEX FAILURE !!!!',7i5)
      GOTO 14
   42 j4=i1+k1
      j5=i2+k2
      j6=i3+k3
      NUM=NUMB
      CALL STD(j4,j5,j6,Laue)
      IF(j4.GT.MAXH.OR.j5.GT.MAXK.OR.j6.GT.MAXL) GOTO 41
      IF(j4.LT.MINH.OR.j5.LT.MINK.OR.j6.LT.MINL) GOTO 41
      LOCK = 1 +j4 +IHR*(j5-MINK + IKR*(j6-MINL))
      M=LOC(LOCK)
      IF(M.EQ.0.OR.M.GT.NOE) GOTO 41
      EHL=E(M)
      NUM=NUM+1
   43 k4=k1+j1
      k5=k2+j2
      k6=k3+j3
      CALL STD(k4,k5,k6,Laue)
      IF(k4.GT.MAXH.OR.k5.GT.MAXK.OR.k6.GT.MAXL) GOTO 41
      IF(k4.LT.MINH.OR.k5.LT.MINK.OR.k6.LT.MINL) GOTO 41
      LOCK = 1 +k4 +IHR*(k5-MINK + IKR*(k6-MINL))
      M=LOC(LOCK)
      IF(M.EQ.0.OR.M.GT.NOE) GOTO 41
      EKL=E(M)
      NUM=NUM+1
   44 CROS=EHK**2 +EHL**2 +EKL**2 -2.
      AA=A*CROS
      BEE=abs(AA)
      IF(AA.ge.0.0.and.BEE.LT.AMIN) GOTO 32
      IF(AA.le.0.0.and.BEE.LT.AMIN2) GOTO 32
      IF(LAUE-2) 23,24,25
   25 IF(LAUE-4) 26,27,28
   28 IF(LAUE.EQ.5.AND.JH(I,1).EQ.JK(I,1).AND.J1.GT.J2) GOTO 32
      IF(LAUE.EQ.5.OR.LAUE.EQ.10.OR.LAUE.EQ.11) GOTO 26
      JT=(LAUE-8)*(LAUE-9)
      IF(JT.EQ.0.AND.JL(I,1).EQ.0.AND.J3.LT.0) GOTO 32
      IF(JT.EQ.0.AND.JH(I,1).EQ.0.AND.JK(I,1).EQ.0.AND.J1.LT.0) GOTO 32
      GOTO 23
C     MONOCLINIC
   24 IF(I1.EQ.0.AND.J1.LT.0) GOTO 32
      IF(I2.EQ.0.AND.J2.LT.0) GOTO 32
      IF(I1.EQ.0.AND.J1.EQ.0.AND.K1.LT.0) GOTO 32
      IF(I2.EQ.0.AND.J2.EQ.0.AND.K2.LT.0) GOTO 32
      GOTO 23
C     ORTHORHOMBIC
   26 IF(I1.EQ.0.AND.J1.LT.0) GOTO 32
      IF(I2.EQ.0.AND.J2.LT.0) GOTO 32
      IF(I3.EQ.0.AND.J3.LT.0) GOTO 32
      IF(I1.EQ.0.AND.J1.EQ.0.AND.K1.LT.0) GOTO 32
      IF(I2.EQ.0.AND.J2.EQ.0.AND.K2.LT.0) GOTO 32
      IF(I3.EQ.0.AND.J3.EQ.0.AND.K3.LT.0) GOTO 32
C     4/M
   27 IF(I3.EQ.0.AND.J3.LT.0) GOTO 32
      IF(I1.EQ.0.AND.I2.EQ.0.AND.J1.LT.0) GOTO 32
   23 arg=cos(p(i)+p(jser)+p(kser)+p(lser)+LMP*.52359878)
      IH=NR(I)+NR(J)+NR(K)+NR(L)
      LPM=MOD(IH,12)
      IRES=0
      IF(LPM.EQ.0) IRES=1
      IF(LPM.EQ.6) IRES=7
      NTRIP=NTRIP+1
      if(AA.eq.0.0) write(*,*)A,CROS,EHK,EKL,EHL
      WRITE(1,39)AA,IRES,LMP,I,JSER,KSER,LSER,arg,NUM
   39 FORMAT(f8.3,6i6,f7.3,i3,3f6.2)   
   32 CONTINUE
   41 CONTINUE
   40 CONTINUE
   30 CONTINUE
   17 CONTINUE
   33 CONTINUE
   16 CONTINUE
   15 CONTINUE
C     
      GOTO 14
   21 WRITE(*,22) IH,IK,IL
   22 FORMAT(' REFLECTION ',3I5,' VIOLATES LATTICE CONDITION ')
      GOTO 14
   13 WRITE(*,12)IH,IK,IL,NE
   12 FORMAT(1X,'*** EXTINCTION $$$',3I5,'   NE = ',I6)
      GOTO 14
    6 WRITE(*,7)I,IDET
    7 FORMAT(' SYMMETRY CARD #',i2,' IS IN ERROR, DETERMINANT = ',i5)
   14 CONTINUE
      CALL CPU_TIME(T1)
      TIME=(T1-T0)/60.
      WRITE(*,29)NTRIP,TIME
   29 FORMAT(I8,' QUARTETS GENERATED IN ',F8.2,' MINUTES ')
      END
      SUBROUTINE STD(IH,IK,IL,LAU)
C     ******************************************************************
C     SUBROUTINE TO DETERMINE STANDARD PARENT FORM
C     ******************************************************************
      INTEGER*4 IH,IK,IL,LAU
      LAUE=IABS(LAU)
      IF(LAUE.GT.3) GOTO 1
      IF(LAUE-2) 2,3,4
C     TRICLINIC SYMMETRY
    2 IF(IH) 5,6,7
    6 IF(IK) 5,8,7
    8 IF(IL) 5,7,7
    5 IH = -IH
      IK = -IK
      IL = -IL
    7 RETURN
C     MONOCLINIC SYMMETRY -  B-AXIS UNIQUE
    3 IK = IABS(IK)
      IF(IH.GT.0) RETURN
      IH = -IH
      IL = -IL
      IF(IH.EQ.0) IL=IABS(IL)
      RETURN
C     ORTHORHOMBIC SYMMETRY
    4 IH = IABS(IH)
      IK = IABS(IK)
      IL = IABS(IL)
      RETURN
    1 IF(LAUE-5) 9,9,10
C     4/M SYMMETRY
    9 IL = IABS(IL)
      IF(IH*IK) 12,13,14
   12 IX = IABS(IH)
      IH = IABS(IK)
      IK = IX
      GOTO 11
   13 IK = IABS(IH) + IABS(IK)
      IH = 0
      GOTO 11
   14 IH = IABS(IH)
      IK = IABS(IK)
   11 IF(LAUE.EQ.4) RETURN
      IX = MIN(IH,IK)
      IK = MAX(IH,IK)
      IH = IX
      RETURN
   10 IF(LAUE-10) 15,16,16
C     BAR-3 SYMMETRY
   15 IF(IH*IK) 17,18,19
   17 IX = IH
      IH = -IK
      IK = IX+IK
      IL = -IL
      GOTO 15
   18 IF(IH.EQ.0) GOTO 20
      IF(IH.GT.0) IL=-IL
      IK=IABS(IH)
      IH=0
      GOTO 21
   20 IF(IK.LT.0) IL=-IL
      IK=IABS(IK)
      GOTO 21
   19 IF(IH.GT.0) GOTO 21
      IH=-IH
      IK=-IK
      IL=-IL
   21 IF(LAUE.EQ.6) RETURN
      IF(LAUE.GE.8) IL=IABS(IL)
      IF(LAUE.EQ.8) RETURN
      IF(IH.LT.IK) RETURN
      IX=IH
      IH=IK
      IK=IX
      IF(LAU.EQ.-7) IL=-IL
      IF(LAU.EQ.-7.AND.IH.EQ.IK) IL=IABS(IL)
      RETURN
C     CUBIC SYMMETRY
   16 IH=IABS(IH)
      IK=IABS(IK)
      IL=IABS(IL)
      IX=MIN(IH,IK)
      MINH=MIN(IX,IL)
   23 IX=IH
      IH=IK
      IK=IL
      IL=IX
      IF(IH.NE.MINH) GOTO 23
      IF(IH.EQ.MINH.AND.IK.NE.MINH.AND.IL.EQ.MINH) GOTO 23
      IF(LAUE.EQ.10) RETURN
      IF(IK.LT.IL) RETURN
      IX=IK
      IK=IL
      IL=IX
      RETURN
      END
