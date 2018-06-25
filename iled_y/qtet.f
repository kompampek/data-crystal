C     TRIPLES GENERATOR
      PARAMETER (IE=2000)
      CHARACTER fil*6
      BYTE JH(IE,24),JK(IE,24),JL(IE,24),NR(IE),NRR(IE),NRT(IE,24)
      INTEGER*2 IH,IK,IL,NS,NE,NP,HH,HK,HL,NA,ISER,JSER,KSER
      INTEGER*2 MAXH,MAXK,MAXL,LAUE,IR(24,9),ITR(24,3),LMP,LSER,NZ
      DIMENSION SCT(16),ATMS(16),E(99999),LOC(500000)
      dimension p(-ie:ie)
      DATA MINK,MINL,NTRIP/2*999,0/
C
      WRITE(*,35)
   35 FORMAT(' ZON/GEN (1,2), NREF, BMIN, X, FILNAME : (2I,F,A4) ',$)
      READ(*,36)NTYPE,NOE,AMIN,CROSS,FIL
   36 FORMAT(2I5,2F12.5,A)
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
      write(*,*)scale
C
      DO 4 I = 1,NEQ
      READ(1,5)T1,HH,HK,HL,T2,KH,KK,KL,T3,LH,LK,LL
    5 FORMAT(3(F10.5,3I5),i5)
      IX=ABS(HH)+ABS(KK)+ABS(LL)
      IY=ABS(HK)+ABS(HL)+ABS(KH)+ABS(KL)+ABS(LH)+ABS(LK)
      IF(IY.EQ.0.AND.IX.NE.3) GOTO 6
      IDET = HH*KK*LL+HK*KL*LH+KH*LK*HL-LH*KK*HL-KH*HK*LL-HH*LK*KL
      IF(ABS(IDET).NE.1) GOTO 6
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
      FIL(5:6)='.E'
      OPEN(1,FILE=FIL,STATUS='OLD',FORM='FORMATTED')
      DO 20 I = 1,99999
      READ(1,*,END=9)IH,IK,IL,NS,NE,a,b,c,NP
C     lattice types are P A B C I F R (obverse hexagonal setting)
      IF(LAT.EQ.1.AND.MOD((IK+IL),2).NE.0) GOTO 21
      IF(LAT.EQ.2.AND.MOD((IH+IL),2).NE.0) GOTO 21
      IF(LAT.EQ.3.AND.MOD((IH+IK),2).NE.0) GOTO 21
      IF(LAT.EQ.4.AND.MOD((IH+IK+IL),2).NE.0) GOTO 21
      IX = MOD((IH+IK),2)+MOD((IH+IL),2)
      IF(LAT.EQ.5.AND.IX.NE.0) GOTO 21
      IF(LAT.EQ.6.AND.MOD((-IH+IK+IL),3).NE.0) GOTO 21
      IF(NE.LE.0) NE=1
      E(I)=0.001*NE
      CALL STD(IH,IK,IL,LAUE)
      MAXH=MAX(IH,MAXH)
      MAXK=MAX(IK,MAXK)
      MAXL=MAX(IL,MAXL)
      MINK=MIN(IK,MINK)
      MINL=MIN(IL,MINL)
      IF(I.GT.NOE) GOTO 20
      p(i)=0.001*np
      p(-i)=-p(i)
      NR(I)=13
      JH(I,1)=IH
      JK(I,1)=IK
      JL(I,1)=IL
   20 CONTINUE
    9 REWIND(1)
      NUMB=I-1
      IHR=MAXH+1
      IKR=MAXK-MINK+1
      DO 8 I=1,NUMB
      READ(1,*)IH,IK,IL
      LOCK=1+IH+IHR*(IK-MINK + IKR*(IL-MINL))
      LOC(LOCK)=I
      IF(I.GT.NOE) GOTO 8
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
      IX=ABS(JH(I,K)-KH) + ABS(JK(I,K)-KK) + ABS(JL(I,K)-KL)
      JX=ABS(JH(I,K)+KH) + ABS(JK(I,K)+KK) + ABS(JL(I,K)+KL)
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
      CLOSE(1)
C
      FIL(6:6)='Q'
      OPEN(1,FILE=FIL,STATUS='UNKNOWN',FORM='FORMATTED')
      TIME = SECNDS(0.0)
C
      DO 15 I=1,NOE-3
      I1=JH(I,1)
      I2=JK(I,1)
      I3=JL(I,1)
      IF(NTYPE.EQ.1.AND.NR(I).EQ.13) GOTO 15
      DO 16 J=I+1,NOE-2
      IF(NTYPE.EQ.1.AND.NR(J).EQ.13) GOTO 16
      A = SCALE*E(I)*E(J)**3
      IF(A.LT.AMIN) GOTO 15
      DO 33 JP=-1,1,2
      DO 17 JJ=1,NRR(J)
      JT=JP*NRT(J,JJ)
      J1=JP*JH(J,JJ)
      J2=JP*JK(J,JJ)
      J3=JP*JL(J,JJ)
      HH=I1+J1
      HK=I2+J2
      HL=I3+J3
      CALL STD(HH,HK,HL,LAUE)
      IF(HH.GT.MAXH.OR.HK.GT.MAXK.OR.HL.GT.MAXL) GOTO 17
      IF(HK.LT.MINK.OR.HL.LT.MINL) GOTO 17
      LOCK = 1 +HH +IHR*(HK-MINK + IKR*(HL-MINL))
      if(lock.le.0) write(*,1)hh,hk,hl,lock,0
      L = LOC(LOCK)
      IF(L.EQ.0) GOTO 17
      EHK=E(L)
      IF(EHK.EQ.0.0.OR.EHK.LT.CROSS) GOTO 17
      DO 30 K=J+1,NOE-1
      A=SCALE*E(I)*E(J)*E(K)**2
      IF(A.LT.AMIN) GOTO 17
      IF(NTYPE.EQ.1.AND.NR(K).EQ.13) GOTO 30
      DO 37 KP=-1,1,2
      DO 38 KK=1,NRR(K)
      KT=KP*NRT(K,KK)
      K1=KP*JH(K,KK)
      K2=KP*JK(K,KK)
      K3=KP*JL(K,KK)
      HH=I1+K1
      HK=I2+K2
      HL=I3+K3
      CALL STD(HH,HK,HL,LAUE)
      IF(HH.GT.MAXH.OR.HK.GT.MAXK.OR.HL.GT.MAXL) GOTO 38
      IF(HK.LT.MINK.OR.HL.LT.MINL) GOTO 38
      LOCK = 1 +HH +IHR*(HK-MINK + IKR*(HL-MINL))
      if(lock.le.0) write(*,1)hh,hk,hl,lock,0
      L = LOC(LOCK)
      IF(L.EQ.0) GOTO 38
      EHL=E(L)
      IF(EHL.EQ.0.0.OR.EHL.LT.CROSS) GOTO 38
      HH=K1+J1
      HK=K2+J2
      HL=K3+J3
      CALL STD(HH,HK,HL,LAUE)
      IF(HH.GT.MAXH.OR.HK.GT.MAXK.OR.HL.GT.MAXL) GOTO 38
      IF(HK.LT.MINK.OR.HL.LT.MINL) GOTO 38
      LOCK = 1 +HH +IHR*(HK-MINK + IKR*(HL-MINL))
      if(lock.le.0) write(*,1)hh,hk,hl,lock,0
      L = LOC(LOCK)
      IF(L.EQ.0) GOTO 38
      EKL=E(L)
      IF(EKL.EQ.0.0.OR.EKL.LT.CROSS) GOTO 38
      HH = -I1 -J1 -K1
      HK = -I2 -J2 -K2
      HL = -I3 -J3 -K3
      N1=HH
      N2=HK
      N3=HL
      CALL STD(HH,HK,HL,LAUE)
      IF(HH.GT.MAXH.OR.HK.GT.MAXK.OR.HL.GT.MAXL) GOTO 38
      IF(HK.LT.MINK.OR.HL.LT.MINL) GOTO 38
      LOCK = 1 +HH +IHR*(HK-MINK + IKR*(HL-MINL))
      if(lock.le.0) write(*,1)hh,hk,hl,lock
      L = LOC(LOCK)
      IF(L.LE.K.OR.L.GT.NOE) GOTO 38
      IF(NTYPE.EQ.1.AND.NR(L).EQ.13) GOTO 38
      A = SCALE*E(I)*E(J)*E(K)*E(L)
      IF(A.LT.AMIN) GOTO 38
      DO 18 LL=1,NRR(L)
      DO 34 LP=-1,1,2
      IX=ABS(LP*JH(L,LL)-N1)+ABS(LP*JK(L,LL)-N2)+ABS(LP*JL(L,LL)-N3)
      IF(IX.NE.0) GOTO 34
      LMP = MOD((JT + KT +LP*NRT(L,LL) +20376),12)
      JSER=J*JP
      KSER=K*KP
      LSER=L*LP
      GOTO 31
   34 CONTINUE
   18 CONTINUE
      WRITE(*,19)L,i1,i2,i3,j1,j2,j3,k1,k2,k3,n1,n2,n3
   19 FORMAT(' SYMMETRY INDEX FAILURE !!!!',i5,4(i4,2i3))
      GOTO 14
   31 IF(LAUE-2) 23,24,25
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
      GOTO 23
C     4/M
   27 IF(I3.EQ.0.AND.J3.LT.0) GOTO 32
      IF(I1.EQ.0.AND.I2.EQ.0.AND.J1.LT.0) GOTO 32
   23 arg=cos(p(i)+p(jser)+p(kser)+p(lser)+LMP*.52359878)
      NTRIP=NTRIP+1
      ISER=I
      WRITE(1,111)A,LMP,NZ,ISER,JSER,KSER,LSER,arg
  111 FORMAT(f8.4,6i5,f6.2)
   32 CONTINUE
   38 CONTINUE
   37 CONTINUE
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
      TIME = SECNDS(TIME)
      WRITE(*,29)NTRIP,TIME
   29 FORMAT(I8,' QUARTETS GENERATED IN ',F8.1,' SECONDS ')
      END
      SUBROUTINE STD(IH,IK,IL,LAU)
C     ******************************************************************
C     SUBROUTINE TO DETERMINE STANDARD PARENT FORM
C     ******************************************************************
      INTEGER*2 IH,IK,IL,LAU
      LAUE=ABS(LAU)
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
    3 IK = ABS(IK)
      IF(IH.GT.0) RETURN
      IH = -IH
      IL = -IL
      IF(IH.EQ.0) IL=ABS(IL)
      RETURN
C     ORTHORHOMBIC SYMMETRY
    4 IH = ABS(IH)
      IK = ABS(IK)
      IL = ABS(IL)
      RETURN
    1 IF(LAUE-5) 9,9,10
C     4/M SYMMETRY
    9 IL = ABS(IL)
      IF(IH*IK) 12,13,14
   12 IX = ABS(IH)
      IH = ABS(IK)
      IK = IX
      GOTO 11
   13 IK = ABS(IH) + ABS(IK)
      IH = 0
      GOTO 11
   14 IH = ABS(IH)
      IK = ABS(IK)
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
      IK=ABS(IH)
      IH=0
      GOTO 21
   20 IF(IK.LT.0) IL=-IL
      IK=ABS(IK)
      GOTO 21
   19 IF(IH.GT.0) GOTO 21
      IH=-IH
      IK=-IK
      IL=-IL
   21 IF(LAUE.EQ.6) RETURN
      IF(LAUE.GE.8) IL=ABS(IL)
      IF(LAUE.EQ.8) RETURN
      IF(IH.LT.IK) RETURN
      IX=IH
      IH=IK
      IK=IX
      IF(LAU.EQ.-7) IL=-IL
      IF(LAU.EQ.-7.AND.IH.EQ.IK) IL=ABS(IL)
      RETURN
C     CUBIC SYMMETRY
   16 IH=ABS(IH)
      IK=ABS(IK)
      IL=ABS(IL)
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
