C     TRIPLES GENERATOR
      PARAMETER (IE=14000)
      CHARACTER fil*6
      BYTE JH(IE,24),JK(IE,24),JL(IE,24),NR(IE),NRR(IE),NRT(IE,24)
      INTEGER*4 IH,IK,IL,NS,NE,NP,HH,HK,HL,LH,LK,LL,NA,I,JSER,KSER,NZ
      INTEGER*4 MAXH,MAXK,MAXL,LAUE,IR(24,9),ITR(24,3),LMP,LAT
      DIMENSION SCT(16),ATMS(16),E(IE),LOC(500000)
      dimension p(-ie:ie)
      DATA MAXH,MAXK,MAXL,MINH,MINK,MINL/3*-99,3*999/
C
      WRITE(*,37)
   37 FORMAT(' FILENAME : ',$)
      READ(*,38)fil
   38 FORMAT(A)     
      WRITE(*,35)
   35 FORMAT(' SG1/ZON/GEN (1,2,3), NREF, AMIN: ',$)
      READ(*,*)NTYPE,NOE,AMIN
   36 FORMAT(2I6,F12.5,A4)
      OPEN(UNIT=1,FILE='sym.dat',STATUS='OLD',FORM='FORMATTED')
      READ(1,1) NEQ,LAUE,ICENT,LAT
    1 FORMAT(5I5,F5.2,f8.3)
      READ(1,2) SCT
      READ(1,2) ATMS
    2 FORMAT(16F5.0)
     
      Z2=0.0  !!!YV
      Z3=0.0   !!!YV
      DO 3 I = 1,16
      Z2 = Z2 + (2-ICENT)*NEQ*ATMS(I)*SCT(I)**2
    3 Z3 = Z3 + (2-ICENT)*NEQ*ATMS(I)*SCT(I)**3
      SCALE = 2.0*Z3/(Z2**1.5)
      write(*,*)'SCALE===',SCALE
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
      DO 20 I = 1,NOE
      READ(1,*,END=9)IH,IK,IL,NS,NE,a,b,c,NP
CCCC      READ(1,*,END=9)IH,IK,IL,NS,NE,a,b,c,d,NP   !!YV
C     lattice types are P A B C I F R (obverse hexagonal setting)
      IF(LAT.EQ.1.AND.MOD((IK+IL),2).NE.0) GOTO 21
      IF(LAT.EQ.2.AND.MOD((IH+IL),2).NE.0) GOTO 21
      IF(LAT.EQ.3.AND.MOD((IH+IK),2).NE.0) GOTO 21
      IF(LAT.EQ.4.AND.MOD((IH+IK+IL),2).NE.0) GOTO 21
      IX = MOD((IH+IK),2)+MOD((IH+IL),2)
      IF(LAT.EQ.5.AND.IX.NE.0) GOTO 21
      IF(LAT.EQ.6.AND.MOD((-IH+IK+IL),3).NE.0) GOTO 21
      E(I)=0.001*NE
      if(I.lt.50) write(*,*),'E====',E(I)   !!YV
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
      IHR=MAXH+1
      IKR=MAXK-MINK+1
      if(ntype.eq.1) fil(6:6)='S'
      if(ntype.eq.2) fil(6:6)='Z'
      if(ntype.eq.3) fil(6:6)='T'
      OPEN(1,FILE=fil,STATUS='UNKNOWN',FORM='FORMATTED')
      DO 8 I=1,NOE
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
      DO 15 I=1,NOE-2
      IF(NTYPE.EQ.2.AND.NR(I).EQ.13) GOTO 15
      IP1=I+1
      NOEM1=NOE-1
      IF(NTYPE.EQ.1) IP1=I
      IF(NTYPE.EQ.1) NOEM1=I
      DO 16 J=IP1,NOEM1
      IF(NTYPE.EQ.2.AND.NR(J).EQ.13) GOTO 16
      A1 = SCALE*E(I)*E(J)
      IF((A1*E(J)).LT.AMIN) GOTO 15
      IR2=1
      if(ntype.eq.1.and.nr(i).eq.13) IR2=-1
      DO 33 LP=-1,IR2,2
      DO 17 JJ=1,NRR(J)
      JT = LP*NRT(J,JJ)
      J1=LP*JH(J,JJ)
      J2=LP*JK(J,JJ)
      J3=LP*JL(J,JJ)
      HH = -JH(I,1) -J1
      HK = -JK(I,1) -J2
      HL = -JL(I,1) -J3
      i1=hh
      i2=hk
      i3=hl
      CALL STD(HH,HK,HL,LAUE)
      IF(HH.GT.MAXH.OR.HK.GT.MAXK.OR.HL.GT.MAXL) GOTO 17
      IF(HH.LT.MINH.OR.HK.LT.MINK.OR.HL.LT.MINL) GOTO 17
      LOCK = 1 +HH +IHR*(HK-MINK + IKR*(HL-MINL))
      if(lock.lt.0) write(*,1)hh,hk,hl,lock
      K = LOC(LOCK)
      IF(NTYPE.NE.1.AND.K.LE.J) GOTO 17
      IF(NTYPE.NE.3.AND.NR(K).EQ.13) GOTO 17
      A = A1*E(K)
      IF(A.LT.AMIN) GOTO 17
      write(*,*)'A,I,J,K,SCALE=====',A,I,J,K,SCALE
      DO 18 KK=1,NRR(K)
      DO 34 MP=-1,1,2
      IX=IABS(MP*JH(K,KK)-I1)+IABS(MP*JK(K,KK)-I2)+IABS(MP*JL(K,KK)-I3)
      IF(IX.NE.0) GOTO 34
      LMP = MOD((JT + MP*NRT(K,KK) +20376),12)
      JSER=J*LP
      KSER=K*MP
      GOTO 31
   34 CONTINUE
   18 CONTINUE
      WRITE(*,19)K,i1,i2,i3,jh(k,1),jk(k,1),jl(k,1)
   19 FORMAT(' SYMMETRY INDEX FAILURE !!!!',7i5)
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
   24 IF(JK(I,1).EQ.0.AND.J2.LT.0) GOTO 32
      IF(JH(I,1).EQ.0.AND.JL(I,1).EQ.0.AND.J1.LT.0) GOTO 32
      GOTO 23
C     ORTHORHOMBIC
   26 IF(JH(I,1).EQ.0.AND.J1.LT.0) GOTO 32
      IF(JK(I,1).EQ.0.AND.J2.LT.0) GOTO 32
C     4/M
   27 IF(JL(I,1).EQ.0.AND.J3.LT.0) GOTO 32
      IF(JH(I,1).EQ.0.AND.JK(I,1).EQ.0.AND.J1.LT.0) GOTO 32
   23 ang=p(i)+p(jser)+p(kser)+LMP*.52359878
      arg = cos(ang)
      ang = 57.296*ang
      ang = mod(ang,360.)
      ng = ang
      NTRIP=NTRIP+1
C     IF(NTYPE.EQ.1.AND.(I*JSER).GT.0.AND.NR(K).EQ.13) GOTO 32
      IF(NTYPE.EQ.1) A=0.5*SCALE*(E(I)*E(I)-1.0)*E(K)
      IF(A.LT.AMIN) GOTO 32
      NA=1000.*A
      IH=NR(I)+NR(J)+NR(K)
      LPM=MOD(IH,12)
      IRES=0
      IF(LPM.EQ.0) IRES=1
      IF(LPM.EQ.6) IRES=7
      WRITE(1,39)A,IRES,LMP,I,JSER,KSER,0,arg,ng
   39 FORMAT(f8.3,2i3,i6,3i7,f8.3,i5)   
   32 CONTINUE
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
      TIME=T1-T0
      WRITE(*,29)NTRIP,TIME
   29 FORMAT(I8,' TRIPLES GENERATED IN ',F8.1,' SECONDS ')
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
