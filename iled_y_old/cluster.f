      CHARACTER*6 FILN
      DIMENSION NAM(5000),IH(1000),IK(1000),IL(1000),NSET(50)
      DIMENSION RMIN(5000),COSAV(5000),P(5000,900),Q(900)
      DIMENSION A(900),B(900)
C      WRITE(*,12)
   12 FORMAT(' NAME.E : ',$)
C      READ(*,13) FILN
   13 FORMAT(A)   
      OPEN(1,FILE='iled.E')
      PI=ACOS(-1.0)
C      WRITE(*,14)
   14 FORMAT(' # PHASES,RMIN THRESHOLD/NAMB SEEDS : ',$)   
      NP=250
      THRES=0.46
C      READ(*,*)NP,THRES
      NSET0=99999
      N=1
      OPEN(2,FILE='NUMBERS')
      DO 15 I=1,50
      READ(2,*,END=23)NSET(I)
      IF(NSET(I).EQ.0) GOTO 23
   15 NSET0=MIN(NSET0,NSET(I))
   23 NSETS=I-1
      DO 2 I=1,NP
    2 READ(1,*)IH(I),IK(I),IL(I)
      CLOSE(1) 
      OPEN(1,FILE='SAVE')
      N=1
      DO 1 I=1,100
      READ(1,*,END=8)NAM(N),NC,RMIN(N),COSAV(N),DEL
      READ(1,6)(P(N,J),J=1,NP)
    6 FORMAT(20F6.2)
      IF(DEL.GT.50.) GOTO 1
   11 FORMAT(I7,I4,2F7.3)
      N=N+1   
    1 CONTINUE
    8 NUM=N-1
      CLOSE(1)
      WRITE(*,17)NSETS,NUM,NSET0
   17 FORMAT(I5,' SETS IN BASES OF ',I5,' TOTAL, NAMIN = ',I6)   
C
      NAV=1
      DO 5 I=1,NUM
      NI=NAM(I)
      IF(NI.NE.NSET0) GOTO 21
      DO 20 J=1,NP
      A(J)=COS(P(I,J))
   20 B(J)=SIN(P(I,J))
      WRITE(*,*)NI,99
      WRITE(10,11)NI
      WRITE(10,6)(P(I,J),J=1,NP)
   21 ID=0
      DO 18 J=1,NSETS
   18 IF(NI.EQ.NSET(J))ID=1
      IF(ID.NE.1) GOTO 5
C      
      DO 10 J=1,NUM
      NJ=NAM(J)
      IF(NI.EQ.NJ) GOTO 10
      PMIN=999.
      DO 4 I0=0,15
      II=I0
      I1=IAND(II,1)
      II=II/2
      I2=IAND(II,1)
      II=II/2
      I3=IAND(II,1)
      I4=II/2
      IF(I4.EQ.0) I4=-1
      S=0
      DO 7 K=1,NP
      IT=I1*IH(K)+I2*IK(K)+I3*IL(K)
      IT=IAND(IT,1)
      DP=COS(P(I,K)-I4*P(J,K)+IT*PI)
      DP=57.296*ACOS(DP)
    7 S=S+DP
      S=S/NP
      IF(S.GT.PMIN) GOTO 4
      PMIN=S
      J1=I1
      J2=I2
      J3=I3
      J4=I4
    4 CONTINUE
C    
      ID=0
      DO 19 K=1,NSETS
   19 IF(NJ.EQ.NSET(K)) ID=1
      IF(PMIN.GT.50.AND.ID.EQ.0) GOTO 10
      IF(ID.EQ.0) WRITE(8,3)NI,NJ,J1,J2,J3,J4,PMIN,0
      IF(ID.EQ.1.AND.NI.LT.NJ) WRITE(9,3)NI,NJ,J1,J2,J3,J4,PMIN,1
    3 FORMAT(2I7,4I3,F7.1,I5)
      IF(ID.EQ.0.OR.NI.NE.NSET0) GOTO 10
      DO 9 K=1,NP
      IT=IH(K)*J1+IK(K)*J2+IL(K)*J3
      IT=IAND(IT,1)
    9 Q(K)=J4*P(J,K)+IT*PI
      WRITE(10,11)NJ
      WRITE(10,6)(Q(K),K=1,NP)
      WRITE(*,16)NI,NJ,PMIN,RMIN(I),COSAV(I)
   16 FORMAT(2I6,F6.1,2F8.3)   
      DO 24 K=1,NP
      A(K)=A(K)+COS(Q(K))
   24 B(K)=B(K)+SIN(Q(K))
      NAV=NAV+1  
   10 CONTINUE
    5 CONTINUE
C    
      DO 25 I=1,NP
   25 Q(I)=ATAN2(B(I),A(I))
      WRITE(11,*)NSET(1),NC,0.0,0.0
      WRITE(11,6) (Q(I),I=1,NP)
      REWIND(10)
      B1=0
      B2=0
      DO 26 I=1,NSETS
      READ(10,*,END=22)NI
      READ(10,6)(P(I,J),J=1,NP)
      A1=0
      A2=0
      DO 28 J=1,NP
      DP=COS(P(I,J)-Q(J))
      DP=57.296*ACOS(DP)
      A1=A1+DP
   28 A2=A2+DP*DP
      B1=B1+A1
      B2=B2+A2
      AV=A1/NP
      SIG=SQRT(A2/NP -AV*AV)
      WRITE(*,*)NI,AV,SIG,NAV 
   26 CONTINUE
   22 MUM=NSETS*NP
      B1=B1/MUM
      SIG=SQRT(B2/MUM -B1*B1)
      WRITE(*,*)B1,SIG   
      END
      
