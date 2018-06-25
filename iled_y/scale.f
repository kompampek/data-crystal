      PARAMETER (L=5000)
      CHARACTER FIL*8
      INTEGER*2 N(L),M(L),IR(L),IT(L),N1(L),N2(L),N3(L)
      DIMENSION A(L),V(L),C(L),W(2),b(100)
      write(*,14)
   14 format(' triples file name :',$)   
      READ(*,1)FIL
    1 FORMAT(A)
      write(*,15)
   15 format(' RANGE =<Tau>+/- X%, NTEST=1 to limit output :',$)    
      read(*,*)scale2,NTEST
      open(1,file='avalues')
      read(1,*)b
      close(1)
      OPEN(1,FILE=FIL,FORM='FORMATTED')
      as=0
      bs=0
      cs=0
      do 16 i=1,99999
      read(1,*,END=17)a1,irr,itt,m1,m2,m3,w,G,cost
      as=as+a1
      cs=cs+a1*cost
      t=0.
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      if(irr.eq.1) T = tanh(0.5*a1)
   16 bs=bs+a1*t
   17 rewind(1)
      cs=cs/as
      bs=bs/as
      scale = cs/bs
      write(*,*)cs,bs,scale
      ND=0
      NUM=100
      NT=NUM/10
      top1=0.
      bot1=0.
      top2=0.
      bot2=0.
      do 8 iter=1,20
      AA=0.
      CS=0.
      I=1
      DO 2 I=1,NUM
      READ(1,*,END=10)A(I),IR(I),IT(I),N1(I),N2(I),N3(I),W,V(I),C(I)
      A(I)=scale*A(I)
      t=0.
      a1=a(i)
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      if(ir(i).eq.1) T = tanh(0.5*a1)
      top2=top2+a1*(t-c(i))**2
      bot2=bot2+a1
      N(I)=1
      AA=AA+A(I)
      CS=CS+C(I)
    2 CONTINUE
      goto 11
   10 NUM=I-1
      NT=NUM/10
      ND=1
   11 CONTINUE
      A1=AA/NUM
      CS=CS/NUM
      Top=scale2*CS
      if(num.eq.100) TOP = (scale2+1)*cs/2.
      DCS=-2.*(scale2-1.)*cs/num
      if(num.eq.100) dcs=dcs/2.
      T=0
      A2=A1*A1
      A3=A1*A2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0)  T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      WRITE(9,9)NUM,A(1),A1,A(NUM),CS,T
      WRITE(*,9)NUM,A(1),A1,A(NUM),CS,T
    9 format(I5,3f6.2,2f6.3) 
      DO 3 I=1,NUM
      II=0
      VM=-99.
      DO 4 J=1,NUM
      IF(N(J).EQ.0) GOTO 4
      IF(V(J).LE.VM) GOTO 4
      II=J
      VM=V(J)
    4 CONTINUE
      M(I)=II
      N(II)=0
    3 CONTINUE
C     FILE SORTED IN REVERSE ORDER (ARRAY M(I))   
      AA=0. 
      CS=0.
      VS=0.
      ws=0.
      NT=NUM/10
      DO 5 I=1,NUM
      dc=dcs
      wt=top+(i-1)*dc
      J=M(I)
      AA=AA+A(J)
      CS=CS+C(J)
      VS=VS+V(J)
      ws=ws+wt
      IF(MOD(I,NT).NE.0) GOTO 5
      AA=AA/NT
      CS=CS/NT
      VS=VS/NT
      ws=ws/nt
      if(NTEST.NE.1) WRITE(9,6)AA,VS,CS,WS
    6 FORMAT(F5.2,3F8.3)
      AA=0.
      CS=0.
      VS=0.
      ws=0.
    5 CONTINUE
C     WRITE NEW WT'D TRIPLES FILE    
      DO 12 I=1,NUM
      J=M(I)
      dc=dcs
      if(num.eq.100) dc=dc/2
      WT=TOP +(I-1)*DC
      if(wt.gt.0.99) wt=0.99
      K=100*WT
      A(J)=b(K)
      top1=top1+a(j)*(c(j)-wt)**2
      bot1=bot1+a(j)
   12 WRITE(10,13)A(J),IR(J),IT(J),N1(J),N2(J),N3(J),C(J),WT
   13 FORMAT(f8.3,2I3,I4,2I6,2f7.3) 
      IF(ND.EQ.1) goto 7
      NUM=NUM+NUM
      if(NUM.gt.3200) num=3200
    8 NT=NUM/10
    7 Rmin=top1/bot1
      R2=top2/bot2
      write(*,*)rmin,r2
      END
      
      
