      parameter mp=600,nd=25000,md=2400
      REAL*4 DET,E,DA,DB,W,ah,bh,ahs,bhs,a1,aa1,a2,aa2,ab,ap,bp
      REAL*4 xa,xb,ya,yb,za,zb,arg,ff,hh,hk,hl
      integer*2 ih(nd),ik(nd),il(nd),l(mp),nt(24,9),ns,ne,np,num(10)
      dimension x(8,mp),y(8,mp),z(8,mp),g(mp),beta(mp,6),F(nd,4),dw(10)
      dimension c(9,5),fj(5),tx(24,3),e(md,md),w(md),da(md),db(md),an(2)
      dimension stl(nd),dl(md),avfo(10),avfc(10),avR(10),ph(nd,2)
      open(1,file='sas.dat')
      read(1,*)a,b,ce,al,be,ga
      read(1,*)NSYM
      do 16 i=1,nsym
      READ(1,17) (tx(i,j),(nt(i,(k-3+3*j)),k=1,3),j=1,3)
   17 format(3(f10.4,3i5))
   16 continue
      read(1,*)nel
      do 5 i=1,nel
      read(1,*) (c(j,i),j=1,9)
      ff=c(1,i)+c(3,i)+c(5,i)+c(7,i)+c(9,i)
    5 write(*,*)i,ff
      close(1)
      open(1,file='XYZ')
      read(1,*)scale
      do 1 i=1,mp
      g(i)=1.
    1 read(1,15,end=2)x(1,i),y(1,i),z(1,i),l(i),beta(i,1)
   15 format(3f10.6,i3,f10.6,i4)
    2 NAT=I-1
      close(1)
      write(*,23)
   23 format(' Binary data file name : ',$)
      read(*,21)an
   21 format(2a4)
      open(1,file=an,form='unformatted')
      do 3 i=1,25000
      read(1,end=4)ih(i),ik(i),il(i),ns,ne,F(i,1),F(i,2)
      stl(i)=0.0001*ns
    3 smax=max(smax,stl(i))
    4 NREF=I-1
      close(1)
      s3=(smax**3)/9.99
      PI=acos(-1.0)
      PI2=2*PI
      write(*,*)NREF,NAT,NSYM
      write(*,31)
   31 format(' MCYC, keep/omit (0,1) weak data, FACT ',$)
      read(*,*)MCYC,NSKIP,FACT
      if(FACT.le.0) FACT=1.
      t0=secnds(0.0)
      DO 14 NCYC=1,MCYC
      do 27 i=1,NAT
      do 27 j=2,NSYM
      x(j,i)=nt(j,1)*x(1,i) +nt(j,2)*y(1,i) +nt(j,3)*z(1,i) +tx(j,1)
      y(j,i)=nt(j,4)*x(1,i) +nt(j,5)*y(1,i) +nt(j,6)*z(1,i) +tx(j,2)
   27 z(j,i)=nt(j,7)*x(1,i) +nt(j,8)*y(1,i) +nt(j,9)*z(1,i) +tx(j,3)
      a1=0
      aa1=0
      ab=0
      a2=0
      aa2=0
      b1=0
      bb1=0
      ba=0
      b2=0
      bb2=0
C     LOOP OVER DATA TO BUILD MATRIX
      do 6 i=1,NREF
      astl=1.+(stl(i)**3)/s3
      jstl=astl
      hh=ih(i)*PI2
      hk=ik(i)*PI2
      hl=il(i)*PI2
      ss=-stl(i)*stl(i)
      do 7 j=1,nel
      d = c(1,j)*exp(ss*c(2,j)) +c(3,j)*exp(ss*c(4,j)) 
    7 fj(j) = d +c(5,j)*exp(ss*c(6,j)) +c(7,j)*exp(ss*c(8,j)) +c(9,j)
      ahs=0.
      bhs=0.
      jj=0
      do 8 j=1,NAT
      ah=0
      bh=0
      ff=fj(l(j))*exp(ss*beta(j,1))
      xa=0
      xb=0
      ya=0
      yb=0
      za=0
      zb=0
      do 9 k=1,NSYM
      arg = hh*x(k,j) + hk*y(k,j) + hl*z(k,j)
      ap=ff*cos(arg)
      bp=ff*sin(arg)
      xa=xa-bp*nt(k,1)
      xb=xb+ap*nt(k,1)
      ya=ya-bp*nt(k,5)
      yb=yb+ap*nt(k,5)
      za=za-bp*nt(k,9)
      zb=zb+ap*nt(k,9)
      ah=ah+ap
    9 bh=bh+bp
      da(jj+1)=xa*hh
      db(jj+1)=xb*hh
      da(jj+2)=ya*hk
      db(jj+2)=yb*hk
      da(jj+3)=za*hl
      db(jj+3)=zb*hl
      da(jj+4)=ss*ah
      db(jj+4)=ss*bh
      jj=jj+4
      ahs=ahs+ah
      bhs=bhs+bh
    8 continue
      jj=jj+1
      da(jj)=ahs
      db(jj)=bhs
      NAT9=jj
      Fc=sqrt(ahs*ahs+bhs*bhs)
      F(i,3)=Fc
      phi=atan2(bhs,ahs)
      a=F(i,1)*cos(phi)/scale
      b=F(i,1)*sin(phi)/scale 
      dah=a-Fc*cos(phi)
      dbh=b-Fc*sin(phi)
      wt=1.
      t1=secnds(0.0)
      do 11 j=1,NAT9
      w(j)=w(j)+wt*(da(j)*dah +db(j)*dbh)
      do 11 k=j,NAT9
   11 e(k,j)=e(k,j)+wt*(da(k)*da(j)+db(k)*db(j))
      t3=t3+secnds(t1)/60.
   30 a1=a1+a
      aa1=aa1+ahs
      ab=ab+ahs*a + bhs*b
      a2=a2+a*a +b*b
      aa2= aa2+ahs*ahs+bhs*bhs
      if(NCYC.ne.MCYC) goto 6
      np = NINT(1000.*phi)
      ns=10000*stl(i)
      Fo=F(i,1)/scale
      write(7)ih(i),ik(i),il(i),ns,ne,Fo,F(i,2),fc,np
   10 format(3i4,4f8.2)
      num(jstl)=num(jstl)+1
      avfo(jstl)=avfo(jstl)+F(i,1)/scale
      avfc(jstl)=avfc(jstl)+F(i,3)
      avR(jstl)=avR(jstl)+abs(F(i,1)/scale -F(i,3))
      dw(jstl)=dw(jstl)+F(i,2)
    6 continue
      do 29 i=1,NAT9
      do 29 j=i+1,NAT9
   29 e(i,j)=e(j,i)
      write(*,*)NAT9,NAT9
      t1=secnds(0.0)
      CALL MATINV(NAT9,NAT9,E,DET)
      t2=t2+secnds(t1)/60.
      amax=0.
      av=0
      do 13 i=1,NAT9
      a=0
      do 12 j=1,NAT9
      a=a+e(j,i)*w(j)
      if(i.eq.NAT9) w(j)=0.
   12 e(j,i)=0
      b1=b1+dl(i)
      bb1=bb1+dl(i)**2
      b2=b2+a
      bb2=bb2+a*a
      ba=ba+a*dl(i)
      dl(i)=a*FACT
      av=av+abs(a)
      amax=max(amax,abs(a))
   13 if(amax.eq.abs(a)) jm=i
      j=0
      do 19 i=1,NAT
      x(1,i)=x(1,i)+dl(j+1)
      y(1,i)=y(1,i)+dl(j+2)
      z(1,i)=z(1,i)+dl(j+3)
      if(abs(dl(j+4)).gt.10) dl(j+4)=dl(j+4)*10./abs(dl(j+4))
      if(abs(dl(j+4)).gt.abs(0.5*beta(i,1)))  dl(j+4)=0.5*dl(j+4)
      beta(i,1)=beta(i,1)+dl(j+4)
      if(beta(i,1).le.1.0) beta(i,1)=1.0
      j=j+4
   19 continue
      scale=scale+0.2*dl(NAT9)
      if(ncyc.ne.MCYC) goto 28
      write(1,15)scale
      do 20 i=1,NAT
   20 write(1,15)x(1,i),y(1,i),z(1,i),l(i),beta(i,1),i
   28 av=av/NAT
      a1=a1/nref
      aa1=aa1/nref
      cc=(ab/nref-a1*aa1)/sqrt((aa2/nref-aa1*aa1)*(a2/nref-a1*a1))
      b1=b1/NAT9
      b2=b2/NAT9
      c2=(ba/NAT9-b1*b2)/sqrt((bb1/NAT9-b1*b1)*(bb2/NAT9-b2*b2))
      if(NCYC.eq.1) c2=0
C     FACT=FACT*(1.+ 0.2*c2)
      write(*,18)NAT,NAT9,cc,DET,av,amax,scale,jm,c2,fact
   18 format(2i5,f8.5,E12.5,3f8.4,i4,2f8.5)
   14 CONTINUE
      write(*,26)
   26 format(' Shell average statistics'/'  NUM  <Fobs>  <Fcal>    <R>   
     $   GOOF')
      do 24 i=1,10
      fo=avfo(i)/num(i)
      fc=avfc(i)/num(i)
      fos=fos+avfo(i)
      fcs=fcs+avfc(i)
      R=avR(i)/avfo(i)
      GOOF=avR(i)/dw(i)
   24 write(*,25)num(i),fo,fc,R,GOOF
   25 format(i5,4f8.3)
      s=fos/fcs
      time=secnds(t0)/60.
      write(*,22)s,time,t2,t3
   22 format(4f8.4)
      end
      SUBROUTINE MATINV (N,M,A,D)
C
C     INVERT MATRIX BY GAUSS-JORDAN ELIMINATION, CALCULATE ITS DETERMINANT.
C
C     A - INPUT MATRIX, WHICH IS REPLACED BY ITS INVERSE
C     N - ORDER OF MATRIX
C     M - PHYSICAL SIZE OF ARRAY A, M .GE. N
C     D - DETERMINANT OF THE INPUT MATRIX
C
C     RETURNS D = 0 TO FLAG A SINGULAR MATRIX.
C
C     MATRIX IS SCALED TO AVOID OVERFLOW OR UNDERFLOW 
C
C     IF N EXCEEDS MAXN, THE WORK VECTORS Q, IK, AND JK MUST BE RE-
C     DIMENSIONED TO A LENGTH OF N.
C
      PARAMETER (MAXN=2400)
      DIMENSION A(MAXN,MAXN),Q(MAXN),IK(MAXN),JK(MAXN)
      REAL*4 Q,AMAX,T,A,D
C
C     STORE RECIPROCAL SQUARE-ROOT DIAGONAL MAGNITUDES FOR SCALING.
C
      DO 10 I=1,N
      T=DSQRT(DABS(A(I,I)))
      IF (T.NE.0) THEN
        Q(I)=1/T
      ELSE
        D=0
        RETURN
      END IF
   10 CONTINUE
C
C     SCALE MATRIX.
C
      DO 11 I=1,N
      DO 11 J=1,N
      A(I,J)=A(I,J)*Q(I)*Q(J)
   11 CONTINUE
C
C     INVERT SCALED MATRIX AND CALCULATE ITS DETERMINANT.
C
      D=1
      DO 100 K=1,N
C
C     FIND LARGEST ELEMENT A(I,J) IN REST OF MATRIX.
C
      AMAX=0
      DO 30 I=K,N
      DO 30 J=K,N
      IF (DABS(A(I,J)).GT.DABS(AMAX)) THEN
        AMAX=A(I,J)
        IK(K)=I
        JK(K)=J
      END IF
   30 CONTINUE
      D=D*AMAX
      IF (D.EQ.0) RETURN
C
C     INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN A(K,K).
C
      I=IK(K)
      IF (I.GT.K) THEN
        DO 50 J=1,N
        T=A(K,J)
        A(K,J)=A(I,J)
        A(I,J)=-T
   50   CONTINUE
      END IF
      J=JK(K)
      IF (J.GT.K) THEN
        DO 60 I=1,N
        T=A(I,K)
        A(I,K)=A(I,J)
        A(I,J)=-T
   60   CONTINUE
      END IF
C
C     ACCUMULATE ELEMENTS OF INVERSE MATRIX.
C
      DO 70 I=1,N
      IF (I.NE.K) A(I,K)=-A(I,K)/AMAX
   70 CONTINUE
      DO 80 I=1,N
      DO 80 J=1,N
      IF (I.NE.K.AND.J.NE.K) A(I,J)=A(I,J)+A(I,K)*A(K,J)
   80 CONTINUE
      DO 90 J=1,N
      IF (J.NE.K) A(K,J)=+A(K,J)/AMAX
   90 CONTINUE
      A(K,K)=1/AMAX
  100 CONTINUE
C
C     RESTORE ORDERING OF MATRIX.
C
      DO 130 K=N,1,-1
      J=IK(K)
      IF (J.GT.K) THEN
        DO 110 I=1,N
        T=A(I,K)
        A(I,K)=-A(I,J)
        A(I,J)=T
  110   CONTINUE
      END IF
      I=JK(K)
      IF (I.GT.K) THEN
        DO 120 J=1,N
        T=A(K,J)
        A(K,J)=-A(I,J)
        A(I,J)=T
  120   CONTINUE
      END IF
  130 CONTINUE
C
C     SCALE INVERSE MATRIX.
C
      DO 140 I=1,N
      DO 140 J=1,N
      A(I,J)=A(I,J)*Q(I)*Q(J)
  140 CONTINUE
C
C     SCALE DETERMINANT, IF POSSIBLE.
C
      T=0
      DO 150 I=1,N
      T=T+DLOG(Q(I))
  150 CONTINUE
      T=DLOG(DABS(D))-2*T
C
C     TEST AGAINST RANGE OF MACHINE ALLOWED MAGNITUDES.
C
C     EIGHT-BIT EXPONENT HAS MAXIMUM MAGNITUDE 2**7 = 128.
C     XMIN = 2**(-128)     = 0.29E-38      LOG(XMIN) = -88.7
C     XMAX = (2**(+128))/2 = 1.70E+38      LOG(XMAX) = +88.0
C
      IF (-87.LT.T.AND.T.LT.+87) D=(D/DABS(D))*EXP(T)
      RETURN
      END

        
