      parameter (mp=4000,na=128)
      character fil*6
      integer*2 ih(mp),ik(mp),il(mp)
      integer*4 iseed
      real*8 rho
      dimension p(mp),e(mp),rho(0:na,0:na,0:na)
      dimension cs(-32768:32768),sn(-32768:32768)
      COMMON /SIGNS/ cs,sn
      COMMON /DEN/ rho,tlow,thi,MPAX,MPAY,MPAZ
      COMMON /REF/ MREF,ih,ik,il,e,p,NSYM,mh,mk,ml,mmk,mml
      data mh,mk,ml,mmk,mml/5*0/
      open(2,file='FFT.dat',status='old',form='formatted')
      read(2,*)MPAX,MPAY,MPAZ
      close(unit=2)
      pi = acos(-1.0)
      pi2=2.*pi
      write(*,3)
    3 format(' DATA FILE NAME ',$)
      read(*,4)fil
    4 format(A)
      CALL CPU_TIME(t0)
C     READ E DATA FILE
      open(unit=1,file=fil,status='old',form='formatted')
      esum=0.
      do 1 i=1,mp
      read(1,*,end=5)ih(i),ik(i),il(i),ns,ne,fo,sf,fc,nar
      p(i)=0.001*nar
      if(p(i).gt.PI) p(i)=p(i)-PI2
      e(i)=0.001*ne
    1 esum=esum+e(i)**2
    5 MREF=i-1
      close(unit=1)
      sig=sqrt(esum/2.)
      write(*,*)sig,88
      do 6 i=1,10
      CALL FORWARD
    6 CALL BACK
      CALL CPU_TIME(t1)
      time=t1-t0
      write(*,*)time
      END
C
      SUBROUTINE SFCLC(NAT,IRAN)
      parameter (mp=4000,npa=128)
      INTEGER*4 ISEED
      real*8 rand
      integer*2 jh(mp),jk(mp),jl(mp)
      integer*2 ih,ik,il,nx(499),ny(499),nz(499),narg
      dimension e(mp),p(mp)
      dimension cs(-32768:32768),sn(-32768:32768)
      COMMON /SIGNS/ cs,sn
      COMMON /REF/ NREF,jh,jk,jl,e,p,NSYM,mh,mk,ml,mmk,mml
      ISEED=IRAN
      call srand(ISEED)
      do 1 i=1,NAT
      x=rand()
      y=rand()
      z=rand()
      nx(i)=65536.*x
      ny(i)=65536.*y
      nz(i)=65536.*z
    1 CONTINUE
      PI=ACOS(-1.0)
      arg=pi/32768.
      do 2 i=-32768,32768
      ar=i*arg
      sn(i)=sin(ar)
    2 cs(i)=cos(ar)
      do 4 i=1,NREF
      ih=jh(i)
      ik=jk(i)
      il=jl(i)
      a=0.
      b=0.
      do 3 j=1,NAT
      narg = ih*nx(j) + ik*ny(j) + il*nz(j)
      a = a + cs(narg)
    3 b = b + sn(narg)
c      p(i)=atan2(b,a)
    4 CONTINUE
      RETURN
      END
C
C      
      SUBROUTINE FORWARD
      parameter (NPAR=129,NPA=NPAR-1,NP=NPAR/2,MP=4000)
      byte nh(-NP:NP,-NP:NP),HH,HK,HL
      real*8 arg,pix,data,r,pa,pb,fa,fb,cs,sn,ff,ph,fsa,fsb
      integer*2 ih(mp),ik(mp),il(mp),jh,jk,jl
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp)
      dimension cs(0:NPA),sn(0:NPA)
      dimension data(2*NPA),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA)
      dimension fa(0:NP,-NP:NP,-NP:NP),fb(0:NP,-NP:NP,-NP:NP)
      COMMON /DEN/ r,tlow,thi,MPAX,MPAY,MPAZ
      COMMON /REF/ NREF,ih,ik,il,e,p,NSYM,mh,mk,ml,mmk,mml
      pix=-1.0
      MPAX1=MPAX-1
      MPAY1=MPAY-1
      MPAZ1=MPAZ-1
      MPAX2=2*MPAX
      MPAY2=2*MPAY
      MPAZ2=2*MPAZ
      PI2=2.0*acos(-1.0)
      pix=2.0*DACOS(pix)/MPAX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPAX1
      arg=-i*pix
      sn(i)=dsin(arg)
   20 cs(i)=dcos(arg)
  177 CONTINUE
C     read reflection data
      do 11 i=1,NREF
      jh=ih(i)
      jk=ik(i)
      jl=il(i)
      ff=e(i)
      ph=p(i)
      nh(jl,jk)=1
      mh=max(jh,mh)
      mk=max(jk,mk)
      ml=max(jl,ml)
      mmk=min(jk,mmk)
      mml=min(jl,mml)
      fa(jh,jk,jl)=ff*dcos(ph)
      fb(jh,jk,jl)=ff*dsin(ph)
   11 continue
      if((mh+mh).ge.MPAX.or.(mk+mk).ge.MPAY.or.(ml+ml).ge.MPAZ) RETURN
      CONTINUE
      VAR=0.
      R3=0.
C     calculate FOURIER
      do 1 ix=0,MPAX1
      do 3 ikk = mmk,mk
      do 4 iz=1,MPAZ2
    4 data(iz)=0.
      do 5 ill = mml,ml
      if(nh(ill,ikk).eq.0) goto 5
      fsa=0.
      fsb=0.
      j1=-ix
      do 2 ihh=0,mh
      j1=(j1+ix)
      j1= iand(j1,MPAX1)
      fsa=fsa+fa(ihh,ikk,ill)*cs(j1)-fb(ihh,ikk,ill)*sn(j1)
    2 fsb=fsb+fa(ihh,ikk,ill)*sn(j1)+fb(ihh,ikk,ill)*cs(j1)
      j1=iand(ill,MPAZ1)
      j1=j1+j1+1
      data(j1)  =fsa
      data(j1+1)=fsb
    5 continue
      CALL FFT(data,MPAZ,1)
      jk=iand(ikk,MPAY1)
      izz=-1
      do 6 iz=0,MPAZ1
      izz=izz+2
      pa(jk,iz)=data(izz)
    6 pb(jk,iz)=data(izz+1)
    3 continue
      do 7 iz=0,MPAZ1
      do 8 ikk=1,MPAY2
    8 data(ikk)=0.
      do 9 ikk = mmk,mk
      jk=iand(ikk,MPAY1)
      j1=jk+jk+1
      data(j1)=pa(jk,iz)
    9 data(j1+1)=pb(jk,iz)
      CALL FFT(data,MPAY,1)
      iyy=-1
      do 10 iy=0,MPAY1
      iyy=iyy+2
      a=data(iyy)
      var=var+a*a
      r3=r3+a**3
   10 r(iy,iz,ix)=a
    7 continue
    1 continue
      var=var/(MPAX*MPAY*MPAZ)
      r3=r3/(MPAX*MPAY*MPAZ)
      sig=sqrt(var)
      write(*,*)sig,r3,99
C     end of FOURIER loops
      RETURN
      END
C     
      SUBROUTINE FFT(DATA,NN,IT)
      PARAMETER (NPAR=129)
C     NPAR MUST BE A POWER OF 2
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,TEMPR,TEMPI,PI2,DATA
      DIMENSION DATA(2*NPAR)
      N=2*NN
      J=1
      PI2 = -1.
      PI2 = IT*2.0*DACOS(PI2)
C     PERFORMS THE BIT REVERSAL
      DO 11 I=1,N,2
         IF(J.GT.I) THEN
         TEMPR=DATA(J)
         TEMPI=DATA(J+1)
         DATA(J)=DATA(I)
         DATA(J+1)=DATA(I+1)
         DATA(I)=TEMPR
         DATA(I+1)=TEMPI
         ENDIF
         M=N/2
    1 IF((M.GE.2).AND.(J.GT.M)) THEN
        J=J-M
        M=M/2
        GOTO 1
        ENDIF
        J=J+M
   11 CONTINUE
      MMAX=2
C      DO THIS LOOP LOG2(NN) TIMES
    2 IF(N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=-PI2/MMAX
        WPR= -2.*DSIN(0.5*THETA)**2
        WPI=DSIN(THETA)
        WR=1.
        WI=0.
      DO 13 M=1,MMAX,2
      DO 12 I=M,N,ISTEP
         J=I+MMAX
         J1=J+1
         I1=I+1
         TEMPR=WR*DATA(J) - WI*DATA(J1)
         TEMPI=WR*DATA(J1) + WI*DATA(J)
         DATA(J)=DATA(I)-TEMPR
         DATA(J1)=DATA(I1)-TEMPI
         DATA(I)=DATA(I)+TEMPR
         DATA(I1)=DATA(I1)+TEMPI
   12 CONTINUE
         WTEMP=WR
         WR=WR*WPR -WI*WPI +WR
         WI=WI*WPR +WTEMP*WPI +WI
   13 CONTINUE
       MMAX=ISTEP
       GOTO 2
       ENDIF
       RETURN
       END
C
      SUBROUTINE BACK
      parameter (NPAR=129,NPA=NPAR-1,mp=4000)
      real*8 arg,pix,DATA,cs,sn
      real*8 r,pa,pb,GFA,GFB,fsa,fsb,fas,fbs
      integer*2 ih(mp),ik(mp),il(mp),HH,HK,HL
      dimension r(0:NPA,0:NPA,0:NPA)
      dimension e(mp),p(mp),cs(0:NPA),sn(0:NPA)
      dimension data(2*NPA),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA)
      dimension GFA(0:NPA,0:NPA,0:NPA),GFB(0:NPA,0:NPA,0:NPA)
      COMMON /DEN/ r,tlow,thi,MPAX,MPAY,MPAZ
      COMMON /REF/NREF,ih,ik,il,e,p,NSYM,mh,mk,ml,mmk,mml
      pix=-1.0
      fact=2./(MPAX*MPAY*MPAZ)
      pi2=2.0*acos(-1.0)
      pix=2.0*DACOS(pix)/MPAY
C     initialize exp(-twopi*x/MPAR), remember positive sign ??
      do 13 i=0,MPAY-1
      arg = i*pix
      cs(i)=dcos(arg)*fact
   13 sn(i)=dsin(arg)*fact
      MPAX1=MPAX-1
      MPAY1=MPAY-1
      MPAZ1=MPAZ-1
      MPAX2=2*MPAX
      MPAY2=2*MPAY
      MPAZ2=2*MPAZ
C
      do 1 jkk= 0,MPAY1
      do 3 jz = 0,MPAZ1
      j0=-1
      do 5 jx = 0,MPAX1
      fsa=0.
      fsb=0
      j1=-jkk
      do 2 jy=0,MPAY1
      rho=r(jy,jz,jx)
      j1=(j1+jkk)
      j1=iand(j1,MPAY1)
      fsa=fsa+rho*cs(j1)
      fsb=fsb+rho*sn(j1)
    2 continue
      j0=j0+2
      data(j0)  =fsa
    5 data(j0+1)=fsb
      CALL FFT(data,MPAX,-1)
      izz=-1
      do 6 jhh=0,MPAX1
      izz=izz+2
      pa(jz,jhh)=data(izz)
    6 pb(jz,jhh)=data(izz+1)
    3 continue
      do 7 jhh=0,mh
      j1=-1
      do 9 jz =0,MPAZ1
      j1=j1+2
      data(j1)  =pa(jz,jhh)
    9 data(j1+1)=pb(jz,jhh)
      CALL FFT(data,MPAZ,-1)
      iyy=-1
      do 10 jll=0,MPAZ1
      iyy=iyy+2
      GFA(jll,jhh,jkk)=data(iyy)
   10 GFB(jll,jhh,jkk)=data(iyy+1)
    7 continue
    1 continue
C    End of FOURIER loops, calc X(I) and CC
      do 14 i=1,10,4
      ihi=iand(ih(i),mpax1)
      ki =iand(ik(i),mpay1)
      li =iand(il(i),mpaz1)
      a=GFA(li,ihi,ki)
      b=GFB(li,ihi,ki)
      fx=sqrt(a*a+b*b)
      px=atan2(b,a)
      write(*,141)i,e(i),p(i),fx,px
  141 format(i5,4f6.2)    
   14 CONTINUE
      RETURN
      END

