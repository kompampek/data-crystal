      PROGRAM RHOMOD
      parameter (mp=50000,NPAR=128,NPA=NPAR-1)
      CHARACTER AN*16
      REAL*8 rho
      integer*2 ih(mp),ik(mp),il(mp),np
      dimension e(mp),p(mp),q(mp),rho(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ rho,e,p,NREF,ih,ik,il,mh,mk,ml,mmk,mml,MX,MY,MZ
      write(*,1)
    1 format(' ENTER FILE NAME : ',$)
      read(*,2)AN
    2 format(a)
      write(*,3)
    3 format(' thres = ',$)  
      read(*,*)thres
      CALL CPU_TIME(t0)
      open(1,file=AN)
      do 4 i=1,50000
      read(1,*,end=5)ih(i),ik(i),il(i),ns,ne,fo,sf,fc,nph
      mh=max(ih(i),mh)
      mk=max(ik(i),mk)
      ml=max(il(i),ml)
      mmk=min(ik(i),mmk)
      mml=min(il(i),mml)
      e(i)=0.001*ne
    4 p(i)=0.001*nph
    5 close(1)
      NREF=i-1
      var=0.
      CALL FORWARD(var,thres)
      write(*,6)Nref,var
    6 format(i6,f10.5)
      do 7 i=1,2
      CALL FORWARD(var,thres)
      CALL BACK
    7 CONTINUE  
      CALL CPU_TIME(t1)
      T1=(T1-T0)/60.
      write(*,*)t1
      END
C      
      SUBROUTINE FORWARD(var,thres)
      parameter (NPAR=128,NPA=NPAR-1,NP=NPAR/2,MP=50000)
      byte nh(-NP:NP,-NP:NP),nt(16,9),HH,HK,HL
      real*8 arg,px,data,r,pa,pb,fa,fb,cs,sn,ff,ph,fsa,fsb
      integer*2 ih(mp),ik(mp),il(mp)
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp),cs(0:NPA),sn(0:NPA)
      dimension data(2*NPA),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),tx(16,3)
      dimension fa(0:NP,-NP:NP,-NP:NP),fb(0:NP,-NP:NP,-NP:NP)
      COMMON /HKL/ r,e,p,NREF,ih,ik,il,mh,mk,ml,mmk,mml,MX,MY,MZ
      if(var.ne.0) goto 12
      do 13 i=0,NP
      do 13 j=-NP,NP
      nh(i,j)=0
      nh(-i,j)=0
      do 13 k=-NP,NP
      fa(i,j,k)=0.
   13 fb(i,j,k)=0.  
      px=-1.0
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MX,MY,MZ
      close(3)
      pi2=2.0*acos(-1.0)
      px=2.0*DACOS(px)/MX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MX-1
      arg=-i*px
      sn(i)=sin(arg)
   20 cs(i)=cos(arg)
C     read reflection data
   12 continue
      MX1=MX-1
      MY1=MY-1
      MZ1=MZ-1
      MX2=2*MX
      MY2=2*MY
      MZ2=2*MZ
      do 11 i=1,NREF
      jh=ih(i)
      jk=ik(i)
      jl=il(i)
      ff=2.*e(i)
      nh(jl,jk)=1
      fa(jh,jk,jl)=ff*cos(p(i))
   19 fb(jh,jk,jl)=ff*sin(p(i))
   11 continue
      if(var.ne.0) goto 24
      do 23 jh =  0,mh
      do 23 jk = mmk,mk
      do 23 jl = mml,ml
      a = fa(jh,jk,jl)
      b = fb(jh,jk,jl)
      var = var + a*a +b*b
   23 CONTINUE
      var=SQRT(var/2.)
      write(*,*)var,99
      RETURN
   24 continue
      back=thres*var
      sr2=0.
      aminr=99.
      amaxr=-99.
      write(*,*)sr2
C     calculate FOURIER
      do 1 ix=0,MX1
      do 3 ikk = mmk,mk
      do 4 iz=1,MZ2
    4 data(iz)=0.
      do 5 ill = mml,ml
      if(nh(ill,ikk).eq.0) goto 5
      fsa=0.
      fsb=0.
      j1=-ix
      do 2 ihh=0,mh
      j1=j1+ix
      j1=iand(j1,MX1)
      fsa=fsa+fa(ihh,ikk,ill)*cs(j1)-fb(ihh,ikk,ill)*sn(j1)
    2 fsb=fsb+fa(ihh,ikk,ill)*sn(j1)+fb(ihh,ikk,ill)*cs(j1)
      j1=iand(ill,MZ1)
      j1=j1+j1+1
      data(j1)  =fsa
      data(j1+1)=fsb
    5 continue
      CALL FFT(data,MZ,1)
      jk=iand(ikk,MY1)
      izz=-1
      do 6 iz=0,MZ1
      izz=izz+2
      pa(jk,iz)=data(izz)
    6 pb(jk,iz)=data(izz+1)
    3 continue
      do 7 iz=0,MZ1
      do 8 ikk=1,MY2
    8 data(ikk)=0.
      do 9 ikk = mmk,mk
      jk=iand(ikk,MY1)
      j1=jk+jk+1
      data(j1)=pa(jk,iz)
    9 data(j1+1)=pb(jk,iz)
      CALL FFT(data,MY,1)
      iyy=-1
      do 10 iy=0,MY1
      iyy=iyy+2
      rho=data(iyy)
      aminr=min(aminr,rho)
      amaxr=max(amaxr,rho)
      sr2=sr2+rho*rho
C      if(rho.lt.back) rho=-rho
   10 r(iy,iz,ix)=rho
    7 continue
    1 continue
      write(77,79)r
   79 format(10f8.1)   
      sr3=sr2/(MX*MY*MZ)
      write(*,*)aminr,amaxr,sr2,sr3,mh,mk,ml,mmk,mml
C     end of FOURIER loops
      RETURN
      END
C     
      SUBROUTINE FFT(DATA,NN,IT)
      PARAMETER (NPAR=128)
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
      parameter (NPAR=128,NPA=NPAR-1,mp=50000)
      real*8 arg,px,a1,a2,aL1,aL2,aa,DATA,CC,cs,sn
      real*8 r,pa,pb,GFA,GFB,fsa,fsb,fas,fbs,G0
      integer*4 hh,hk,hl
      integer*2 ih(mp),ik(mp),il(mp),nt(16,9)
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp),tx(16,3),cs(0:NPA)
      dimension data(2*NPA),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),sn(0:NPA)
      dimension GFA(0:NPA,0:NPA,0:NPA),GFB(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ r,e,p,NREF,ih,ik,il,mh,mk,ml,mmk,mml,MX,MY,MZ
      px=-1.0
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MPAX,MPAY,MPAZ
      fact=2/(MPAX*MPAY*MPAZ)
      close(3)
      pi2=2.0*acos(-1.0)
      px=2.0*DACOS(px)/MPAY
C     initialize exp(-twopi*x/MPAR), remember positive sign ??
      do 13 i=0,MPAY-1
      arg = i*px
      cs(i)=dcos(arg)*fact
   13 sn(i)=dsin(arg)*fact
      MPAX1=MPAX-1
      MPAY1=MPAY-1
      MPAZ1=MPAZ-1
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
      j1=j1+jkk
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
      do 7 jhh=0,MPAX1
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
C      end of FOURIER loops
      aa=0
      a1=0
      a2=0
      aL1=0
      aL2=0
      sum=0.
      do 14 i=1,10
      aL = e(i)*cos(p(i))
      bE = e(i)*sin(p(i))
      i1=iand(il(i),mpaz1)
      i2=iand(ih(i),mpax1)
      i3=iand(ik(i),mpay1)
      a= GFA(i1,i2,i3)
      b= GFB(i1,i2,i3)
      ph = atan2(b,a)
      del = cos(p(i)-ph)
      del = acos(del)
      sum=sum+del
      a1=a1+a
      a2=a2+a*a+b*b
      aL1=aL1+aL
      aL2=aL2+e(i)*e(i)
      aa=aa +a*aL +b*be
      p(I)=ph
      write(*,*)i,aL,bE,a,b
   14 CONTINUE
      sum=57.296*sum/NREF
      a1=a1/NREF
      aL1=aL1/NREF
      a2=dsqrt((a2/NREF)*(aL2/NREF))
      CC = (aa/NREF-a1*aL1)/a2
      write(*,4)sum,cc
      write(21,4)sum,cc
    4 format(f6.1,f7.3,$)  
      RETURN
      END
