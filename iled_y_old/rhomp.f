      PROGRAM RHOMOD 
      parameter (mp=10000,NPAR=129,NPA=NPAR-1)
      REAL*8 rho
      integer*2 ih(mp),ik(mp),il(mp),np
      dimension e(mp),p(mp),q(mp),rho(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ rho,e,p,NREF,ih,ik,il,mh,mk,ml,mmk,mml
      write(*,3)
    3 format(' THRES value = ',$)  
      read(*,*)THRES
      CALL CPU_TIME(t0)
      open(1,file='iled.E')
      do 4 i=1,840
      read(1,*,end=5)ih(i),ik(i),il(i),ns,ne,fo,sf,fc,np
      e(i)=0.001*ne
      p(i)=0.001*np
    4 q(i)=0.001*np
    5 close(1)
      NREF=i-1
      var=0.
      CALL FORWARD(THRES,var)
      write(*,9)Nref,var
    9 format(i6,f10.5)
      CALL FORWARD(THRES,var)
      CALL BACK
      CALL CPU_TIME(t1)
      T1=(T1-T0)/60.
      write(*,*)t1
      END
C      
      SUBROUTINE FORWARD(thres,var)
      parameter (NPAR=129,NPA=NPAR-1,NP=NPAR/2,MP=10000)
      byte nh(-NP:NP,-NP:NP),nt(16,9),HH,HK,HL,Iset(mp)
      real*8 arg,px,data,r,pa,pb,fa,fb,cs,sn,ff,ph,fsa,fsb
      integer*2 ih(mp),ik(mp),il(mp),Jset(mp)
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp),cs(0:NPA),sn(0:NPA)
      dimension data(2*NPA),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),tx(16,3)
      dimension fa(0:NP,-NP:NP,-NP:NP),fb(0:NP,-NP:NP,-NP:NP)
      COMMON /HKL/ r,e,p,NREF,ih,ik,il,mh,mk,ml,mmk,mml
      DATA ISET/mp*0/
      open(3,file='BASIS')
      read(3,22)lp
      read(3,22)(jset(i),i=1,lp)
   22 format(20i5)
      CLOSE(3)  
      do 21 i=1,lp
      j=JSET(i)
      if(var.eq.0) write(*,*)J
   21 ISET(j)=1
      lp=NP
      do 12 k=  0,lp
      do 12 j=-lp,lp
      do 12 i=-lp,lp
      fb(k,j,i)=0.
   12 fa(k,j,i)=0.
      mh=0
      mk=0
      ml=0
      mmk=0
      mml=0
      px=-1.0
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MX,MY,MZ
      MX1=MX-1
      MY1=MY-1
      MZ1=MZ-1
      MX2=2*MX
      MY2=2*MY
      MZ2=2*MZ
      read(3,*)nsym,VOL
      do 16 i=1,nsym
      read(3,17)tx(i,1),HH,HK,HL,tx(i,2),KH,KK,KL,tx(i,3),LH,LK,LL
   17 format(3(f10.4,3i5))
      IDET = HH*KK*LL+HK*KL*LH+KH*LK*HL-LH*KK*HL-KH*HK*LL-HH*LK*KL
      nt(i,1) = (KK*LL-LK*KL)/IDET
      nt(i,2) =-(KH*LL-LH*KL)/IDET
      nt(i,3) = (KH*LK-KK*LH)/IDET
      nt(i,4) =-(HK*LL-HL*LK)/IDET
      nt(i,5) = (HH*LL-LH*HL)/IDET
      nt(i,6) =-(HH*LK-LH*HK)/IDET
      nt(i,7) = (HK*KL-KK*HL)/IDET
      nt(i,8) =-(HH*KL-HL*KH)/IDET
   16 nt(i,9) = (HH*KK-HK*KH)/IDET
      close(3)
      pi2=2.0*acos(-1.0)
      px=2.0*DACOS(px)/MX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MX1
      arg=-i*px
      sn(i)=dsin(arg)
   20 cs(i)=dcos(arg)
C     read reflection data
      do 11 i=1,NREF
      if(ISET(i).eq.0) goto 11
      ihh=ih(i)
      ikk=ik(i)
      ill=il(i)
      ff=2.*e(i)/VOL
      do 19 j=1,nsym
      jh = ihh*nt(j,1) + ikk*nt(j,2) + ill*nt(j,3)
      jk = ihh*nt(j,4) + ikk*nt(j,5) + ill*nt(j,6)
      jl = ihh*nt(j,7) + ikk*nt(j,8) + ill*nt(j,9)
      ph2 =(ihh*tx(j,1) + ikk*tx(j,2) + ill*tx(j,3))*PI2 + p(i)
      ph = mod(ph2,PI2)
      if(jh) 13,14,15
   14 if(jk) 13,18,15
   18 if(jl) 13,15,15
   13 jh=-jh
      jk=-jk
      jl=-jl
      ph=-ph
   15 nh(jl,jk)=1
      mh=max(jh,mh)
      mk=max(jk,mk)
      ml=max(jl,ml)
      mmk=min(jk,mmk)
      mml=min(jl,mml)
      fa(jh,jk,jl)=ff*dcos(ph)
   19 fb(jh,jk,jl)=ff*dsin(ph)
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
      write(*,*)var
      RETURN
   24 continue
C     calculate FOURIER
      rmax=-999.
      rmin=999.
      rss = 0
      IXYZ=MX*MY*MZ
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
      rmax=max(rmax,rho)
      rmin=min(rmin,rho)
      rss=rss+rho*rho
      if(rho.le.0.0) rho=0.
      sig=rho/var
      xx=exp(sig)
      yy=1./xx
      wt=(xx-yy)/(xx+yy)
   10 r(iy,iz,ix)=rho*wt
    7 continue
    1 continue
      rss=sqrt(rss/IXYZ)
      write(*,*)99,rss,rmax,rmin
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
      parameter (NPAR=129,NPA=NPAR-1,mp=10000)
      real*8 arg,px,a1,a2,aL1,aL2,aa,DATA,CC,cs,sn
      real*8 r,pa,pb,GFA,GFB,fsa,fsb,fas,fbs,G0
      integer*4 hh,hk,hl
      integer*2 ih(mp),ik(mp),il(mp),nt(16,9)
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp),tx(16,3),cs(0:NPA)
      dimension data(2*NPA),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),sn(0:NPA)
      dimension GFA(0:NPA,0:NPA,0:NPA),GFB(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ r,e,p,NREF,ih,ik,il,mh,mk,ml,mmk,mml
      px=-1.0
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MPAX,MPAY,MPAZ
      read(3,*)nsym,VOL
      fact=2*VOL/(MPAX*MPAY*MPAZ)
      do 16 i=1,nsym
      read(3,11)tx(i,1),HH,HK,HL,tx(i,2),KH,KK,KL,tx(i,3),LH,LK,LL
   11 format(3(f10.4,3i5))
      IDET = HH*KK*LL+HK*KL*LH+KH*LK*HL-LH*KK*HL-KH*HK*LL-HH*LK*KL
      nt(i,1) = (KK*LL-LK*KL)/IDET
      nt(i,2) =-(KH*LL-LH*KL)/IDET
      nt(i,3) = (KH*LK-KK*LH)/IDET
      nt(i,4) =-(HK*LL-HL*LK)/IDET
      nt(i,5) = (HH*LL-LH*HL)/IDET
      nt(i,6) =-(HH*LK-LH*HK)/IDET
      nt(i,7) = (HK*KL-KK*HL)/IDET
      nt(i,8) =-(HH*KL-HL*KH)/IDET
   16 nt(i,9) = (HH*KK-HK*KH)/IDET
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
      aa=0.
      a1=0.
      a2=0.
      aL1=0.
      aL2=0.
      sum=0.
      do 14 i=1,NREF
      i1=iand(il(i),mpaz1)
      i2=iand(ih(i),mpax1)
      i3=iand(ik(i),mpay1)
      a= GFA(i1,i2,i3)
      b= GFB(i1,i2,i3)
      ramp=sqrt(a*a+b*b)/2.
      ph = atan2(b,a)
      del = acos(cos(p(i)-ph))
      sum=sum+del
      a1=a1 +ramp
      a2=a2 +ramp*ramp
      aL1=aL1 +e(i)
      aL2=aL2 +e(i)*e(i)
      aa=aa +ramp*e(i)
      write(42,422)i,e(i),ramp,del
  422 format(i5,3f8.3)   
   14 CONTINUE
      sum=57.296*sum/NREF
      a1=a1/NREF
      aL1=aL1/NREF
      a2=a2/NREF
      aL2=aL2/NREF
      aa=aa/NREF
      a3=dsqrt((a2-a1*a1)*(aL2-aL1*aL1))
      write(*,*)a1,a2,AL1,AL2,aa,a3
      CC = (aa-a1*aL1)/a3
      write(*,4)sum,cc
    4 format(f6.1,f7.3)  
      RETURN
      END
