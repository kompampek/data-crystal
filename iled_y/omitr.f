C     LIMITS 10K unique data, 128**3 grid 
      character  AN*12
      REAL*8 CC,CC1,CC2
      parameter (mp=150,NPAR=128,NPA=NPAR-1)
      integer*2 ih(mp),ik(mp),il(mp),np
      dimension e(mp),q(mp),p(mp),rho(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ rho,NREF,ih,ik,il,e,p,mh,mk,ml,mmk,mml
      write(*,2)
    2 format(' INPUT {hk0} file name: ')  
      read(*,1)AN
      write(*,1)AN
    1 format(A)  
      open(1,file=AN)
      write(*,3)
    3 format(' INPUT CODE file name : ')  
      read(*,1)AN
      open(77,file=AN)
      read(77,*)NSETS,nbits
      write(*,11)
   11 format('NFIXED, #SIG_THRESHOLD  ',$)
      read(*,*)NFIX,THRES
      NREF=NFIX+nbits
      do 4 i=1,NREF
      read(1,*,end=5)ih(i),ik(i),il(i),e(i),q(i)
      p(i)=q(i)
    4 CONTINUE
    5 close(1)
      mref=i-1
      write(*,*)97,NFIX,NREF,mref,thres
      t0=secnds(0.0)
      va=-99.99
      CALL FORWARD(va)
      var=THRES*va
      CALL FORWARD(VAR)
      CALL BACK(VAR,CC,CC2)
      write(*,7)0,cc,0.0,0.0,cc2
      write(22,7)0,cc,0.0,0.0,cc2
    7 format(i9,f9.5,2f6.1,f9.5)
C    
      do 6 i=1,NSETS
      read(77,*)nrob
      nr=nrob
      do 8 j=NFIX+1,NREF
      k=iand(nr,1)
      nr=nr/2
      p(j)=q(j)+k*3.141592
    8 continue
      bs=999.
      jx=0
      do 9 i1=0,1
      do 9 i2=0,1
      do 9 i3=0,1
      do 9 i4=-1,1,2
      as=0.
      do 10 k=1,NREF
      ix= ih(k)*i1+ik(k)*i2+il(k)*i3
      ix=iand(ix,1)
      dp=cos(q(k)+i4*p(k)+ix*3.141592)
   10 as=as+acos(dp)
      as=57.296*as/NREF
      if(jx.eq.0) ang0=as
      jx=jx+1
    9 if(as.lt.bs) bs=as
      call forward(var)
      call back(VAR,CC,CC2)
      write(*,7) nrob,cc,ang0,bs,cc2
      write(22,7)nrob,cc,ang0,bs,cc2
    6 CONTINUE
      t1=secnds(t0)/60
      write(*,*)t1
      write(22,*)t1
      end
C      
      SUBROUTINE FORWARD(var)
      parameter (NPAR=128,NPA=NPAR-1,NP=NPAR/2,MP=150)
      byte nh(-NP:NP,-NP:NP),nt(16,9),HH,HK,HL
      real*8 arg,px,data
      integer*2 ih(mp),ik(mp),il(mp),jh,jk,jl
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp),cs(0:NPA),sn(0:NPA)
      dimension data(2*NPAR),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),tx(16,3)
      dimension fa(0:NP,-NP:NP,-NP:NP),fb(0:NP,-NP:NP,-NP:NP)
      COMMON /HKL/ r,NREF,ih,ik,il,e,p,mh,mk,ml,mmk,mml
      data mh,mk,ml,mmk,mml/5*0/
      px=-1.0
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MPAX,MPAY,MPAZ
      MPAX1=MPAX-1
      MPAY1=MPAY-1
      MPAZ1=MPAZ-1
      MPAX2=2*MPAX
      MPAY2=2*MPAY
      MPAZ2=2*MPAZ
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
      px=2.0*DACOS(px)/MPAX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPAX1
      arg=-i*px
      sn(i)=dsin(arg)
   20 cs(i)=dcos(arg)
C     read reflection data
      do 11 i=1,NREF
      ihh=ih(i)
      ikk=ik(i)
      ill=il(i)
      ff=2.*e(i)/VOL
      phase=p(i)
      do 19 j=1,nsym
      jh = ihh*nt(j,1) + ikk*nt(j,2) + ill*nt(j,3)
      jk = ihh*nt(j,4) + ikk*nt(j,5) + ill*nt(j,6)
      jl = ihh*nt(j,7) + ikk*nt(j,8) + ill*nt(j,9)
      ph =(ihh*tx(j,1) + ikk*tx(j,2) + ill*tx(j,3))*PI2
      ph = mod((ph+phase),PI2)
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
      mum=mum+1
      fa(jh,jk,jl)=ff*cos(ph)
   19 fb(jh,jk,jl)=ff*sin(ph)
   11 continue
      if((mh+mh).ge.MPAX.or.(mk+mk).ge.MPAY.or.(ml+ml).ge.MPAZ) RETURN
      if(var.ne.-99.99) goto 24
      var=0.
      do 23 jl = mml,ml
      do 23 jk = mmk,mk
      do 23 jh =   0,mh
      a = fa(jh,jk,jl)
      b = fb(jh,jk,jl)
   23 var = var + a*a+b*b
      var=SQRT(var/2.)
      write(*,*)var
      RETURN
   24 continue
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
      j1=j1+ix
      j1=iand(j1,MPAX1)
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
   10 r(iy,iz,ix)=a
    7 continue
    1 continue
C     end of FOURIER loops
      return
      end
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
      SUBROUTINE BACK(VAR,CC,CC2)
      parameter (NPAR=128,NPA=NPAR-1,mp=150 )
      real*8 arg,px,a1,a2,aL1,aL2,aa,DATA,CC,CC1,CC2
      integer*2 ih(mp),ik(mp),il(mp),nt(16,9),HH,HK,HL
      integer*2 hi(16),ki(16),li(16)
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp),tx(16,3),cs(0:NPA)
      dimension data(2*NPAR),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),sn(0:NPA)
      dimension GFA(0:NPA,0:NPA,0:NPA,2),GFB(0:NPA,0:NPA,0:NPA,2)
      COMMON /HKL/ r,NREF,ih,ik,il,e,p,mh,mk,ml,mmk,mml
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
      MPAX2=2*MPAX
      MPAY2=2*MPAY
      MPAZ2=2*MPAZ
C     jm=1 (modified,Xh), jm=2 (mask)
      do 8 jm=1,2
      do 1 jkk= 0,MPAY1
      do 3 jz = 0,MPAZ1
      j0=-1
      do 5 jx = 0,MPAX1
      fsa=0.
      fsb=0
      j1=-jkk
      do 2 jy=0,MPAY1
      rho=r(jy,jz,jx)
      jb = j1+jkk
      j1=iand(jb,MPAY1)
      if(rho.lt.0.0) rho=0.
      if(rho.eq.0.0) goto 2
      tex=1.0
      if(rho.lt.var) tex=rho/var
C      tex=tex*tex
      rho=rho*tex
      if(jm.eq.2) rho=tex
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
      GFA(jll,jhh,jkk,jm)=data(iyy)
   10 GFB(jll,jhh,jkk,jm)=data(iyy+1)
    7 continue
    1 continue
    8 CONTINUE
C     end of FOURIER loops
      aa=0
      a1=0
      a2=0
      aL1=0
      aL2=0
      do 14 i=1,NREF
      hh=ih(i)
      hk=ik(i)
      hl=il(i)
      it=iand(hh,1)+iand(hk,1)+iand(hl,1)
      nwt=0
      if(it.eq.0) nwt=1
      fas=0
      fbs=0
      IT=0
      do 12 inv=-1,1,2
      ph=inv*p(i)
      do 15 j=1,nsym
      KH = (hh*nt(j,1) + hk*nt(j,2) + hl*nt(j,3))*INV
      KK = (hh*nt(j,4) + hk*nt(j,5) + hl*nt(j,6))*INV
      KL = (hh*nt(j,7) + hk*nt(j,8) + hl*nt(j,9))*INV
      phj= (hh*tx(j,1) + hk*tx(j,2) + hl*tx(j,3))*PI2 + ph
      if(IT.eq.0) goto 4
      DO 17 K=1,IT
      IX=ABS(HI(K)-KH)+ABS(KI(K)-KK)+ABS(LI(K)-KL)
   17 IF(IX.EQ.0) GOTO 15
    4 IT=IT+1
      HI(IT)=kh
      KI(IT)=kk
      LI(IT)=kl      
C     fsym = Fobs(j), GF,2 = mask
      fsa= e(i)*cos(phj)
      fsb= e(i)*sin(phj)
      LH = iand((hh-KH),MPAX1)
      LK = iand((hk-KK),MPAY1)
      LL = iand((hl-KL),MPAZ1)
      fas=fas + fsa*GFA(LL,LH,LK,2)-fsb*GFB(LL,LH,LK,2)
      fbs=fbs + fsa*GFB(LL,LH,LK,2)+fsb*GFA(LL,LH,LK,2)
   15 continue
   12 CONTINUE
C     fss = Fext, fs = correction terms, fss-fs = Fext(omit)
      ja = iand(hl,MPAZ1)
      jb = iand(hh,MPAX1)
      jc = iand(hk,MPAY1)
      a=GFA(ja,jb,jc,1) -fas/VOL
      b=GFB(ja,jb,jc,1) -fbs/VOL
      XH = sqrt(a*a+b*b)
      phom = 57.29578*atan2(b,a)
C      write(89,88)ih(i),ik(i),il(i),XH,phom
   88 format(3i4,f9.3,f8.1)   
      a1=a1+nwt*a
      a2=a2+a*a+b*b
      aL=e(i)*cos(p(i))
      be=e(i)*sin(p(i))
      aL1=aL1+nwt*aL
      aL2=aL2+aL*aL+be*be
   14 aa=aa +a*aL +b*be
      a1=a1/NREF
      aL1=aL1/NREF
      a3=dsqrt((a2/NREF-a1*a1)*(aL2/NREF-aL1*aL1))
      a4=dsqrt((a2/NREF)*(aL2/NREF))
      CC = (aa/NREF-a1*aL1)/a3
      CC2= (aa/NREF)/a4
      return
      end
