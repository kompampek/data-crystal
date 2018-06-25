      PROGRAM CCPAT
      CHARACTER AN*10
      parameter (mp=50000,NPAR=128,NPA=NPAR-1)
      integer*2 ih(mp),ik(mp),il(mp),ns,ne
      dimension e(mp),f(mp),p(mp),rho(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ MREF,ih,ik,il,e,f,p,rho,mh,mk,ml,mmk,mml
      COMMON /COR/MXH,MXK,MXL,NH,NK,NL
      WRITE(*,2)
    2 FORMAT(' ENTER DATA FILE NAME : ',$)
      READ(*,3)AN
    3 FORMAT(A)
      WRITE(*,6)
    6 FORMAT(' sig(rho)threshold, FOUR/PAT/PAT-orig (0,1,2) ',$)
      READ(*,*)THRES,NMAP
      t=secnds(0.0)
      OPEN(1,file=AN,form='formatted')
      i=1
    1 read(1,*,end=4)ih(i),ik(i),il(i),ns,ne,a,b,c,np
      MXH=MAX(MXH,ih(i))
      MXK=MAX(MXK,ik(i))
      MXL=MAX(MXL,il(i))
      NH=MIN(NH,ih(i))
      NK=MIN(NK,ik(i))
      NL=MIN(NL,il(i))
      e(i)=0.001*ne
      p(i)=0.001*np
      if(NMAP.NE.0) e(i)=E(I)**2
      if(NMAP.eq.2) e(i)=e(i)-1.0
      i=i+1
      if(i.gt.100000) goto 4
      goto 1
    4 close(1)
      MREF=i-1
      DO 5 I=1,1
      CALL FORWARD(var,NMAP)
      va=var*THRES
      CALL BACK(CC,va,NMAP)
      t3=secnds(t)
    5 write(*,7)THRES,CC,t3
    7 format('THRES = ',f8.2,' CC = ',f8.5,' TIME = ',f5.0)
      END
C      
      SUBROUTINE FORWARD(var,NMAP)
      parameter (NPAR=128,NPA=NPAR-1,NP=NPAR/2,MP=50000)
      byte nh(-NP:NP,-NP:NP),nt(16,9),HH,HK,HL
      real*8 arg,px,ph,vip,vas
      integer*2 ih(mp),ik(mp),il(mp),jh,jk,jl
      dimension r(0:NPA,0:NPA,0:NPA),e(mp)
      dimension data(2*NPAR),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),f(mp)
      dimension tx(16),ty(16),tz(16),p(mp)
      complex fs,cx(0:NPA),fa(0:NP,-NP:NP,-NP:NP)
      COMMON /HKL/ NREF,ih,ik,il,e,f,p,r,mh,mk,ml,mmk,mml
      data mh,mk,ml,mmk,mml,mum/6*0/
      do 23 jh=0,NP
      do 23 jk=-NP,NP
      do 23 jl=-NP,NP
   23 fa(jh,jk,jl)=0.   
      px=-1.0
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MPAX,MPAY,MPAZ
      MPAX1=MPAX-1
      MPAY1=MPAY-1
      MPAZ1=MPAZ-1
      MPAY2=2*MPAY
      MPAZ2=2*MPAZ
      read(3,*)nsym,VOL
      do 16 i=1,nsym
      read(3,17)tx(i),HH,HK,HL,ty(i),KH,KK,KL,tz(i),LH,LK,LL
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
      write(*,*)MPAX,MPAY,MPAZ,nsym,VOL
      px=-1.0
      px=2.0*DACOS(px)/MPAX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPAX-1
      arg=-i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     read reflection data
      do 21 i=1,NREF
      ihh=ih(i)
      ikk=ik(i)
      ill=il(i)
      ff=e(i)
      do 19 j=1,nsym
      jh = ihh*nt(j,1) + ikk*nt(j,2) + ill*nt(j,3)
      jk = ihh*nt(j,4) + ikk*nt(j,5) + ill*nt(j,6)
      jl = ihh*nt(j,7) + ikk*nt(j,8) + ill*nt(j,9)
      tt = ihh*tx(j) +ikk*ty(j) +ill*tz(j)
      tt = mod(tt,1.0)
      ph = p(i) +tt*6.2831853
      if(jh) 13,14,15
   14 if(jk) 13,18,15
   18 if(jl) 13,15,15
   13 jh=-jh
      jk=-jk
      jl=-jl
      ph=-ph
   15 nh(jk,jl)=1
      mh=max(jh,mh)
      mk=max(jk,mk)
      ml=max(jl,ml)
      mmk=min(jk,mmk)
      mml=min(jl,mml)
      if(fa(jh,jk,jl).eq.0.0) mum=mum+1
      if(NMAP.eq.0) fa(jh,jk,jl)=ff*dcmplx(dcos(ph),dsin(ph))
   19 if(NMAP.ne.0) fa(jh,jk,jl)=ff
   21 continue
      write(*,*)mh,mk,ml,mmk,mml,mum
      vip=0
      nip=0
      do 22 jh = 0,mh
      do 22 jk = mmk,mk
      do 22 jl = mml,ml
      fs=fa(jh,jk,jl)
      a=real(fs)
      b=imag(fs)
      ff=sqrt(a*a+b*b)
      vip=vip+ff*ff
   22 if(ff.ne.0) nip=nip+1
      vip=SQRT(vip)/sqrt(2.0)
      if((mh+mh).gt.MPAX.or.(mk+mk).gt.MPAY.or.(ml+ml).gt.MPAZ) RETURN
      write(*,*)nip,vip
C     calculate PATTERSON
      do 1 ix=0,MPAX1
      do 3 ikk = mmk,mk
      do 4 iz=1,MPAZ2
    4 data(iz)=0.
      do 5 ill = mml,ml
      if(nh(ikk,ill).eq.0) goto 5
      fs=0.
      do 2 ihh=0,mh
      imp=ix*ihh
      imp =iand(imp,MPAX1)
    2 fs=fs+fa(ihh,ikk,ill)*cx(imp)
      j1=iand(ill,MPAZ1)
      j1=j1+j1+1
      data(j1)  =real(fs)
      data(j1+1)=imag(fs)
    5 continue
      call FFT(data,MPAZ,1)
      jk=iand(ikk,MPAY1)
      do 6 iz=0,MPAZ1
      izz=iz+iz+1
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
      call FFT(data,MPAY,1)
      do 10 iy=0,MPAY1
      iyy=iy+iy+1
   10 r(iy,iz,ix)=data(iyy)
    7 continue
    1 continue
C     end of FOURIER loops
      vas=0.
      sas=0.
      RMAX=-999.
      RMIN=999.
      do 12 ix=0,MPAX1
      do 12 iz=0,MPAZ1
      do 12 iy=0,MPAY1
      a=r(iy,iz,ix)
      sas=sas+a
      RMIN=MIN(RMIN,A)
      RMAX=MAX(RMAX,A)
      if(a.ne.rmax) goto 12
      jh=ix
      jk=iy
      jl=iz
   12 vas=vas+a*a
      IXYZ=MPAX*MPAZ*MPAY
      sas=sas/IXYZ
      VAR=sqrt(vas/IXYZ)
      T1=RMAX/VAR
      T2=RMIN/VAR
      write(*,*)jh,jk,jl,rmax,sas,IXYZ
      write(*,11)VAR,T1,T2
   11 FORMAT(' VAR,MAX/VAR,MIN/VAR = ',5f10.3)
      RETURN
      END
C     
      SUBROUTINE FFT(DATA,NN,IT)
      PARAMETER (NPAR=128)
C     NPAR MUST BE A POWER OF 2
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI2
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
         TEMPR=SNGL(WR)*DATA(J) -SNGL(WI)*DATA(J+1)
         TEMPI=SNGL(WR)*DATA(J+1) +SNGL(WI)*DATA(J)
         DATA(J)=DATA(I)-TEMPR
         DATA(J+1)=DATA(I+1)-TEMPI
         DATA(I)=DATA(I)+TEMPR
         DATA(I+1)=DATA(I+1)+TEMPI
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
      SUBROUTINE BACK(CC,var,NMAP)
      parameter (NPAR=128,NPA=NPAR-1,mp=50000)
      real*8 arg,px
      integer*2 ih(mp),ik(mp),il(mp),nt(16,9),HH,HK,HL
      integer*2 hi(16),ki(16),li(16)
      dimension r(0:NPA,0:NPA,0:NPA),e(mp),p(mp)
      dimension data(2*NPAR),pa(0:NPA,0:NPA),pb(0:NPA,0:NPA),f(mp)
      complex fs,cx(0:NPA),GF(0:NPA,0:NPA,0:NPA)
      COMMON /HKL/ NREF,ih,ik,il,e,f,p,r,mh,mk,ml,mmk,mml
      COMMON /COR/MXH,MXK,MXL,NH,NK,NL
      px=-1.0
      open(3,file='cell')
      read(3,*)ax,bx,cz
      close(3)
      OPEN(3,file='FFT.dat',status='old',form='formatted')
      read(3,*)MPAX,MPAY,MPAZ
      read(3,*)nsym,VOL
      do 16 i=1,nsym
      read(3,11)tx,HH,HK,HL,ty,KH,KK,KL,tz,LH,LK,LL
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
      px=2.0*DACOS(px)/MPAY
C     initialize exp(-twopi*x/MPAR), remember positive sign ??
      do 13 i=0,MPAY-1
      arg = i*px
   13 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     calculate structure factors
      MPAX1=MPAX-1
      MPAY1=MPAY-1
      MPAZ1=MPAZ-1
      MPAX2=2*MPAX
      MPAY2=2*MPAY
      MPAZ2=2*MPAZ
      fact=2./(MPAX*MPAY*MPAZ)
C
      INK=0
      do 1 jkk= 0,MPAY1
      do 3 jz = 0,MPAZ1
      do 5 jx = 0,MPAX1
      fs=0.
      do 2 jy=0,MPAY1
      rho=r(jy,jz,jx)
      if(rho.lt.var) rho=0.
      if(jkk.eq.0.and.rho.eq.0.0) ink=ink+1
      imp=jkk*jy
      imp=iand(imp,MPAY1)
    2 fs=fs +rho*cx(imp)
      fs=fs*fact
      j1=iand(jx,MPAX1)
      j1=j1+j1+1
      data(j1)  =real(fs)
    5 data(j1+1)=imag(fs)
      call FFT(data,MPAX,-1)
      jk=iand(jz,MPAZ1)
      do 6 jhh=0,MPAX1
      izz=jhh+jhh+1
      pa(jk,jhh)=data(izz)
    6 pb(jk,jhh)=data(izz+1)
    3 continue
      do 7 jhh=0,MPAX1
      do 9 jz =0,MPAZ1
      jk=iand(jz,MPAZ1)
      j1=jk+jk+1
      data(j1)  =pa(jk,jhh)
    9 data(j1+1)=pb(jk,jhh)
      call FFT(data,MPAZ,-1)
      do 10 jll=0,MPAZ1
      iyy=jll+jll+1
   10 GF(jhh,jkk,jll)=cmplx(data(iyy),data(iyy+1))
    7 continue
    1 continue
      IXYZ=MPAX*MPAY*MPAZ
      F000 = FLOAT(IXYZ-ink)/FLOAT(ixyz)
      write(*,*)F000,77
C     end of FOURIER loops
      eo=0
      ec=0
      eo2=0
      ec2=0
      eoec=0
      OPEN(13,file='DATA.OUT')
C      
      do 14 i=1,NREF
      imp=iand(ih(i),mpax1)
      jmp=iand(ik(i),mpay1)
      kmp=iand(il(i),mpaz1)
      fs = GF(imp,jmp,kmp)
      a=real(fs)
      b=imag(fs)
      if(NMAP.eq.0) a = sqrt(a*a+b*b)
      b=e(i)
      eo=eo+b
      ec=ec+a
      eo2=eo2 +b*b
      ec2=ec2 +a*a
      eoec=eoec +a*b
      Xh= fs -F000*e(i)
      GF(imp,jmp,kmp) = -999.
      write(13,4)i,ih(i),ik(i),il(i),b,a,Xh
    4 format(4i5,3f8.3)
   14 continue
      eo=eo/NREF
      ec=ec/NREF
      cc =(eoec/NREF-eo*ec)/sqrt((eo2/NREF-eo*eo)*(ec2/NREF-ec*ec))
      maxh=mh*1.3
      maxk=mk*1.3
      maxl=ml*1.3
      do 8 jh=0,maxh
      do 8 jk=0,maxk
      do 8 jl=0,maxl
      stl = 0.5*SQRT((jh/ax)**2 +(jk/bx)**2 +(jl/cz)**2)
      if(stl.gt.0.6) goto 12
      fs=gf(jh,jk,jl)
      aa=real(fs)
      bb=imag(fs)
      if(NMAP.eq.0) aa=sqrt(aa*aa+bb*bb)
      del=abs(999.+aa)
      if(del.gt.0.01)write(37,15)jh,jk,jl,stl,aa
   15 format(3i5,f8.4,f10.3)   
   12 CONTINUE
    8 CONTINUE   
      RETURN
      END
