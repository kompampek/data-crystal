C     PROGRAM RHOMin FLIP
      PARAMETER (NPAR=128,NPA=NPAR-1,NP=63,NRF=50000)
      BYTE nh(0:NP,-NP:NP),ih(NRF),ik(NRF),il(NRF)
      BYTE mh,mk,ml,mmk,mml,jh,jk,jl
      dimension E(NRF),p(NRF)
      complex f(0:NP,-NP:NP,-NP:NP)
      COMMON/hkl/MX,MY,MZ,ih,ik,il,e,p,f,nh,mh,mk,ml,mmk,mml
      data mh,mk,ml,mmk,mml,var/3*-99,2*99,0.0/
      open(1,file='fft.dat')
      read(1,*)MX,MY,MZ,IHORZ,IVERT,ISECT
      close(1)
      do 25 i=0,NP
      do 25 j=-NP,NP
      nh(i,j)=0
      do 25 k=-NP,NP
      f(i,j,k)=0.
   25 continue
      open(unit=1,file='P1.E',status='old',form='formatted')
      do 26 i=1,50000
      read(1,*,end=27)jh,jk,jl,ns,ne,ff,sf,fc,nph
      call SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
      p(i)=0.001*nph
      ff = 0.001*ne
      E(i)=ff
      if(jh.gt.0) goto 29
      if(jh.eq.0) goto 28
      jh=-jh
      jk=-jk
      jl=-jl
      p(i)=-p(i)
      goto 29
   28 if(jk.gt.0) goto 29
      if(jk.eq.0) goto 30
      jk=-jk
      jl=-jl
      p(i)=-p(i)
      goto 29
   30 if(jl.gt.0) goto 29
      jl=-jl
      p(i)=-p(i)
   29 nh(jh,jk)=1
      mh=max(jh,mh)
      mk=max(jk,mk)
      ml=max(jl,ml)
      mmk=min(jk,mmk)
      mml=min(jl,mml)
      ih(i)=jh
      ik(i)=jk
      il(i)=jl
      var=var+ff*ff
   26 f(jh,jk,jl) = ff*cmplx(cos(p(i)),sin(p(i)))
   27 nref=i-1
      var=sqrt(var/2.)
      write(*,*)var,mh,mk,ml,mmk,mml,nref
      close(1)
      read(*,*)AMOD,ncycle
      do 1 i=1,ncycle
      call forward(AMOD,amax,amin)
C     flip all grid point intensities less than AMOD*sig(RHO)
      call back(sum)
    1 continue
      end
C
      SUBROUTINE FORWARD(AMOD,AMAX,AMIN)
      parameter (NPAR=128,NPA=NPAR-1,NP=63,NRF=50000)
      real*8 arg,px,pi2
      BYTE ih,ik,il,lh(NRF),lk(NRF),ll(NRF),mh,mk,ml,mmk,mml
      byte nh(0:NP,-NP:NP)
      integer*2 nm(0:NPA)
      dimension r(0:NPA,0:NPA,0:NPA),E(NRF),p(NRF)
      complex f(0:NP,-NP:NP,-NP:NP),ppp(0:NPA,0:NPA,0:NPA)
      complex cx(0:NPA),fx(0:NPA),pp(0:NPA,0:NPA)
      COMMON /hkl/MX,MY,MZ,lh,lk,ll,e,p,f,nh,mh,mk,ml,mmk,mml
      COMMON /RHO/ r
      data pi2/6.283185308/
      MPA=MX-1
      px=pi2/MX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPA
      arg=-i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     set up bit inversion arrays
      nm(0)=0
      NX = alog(float(MPA))/alog(2.0) +.1
      mm=1
      do 23 i=1,NX
      do 24 j=0,mm-1
      nm(j)=2*nm(j)
   24 nm(j+mm) = nm(j)+1
   23 mm=mm+mm
C     read reflection data
C     calculate FOURIER
      do 1 ih=0,mh
      do 2 iz=0,MPA
      do 2 iy=0,MPA
    2 pp(iy,iz)=0.0
      do 3 ik = mmk,mk
      if(nh(ih,ik).eq.0) goto 3
      do 4 iz=0,MPA
    4 fx(iz)=0.
      do 5 il = mml,ml
      ill=iand(il,MPA)
    5 fx(nm(ill))=f(ih,ik,il)
      call FFT(MX,fx,cx,NX)
      write(97,*)(fx(i),i=0,MPA)
      jk=iand(ik,MPA)
      do 6 iz=0,MPA
    6 pp(jk,iz) = fx(iz)
    3 continue
      do 7 iz=0,MPA
      do 8 ik=0,MPA
    8 fx(ik)=0.
      do 9 ik = mmk,mk
      jk=iand(ik,MPA)
    9 fx(nm(jk))=pp(jk,iz)
      write(98,*)(fx(i),i=0,MPA)
      call FFT(MX,fx,cx,NX)
      do 10 iy=0,MPA
   10 ppp(ih,iy,iz)=fx(iy)
    7 continue
    1 continue
      var2=0.
      rhmax=-99.
      rhmin=99.
      do 11 iz=0,MPA
      do 11 iy=0,MPA
      do 12 ih=0,MPA
   12 fx(ih)=0.
      do 13 ih=0,mh
   13 fx(nm(ih))=ppp(ih,iy,iz)
      call FFT(MX,fx,cx,NX)
      do 14 ix=0,MPA
      rho=fx(ix)
      rhmin=min(rho,rhmin)
      rhmax=max(rho,rhmax)
      var2=var2+rho*rho
   14 r(ix,iy,iz)=rho
   11 continue
C     end of FOURIER loops
      var2=SQRT(var2/(MPA*MPA*MPA))
      write(*,*)MX,MY,MZ,var2,99,rhmax,rhmin
      return
      end
C     
C     RADIX-2 FFT - ALL F(H)*EXP(-ip) ARE TREATED AS -F(H)
      SUBROUTINE FFT(N,F,CX,NX)
      PARAMETER (NPAR=128, NPA=NPAR-1)
      COMPLEX F(0:NPA),CX(0:NPA),A
      I1=1
      M2=N/2
      DO 1 I=1,NX
      I0=I1
      I2=I1-1
      I1=I1+I1
      DO 2 J=0,N-1,I1
      A=F(J+I0)
      F(J+I0)=F(J)-A
      F(J)=F(J)+A
      I4=0
      DO 3 K=J+1,J+I2
      I4=I4+M2
      A = F(K+I0)*CX(I4)
      F(K+I0)=F(K)-A
    3 F(K)=F(K)+A
    2 CONTINUE
    1 M2=M2/2
      RETURN
      END
C
      SUBROUTINE SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
      BYTE jh,jk,jl,nv(3)
      nv(IHORZ) = jh
      nv(IVERT) = jk
      nv(ISECT) = jl
      jh = nv(1)
      jk = nv(2)
      jl = nv(3)
      return
      end
C
      SUBROUTINE BACK(SUM)
      parameter (NPAR=128,NPA=NPAR-1,NP=63,NRF=50000)
      real*8 arg,px,pi2
      BYTE nh(0:NP,-NP:NP),lh(NRF),lk(NRF),ll(NRF),ih,ik,il
      BYTE mh,mk,ml,mmk,mml
      integer*2 nm(0:NPA)
      dimension r(0:NPA,0:NPA,0:NPA),E(NRF),p(NRF)
      complex cx(0:NPA),fx(0:NPA),pp(0:NPA,0:NPA)
      complex fcal,ppp(0:NPA,0:NPA,0:NPA),f(0:NP,-NP:NP,-NP:NP)
      COMMON /hkl/MX,MY,MZ,lh,lk,ll,e,p,f,nh,mh,mk,ml,mmk,mml
      COMMON /RHO/ r
      data pi2/6.283185308/
      MPA=MX-1
      px=pi2/MX
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPA
      arg= i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     set up bit inversion arrays
      NX = alog(float(MPA))/alog(2.0) +.1
      mm=1
      nm(0)=0
      do 23 i=1,NX
      do 24 j=0,mm-1
      nm(j)=2*nm(j)
   24 nm(j+mm) = nm(j)+1
   23 mm=mm+mm
      fact = 1./(MX**3)
      write(*,*)NX,fact,mpa
C     calculate FOURIER
      do 1 ih=0,MPA
      do 2 iz=0,MPA
      do 2 iy=0,MPA
    2 pp(iy,iz)=0.0
      do 3 ik = 0,MPA
      do 4 iz = 0,MPA
    4 fx(iz)=0.
      do 5 il = 0,MPA
      ill = iand(il,MPA)
    5 fx(nm(ill))=r(ih,ik,il)
      call FFT(MX,fx,cx,NX)
      jk=iand(ik,MPA)
      do 6 iz=0,MPA
    6 pp(jk,iz) = fx(iz)
    3 continue
      do 7 iz=0,MPA
      do 8 ik=0,MPA
    8 fx(ik)=0.
      do 9 ik = 0,MPA
      jk=iand(ik,MPA)
    9 fx(nm(jk))=pp(jk,iz)
      call FFT(MX,fx,cx,NX)
      do 10 iy=0,MPA
   10 ppp(ih,iy,iz)=fx(iy)
    7 continue
    1 continue
      do 11 iz=0,MPA
      do 11 iy=0,MPA
      do 12 ih=0,MPA
   12 fx(ih)=0.
      do 13 ih=0,MPA
   13 fx(nm(ih))=ppp(ih,iy,iz)
      call FFT(MX,fx,cx,NX)
      do 14 ix=0,MPA
   14 ppp(ix,iy,iz)=fx(ix)
   11 continue
C     end of inverse FOURIER loops
      do 15 i=1,50,4
      ih=iand(lh(i),MPA)
      ik=iand(lk(i),MPA)
      il=iand(ll(i),MPA)
      fcal=ppp(ih,ik,il)
      a=fcal
      b=imag(fcal)
      fc=2.*sqrt(a*a+b*b)*fact
      ph=atan2(b,a)
      fcal=f(lh(i),lk(i),ll(i))
      a=fcal
      b=imag(fcal)
      fo=sqrt(a*a+b*b)
      ph0=atan2(b,a)
   15 write(*,16)lh(i),lk(i),ll(i),e(i),p(i),fc,ph,fo,ph0
   16 format(3i4,6f8.3)   
      return
      end
