C     PROGRAM RHOMOD
      parameter (NPAR=128,NPA=NPAR-1,NP=40)
      dimension r(0:NPA,0:NPA,0:NPA)
      COMMON /RHO/ r,amax,amin,sig,amod
      read(*,*)AMOD,ncycle
      do 1 i=1,ncycle
      call forward(av)
      call back(sum)
      write(*,2)amax,amin,sig,sum,av
    2 format(' RHOMX/MIN/SIG = ',2f6.2,f9.6,' e/A**3, DELPHI = ',2f7.2)
    1 continue
      end
C
      SUBROUTINE FORWARD(av)
      parameter (NPAR=128,NPA=NPAR-1,NP=40)
      byte nh(0:NP,-NP:NP),HH,HK,HL,nt(16,9)
      real*8 arg,px,pi2,sar,var,a
      integer*2 jh,jk,jl,ih,ik,il,ns,ne,nph,nm(0:NPA)
      dimension r(0:NPA,0:NPA,0:NPA),tx(16,3)
      complex f(0:NP,-NP:NP,-NP:NP),ppp(0:NPA,0:NPA,0:NPA)
      complex cx(0:NPA),fx(0:NPA),pp(0:NPA,0:NPA)
      COMMON /RHO/ r,amax,amin,sig,amod
      data pi2/6.283185308/
      open(unit=1,file='fft.dat',status='old',form='formatted')
      READ(1,*)MPAR,MPAX,IHORZ,IVERT,ISECT
      close(unit=1)
      MPA=MPAR-1
      px=pi2/MPAR
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPA
      arg=-i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     set up bit inversion arrays
      nm(0)=0
      NX = alog(float(MPAR))/alog(2.0) +.1
      mm=1
      do 23 i=1,NX
      do 24 j=0,mm-1
      nm(j)=2*nm(j)
   24 nm(j+mm) = nm(j)+1
   23 mm=mm+mm
C     read reflection data
      sar=0.
      num=0
      open(unit=1,file='P1.E',status='old',form='formatted')
      do 26 i=1,999999
      read(1,*,end=27)ih,ik,il,ns,ne,ff,sf,fc,nph
      ph=0.001*nph
      ff=0.002*ne
      call SWAP(ih,ik,il,IHORZ,IVERT,ISECT)
      f(ih,ik,il) = ff*cmplx(cos(ph),sin(ph))
      sar=sar+f(ih,ik,il)**2
   26 continue
   27 close(1)
      num=i-1
      var=sqrt(2.*sar)
      if(num.ne.0) write(*,*)num,var
      if(num.ne.0) sig=var
C
C     calculate FOURIER
C
      do 1 ih=0,mh
      write(*,*)1
      do 2 iz=0,MPA
      do 2 iy=0,MPA
    2 pp(iy,iz)=0.0
      write(*,*)11
      do 3 ik = mmk,mk
      if(nh(ih,ik).eq.0) goto 3
      do 4 iz=0,MPA
    4 fx(iz)=0.
      do 5 il = mml,ml
      jl=iand(il,MPA)
    5 fx(nm(jl))=f(ih,ik,il)
      write(*,*)2
      CALL FFT(MPAR,fx,cx,NX)
      jk=iand(ik,MPA)
      do 6 iz=0,MPA
    6 pp(jk,iz) = fx(iz)
    3 continue
      write(*,*)3
      do 7 iz=0,MPA
      do 8 ik=0,MPA
    8 fx(ik)=0.
      do 9 ik = mmk,mk
      jk=iand(ik,MPA)
    9 fx(nm(jk))=pp(jk,iz)
      CALL FFT(MPAR,fx,cx,NX)
      do 10 iy=0,MPA
   10 ppp(ih,iy,iz)=fx(iy)
    7 continue
    1 continue
    
      write(*,*)88
      do 11 iz=0,MPA
      do 11 iy=0,MPA
      do 12 ih=0,MPA
   12 fx(ih)=0.
      do 13 ih=0,mh
   13 fx(nm(ih))=ppp(ih,iy,iz)
      call FFT(MPAR,fx,cx,NX)
      do 14 ix=0,MPA
   14 r(ix,iy,iz)=fx(ix)
   11 continue
C     end of FOURIER loops, DENSITY MODIFICATION ?
      write(*,*)99
      amax=-999.
      amin=9999.
      var=0.
      av=0.
      do 18 iz=0,MPA
      do 18 iy=0,MPA
      do 18 ix=0,MPA
      a=r(ix,iy,iz)
      if(a.lt.amin) amin=a
      if(a.gt.amax) amax=a
c      aa=(a+.4)**2
C      var=var+aa
C
C     DENSITY MODIFICATION OCCURS HERE !
C
c      prob = (1.-exp(aa*var2))
C      av=av+prob
C      r(ix,iy,iz) = a*prob
      if(a.lt.amod) r(ix,iy,iz) = 0.0
   18 continue
      write(*,*)77
      av=av/(MPAR**3)
c      sig=sqrt(var/(MPAR**3))
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
      integer*2 jh,jk,jl,nv(3)
      nv(IHORZ) = jh
      nv(IVERT) = jk
      nv(ISECT) = jl
      jh = nv(1)
      jk = nv(2)
      jl = nv(3)
      RETURN
      END
C
      SUBROUTINE BACK(SUM)
      parameter (NPAR=128,NPA=NPAR-1,NH=NPA)
      real*8 arg,px,pi2
      integer*2 nm(0:NPA),ih,ik,il,ns,ne,np,mp
      dimension r(0:NPA,0:NPA,0:NPA)
      complex cx(0:NPA),fx(0:NPA),pp(0:NPA,0:NPA)
      complex fcal,ppp(0:NH,0:NH,0:NH)
      COMMON /RHO/ r,amax,amin,sig,amod
      data pi2/6.283185308/
      open(unit=1,file='fft.dat',status='old',form='formatted')
      READ(1,*)MPAR,MPAX,IHORZ,IVERT,ISECT
      READ(1,*) NSYM,VOL
      close(unit=1)
      MPA=MPAR-1
      px=pi2/MPAR
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPA
      arg= i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     set up bit inversion arrays
      NX = alog(float(MPAR))/alog(2.0) +.1
      mm=1
      do 23 i=1,NX
      do 24 j=0,mm-1
      nm(j)=2*nm(j)
   24 nm(j+mm) = nm(j)+1
   23 mm=mm+mm
      fact = VOL/(MPAR**3)
C     calculate FOURIER
      do 1 ih=0,MPA
      do 2 iz=0,MPA
      do 2 iy=0,MPA
    2 pp(iy,iz)=0.0
      do 3 ik = 0,MPA
      do 4 iz=0,MPA
    4 fx(iz)=0.
      do 5 il = 0,MPA
      jl=iand(il,MPA)
    5 fx(nm(jl))=fact*r(ih,ik,il)
      call FFT(MPAR,fx,cx,NX)
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
      call FFT(MPAR,fx,cx,NX)
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
      call FFT(MPAR,fx,cx,NX)
      do 14 ix=0,MPA
   14 ppp(ix,iy,iz)=fx(ix)
   11 continue
C     end of inverse FOURIER loops
C
C     lines added to check if phases converged towared LS-ref values
C
      open(1,file='fourr.tmp',form='formatted')
      rad=180./3141
      sum=0.
      do 18 i=1,999999
      read(1,*,end=19)ih,ik,il,ns,ne,a,b,c,np
      jk=iand(ik,MPA)
      jl=iand(il,MPA)
      fcal = ppp(ih,jk,jl)
      aa =  real(fcal)
      bb = aimag(fcal)
      ix=ih*ik*il
      if(ix.ne.0) goto 17
C
C     reset zonal phases for space group P212121
C
      if(abs(aa).lt.abs(bb)) aa=0.
      if(abs(bb).lt.abs(aa)) bb=0.
   17 mp = 1000.*atan2(bb,aa)
      nd=mod(abs(np-mp),6283)
      if(nd.gt.3141) nd=6283-nd
      sum=sum+rad*nd
   18 CONTINUE
   19 num=i-1
      close(1)
      close(2)
      sum=sum/num
      return
      end
