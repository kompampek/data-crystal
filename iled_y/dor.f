      PROGRAM SnB3
C     COMPARE 3 DIFFERENT Rmin TARGET FUNCTIONS, SAVE PHASE SETS
      parameter (mp=9000,mt=200000,mp2=18001,mt3=3*mt)
      character fil*6
      integer*4 IPR(9,3,3),ih(mp),ik(mp),il(mp),nr(mp)
      integer*4 iseed,MM(0:mp),n0,nm(mt3),M(mp),NN(3,mt)
      real*8 x,y,z,rand
      dimension p(-mp:mp),tt(mt),an(mt),e(mp),TR(9,3)
      dimension tnt(mt),r(4),true(mp),Q(mp)
      dimension cs(0:65536),sn(0:65536),xx(1500,3)
      COMMON /SIGNS/ cs,sn,xx,NLIST
      COMMON /HKL/ e,p,tr,NREF,NTOTAL,ih,ik,il,NSYM,ipr
      data M,MM,tmap,pksum,tSF/mp2*0,3*0.0/
C     extensions E=.eee, T=.trp, A=.amb, Q=.nqt, P=.pks, X=.xyz     
      open(unit=2,file='fft0.dat',status='old',form='formatted')
      read(2,*)NSYM
      read(2,*)NSYM,VOL
      do 4 i=1,NSYM
    4 read(2,19) (tr(i,j),(ipr(i,j,k),k=1,3),j=1,3)
   19 format(3(f10.5,3i5))
      close(unit=2)
      pi = acos(-1.0)
      pi2=2.*pi
      ARG=pi/32768.
      do 32 i=0,65536
      ar=i*arg
      sn(i)=sin(ar)
   32 cs(i)=cos(ar)
      pi6=pi/6.
      pion2=pi/2.
      write(*,27)
   27 format(' Enter 4-letter registry name : ',$)
      read(*,25)fil
   25 format(A)
      fil(5:5)='.'
      fil(6:6)='A'
      open(unit=3,file=fil,status='unknown',form='formatted')
      write(*,11)
   11 format(' IRAN,NAT,PS,CYC,SETS,LIST,MAT,NST,NOP(0,1),ESC:',$)
C     ix odd > 0 to scan random phase sets
      read(*,*)id,NATOM,iter,mcyc,msets,list,MATOM,NST,NOP,ESC
      write(3,24)id,NATOM,iter,mcyc,msets,list,MATOM,NST,NOP,ESC
   24 format(9i5,f5.1)
      if((MATOM+NST).eq.0) write(*,213)
  213 format(' UNFORTUNATELY BOTH MAT & NST EQUAL ZERO')
      write(*,214)
  214 format(' DO PHI INCREMENTS START AT PHIZER0 (Y/N = 0,1)')
      read(*,*)NICK      
      if(id.lt.0) open(33,file='IRAN-NOS')
      ix=abs(id)
      iran=ix
      ncyc=0
C     READ SYM FILE
      fil(6:6)='E'
      open(unit=1,file=fil,status='old',form='formatted')
      do 1 i=1,mp
      read(1,*,end=29)ih(i),ik(i),il(i),ns,ne,a,b,c,np
      nr(i)=13
      true(i)=0.001*np
      do 22 j=2,NSYM
      i1=(ih(i)*ipr(j,1,1)+ik(i)*ipr(j,2,1)+il(i)*ipr(j,3,1))
      i2=(ih(i)*ipr(j,1,2)+ik(i)*ipr(j,2,2)+il(i)*ipr(j,3,2))
      i3=(ih(i)*ipr(j,1,3)+ik(i)*ipr(j,2,3)+il(i)*ipr(j,3,3))
      if((abs(ih(i)+i1)+abs(ik(i)+i2)+abs(il(i)+i3)).ne.0) goto 22
      k = NINT(12*mod(abs(i1*tr(j,1)+i2*tr(j,2)+i3*tr(j,3)),1.0))
      if(k.gt.6) k=12-k
      nr(i)=k
      goto 1
   22 continue
    1 e(i)=0.001*ne
   29 close(unit=1)
      NTOTAL=i-1
C     READ TRIPLES FILE
      fil(6:6)='T'
      open(unit=1,file=fil,status='old',form='formatted')
      NREF=0
      bot=0.
      trip=0.
      ntrip=0
      do 33 i=1,mt
      read(1,*,end=2)a1,m1,m2,n1,n2,n3
      if(m1.eq.7) goto 33
      ntrip=ntrip+1
      nn(1,ntrip)=n1
      nn(2,ntrip)=n2
      nn(3,ntrip)=n3
      t=0.
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      if(m1.eq.1) T = tanh(0.5*a1)
      tnt(ntrip)= M2*pi6
      tt(ntrip)=t
      n2=iabs(n2)
      n3=iabs(n3)
      an(ntrip)=A1
      bot=bot + A1
      trip=trip + A1*t
      if(n3.gt.NREF) NREF=n3 
      m(n1)=m(n1)+1
      m(n2)=m(n2)+1
      m(n3)=m(n3)+1
   33 continue
    2 trip=trip/bot
      close(unit=1)
      esq=0.
      do 34 i=1,NREF
      esq=esq+e(i)*e(i)
   34 e(i)=2.*e(i)/VOL
      var=NSYM*esq/(VOL*VOL)
      n0=0
      do 12 i=1,NREF
      n0=n0+m(i)
   12 mm(i)=n0
      do 13 i=1,ntrip
      do 14 j=1,3
      n1=iabs(nn(j,i))
      m(n1)=m(n1)-1
      n0=mm(n1)-m(n1)
   14 nm(n0)=i
   13 continue
      write(*,36)NREF,trip,ntrip
   36 format(i5,' phs, Cosav =',f6.3,' for ',i6,' trips',//
     *'   Amb# Ncyc  Rmin  Cosav  <|dp|> '/)
CC
      call cpu_time(T0)
C
C     LOOP OVER THE NUMBER OF DIFFERENT TRIALS
C
      do 15 nsets=1,msets
      axe=99.
      if(id.lt.0) read(33,*)iran
      NLIST=0
      if(nsets.eq.msets) NLIST=1
      MAD=NST
C
C     GENERATE RANDOM ATOMS or read XYZ file
C
      IF(MATOM.eq.0) GOTO 31
      open(unit=4,file='FRAG.xyz',status='old',form='formatted')
      NEND = MATOM*(IRAN-1)/2
      DO 28 I=1,NEND
   28 read(4,*)X
      DO 26 I=1,MATOM
   26 read(4,*)xx(i,1),xx(i,2),xx(i,3)
      CLOSE(4)
      MAD=MATOM
      GOTO 30
   31 ix=iran
      do 7 i=1,NST
      iy=899*ix
      iy=iand(iy,65535)
      if(iy.lt.0) iy=iy+32768
      iseed=iy
      call srand(iseed)
      x=rand()
      y=rand()
      z=rand()
      xx(i,1)=x
      xx(i,2)=y
      xx(i,3)=z
    7 ix=iy
   30 call SFCLC(MAD,MAD)
      close(unit=1)
C
C     THIS IS THE MAIN CYCLES LOOP
C
      do 10 ncyc=1,mcyc
      irest=mcyc-ncyc
      do 21 nit=0,iter
      if(nit.eq.0) goto 20
C     PARMETER SHIFT - INNER LOOP 16
      do 16 i=1,NREF
      do 17 j=1,4
   17 r(j)=0.
      iseed = i+nit
      call srand(iseed)
      x=rand()
      x=NICK*x
      phi0=6.28318*x
      INKPHI = 65535*x
      do 18 ij=mm(i-1)+1,mm(i)
      k=nm(ij)
      ank=an(k)
      ttk=tt(k)
      n2=nn(2,k)
      n3=nn(3,k)
      narg = (p(nn(1,k))+p(n2)+p(n3)+tnt(k))/arg
      if(i.eq.-n2.or.i.eq.-n3) narg=-narg
      narg=narg+inkphi
      do 9 j=1,4
      narg=narg+16384
      narg=iand(narg,65535)
      cost=cs(narg)
      sint2=0.
      if(NOP.eq.1.and.nr(i).eq.13) sint2=sn(narg)**2
    9 if(NOP.ne.2) r(j)=r(j)+ank*(sint2 + (cost-ttk)**2 )
   18 CONTINUE
      rmin=r(1)
      ind=1
C     NOP=0 Rmin, NOP=1 A*((cos2-t)**2 +sin**2) , NOP=2 statistical    
      do 8 j=1,3
      if(irest.lt.3.and.nr(i).ne.13.and.j.ne.2) goto 8
      if(r(j).gt.rmin) goto 8
      rmin=r(j)
      ind=j
    8 continue
      p(i)=p(i)+phi0+ind*pion2
      p(-i)=-p(i)
   16 continue
      goto 21
C     CALCULATE Rmin, COSAV, TRIPLES Phase error
   20 continue
      Rold=0.
      Cosav=0.
      do 3 i=1,ntrip
      a1=an(i)
      narg=(p(nn(1,i))+p(nn(2,i))+p(nn(3,i))+tnt(i))/arg
      narg=iand(narg,65535)
      cost=cs(narg)
      Cosav=Cosav+a1*cost
      cost=tt(i)-cost
    3 Rold=Rold + a1*cost*cost
      Rold=Rold/bot
      Cosav = Cosav/bot
      kcyc=ncyc
      a1=999.
      do 35 i=0,15
      i1=iand(i,1)
      i0=i/2
      i2=iand(i0,1)
      i0=i0/2
      i3=iand(i0,1)
      i4=i0/2
      i4=i4+i4-1
      a2=0.
      do 355 j=1,NREF
      it=abs(ih(j)*i1+ik(j)*i2+il(j)*i3)
      it=iand(it,1)
      dp=cos(true(j)-i4*p(j) +it*PI)
  355 a2=a2+acos(dp)
      a2=57.296*a2/nref
      if(a2.gt.a1) goto 35
      a1=a2
      j1=i1
      j2=i2
      j3=i3
      j4=i4
   35 CONTINUE   
      angle = a1
      if(ncyc.gt.10) axe=min(axe,angle)
      if(angle.gt.axe) goto 356
      do 357 i=1,NREF
      it=abs(ih(i)*j1+ik(i)*j2+il(i)*j3)
      it=iand(it,1)
  357 q(i)=j4*p(i)+it*PI
  356 MIST=1
      if(list.ne.0) mist=mod(ncyc,10)*LIST
      if(mist.EQ.0.) write(*,5)iran,ncyc,Rold,Cosav,angle,axe
    5 format(i7,i4,2f7.3,f6.1,f6.1)
      if(angle.lt.ESC) goto 789
   21 CONTINUE 
      if(NLIST.eq.1.and.ncyc.eq.mcyc) NLIST=2
      call cpu_time(t1)
      CALL EDMAP(pksum,VAR)
      call cpu_time(t2)
      tmap=tmap+t2-t1
      MAD=NATOM
      call SFCLC(NATOM,MAD)
      call cpu_time(t3)
      tSF=tSF+t3-t2
   10 CONTINUE
  789 write(3,5)iran,kcyc,rold,cosav,angle,axe
      if(axe.lt.65.0) goto 790
      write(13,5)iran,kcyc,rold,cosav,angle,axe
      write(13,37) (q(i),i=1,NREF)
   37 format(20f6.2)   
  790 iran=iran+2
   15 continue
      call cpu_time(TT0)
      tsum=TT0-T0
      if(tsum.lt.0.0) tsum=tsum+1440.
      tmap=tmap-pksum
      tPS=(tsum-tmap-pksum-tSF)/60.
      tsum=tsum/60.
      tmap=tmap/60.
      pksum=pksum/60.
      tSF=tSF/60.
      write(*,6)tsum,tPS,tSF,tmap,pksum
      write(3,6)tsum,tPS,tSF,tmap,Pksum
    6 format(' Total time =',f9.3,' Minutes'/' PS =',f9.3,' SF =',f9.3, 
     *' FR =',f9.3,' PPK =',f9.3,' Minutes')
      CLOSE(UNIT=3)
      END
C
      SUBROUTINE SFCLC(NAT,NAT2)
      parameter (mp=9000)
      integer*4 jh(mp),jk(mp),jl(mp),IPR(9,3,3),nx(4500,16)
      dimension e(mp),p(-mp:mp),TR(9,3)
      dimension cs(0:65536),sn(0:65536),xxx(1500,3)
      COMMON /SIGNS/ cs,sn,xxx,NLIST
      COMMON /HKL/ e,p,tr,NREF,NTOTAL,jh,jk,jl,NSYM,ipr
      ND = NAT-NAT2
      do 1 ii=1,NAT2
      i=ii+ND
      xx=xxx(i,1)
      yy=xxx(i,2)
      zz=xxx(i,3)
      i3=3*i
      i2=i3-1
      i1=i3-2
      do 2 j=1,NSYM
      x = ipr(j,1,1)*xx +ipr(j,1,2)*yy +ipr(j,1,3)*zz +tr(j,1)
      nx(i1,j)=65536.*amod(x,1.0)
      y = ipr(j,2,1)*xx +ipr(j,2,2)*yy +ipr(j,2,3)*zz +tr(j,2)
      nx(i2,j)=65536.*amod(y,1.0)
      z = ipr(j,3,1)*xx +ipr(j,3,2)*yy +ipr(j,3,3)*zz +tr(j,3)
    2 nx(i3,j)=65536.*amod(z,1.0)
    1 continue
      NAT3=3*NAT2
      do 5 i=1,NREF
      ih=jh(i)
      ik=jk(i)
      il=jl(i)
      aa=0.
      bb=0.
      do 3 k=1,NSYM
      do 4 j=1,NAT3,3
      narg = ih*nx(j,k) + ik*nx(j+1,k) + il*nx(j+2,k)
      narg=iand(narg,65535)
      aa = aa + cs(narg)
    4 bb = bb + sn(narg)
    3 continue
      p(i)=atan2(bb,aa)
    5 p(-i)=-p(i)
      RETURN
      END
C
C     3-D SECTIONAL FFT
      SUBROUTINE EDMAP(pks,VAR)
C     RADIX-2 FFT  -  PLUS PEAK LISTING SUBROUTINES
      parameter (NPAR=256,NPA=NPAR-1,NP=80,mp=9000)
      byte nt(24,9),nh(-NP:NP,-NP:NP),HH,HK,HL
      real*8 arg,pi2,px,dp
      integer*4 ih(mp),ik(mp),il(mp),nmy(0:NPA),nmz(0:NPA),IPR(9,3,3)
      dimension r(-1:NPAR,-1:NPAR,3),t(4,4),tx(24,3),TR(9,3),tq(4,4)
      dimension e(mp),ph(-mp:mp)
      complex f(0:NP,-NP:NP,-NP:NP),cx(0:NPA),cy(0:NPA),cz(0:NPA)
      complex p(0:NPA,0:NPA),fx(0:NPA),fhkl,fs
      COMMON /WORK/ fx,cx,cy,cz
C      COMMON /PAD/ nmy,nmz,MPX1,MPY1,MPZ1,NH,NT
      COMMON /PKS/ r,tq,peak,MPX,MPAX,MPY,MPZ,IHORZ,IVERT,ISECT,NATOM
      COMMON /HKL/ e,ph,tr,NREF,NTOTAL,ih,ik,il,NSYM,ipr
      data mh,mk,ml,mmk,mml/5*0/
      data t/18.,4*-10.,11.,2*4.,-10.,4.,11.,4.,-10.,2*4.,11./
      npeaks=0
      rlast=0.5*var
      peak=0.5*var
      f(0,0,0)=0.1
c      if(t(1,1).ne.18.) goto 21
      open(unit=1,file='fft0.dat',status='old',form='formatted')
      READ(1,*)MPX,MPAX,MPY,MPZ,IHORZ,IVERT,ISECT,NATOM
      READ(1,*)NSYM,VOL
      MPX1=MPX-1
      MPY1=MPY-1
      MPZ1=MPZ-1
      do 16 i=1,nsym
      READ(1,17)tx(i,1),HH,HK,HL,tx(i,2),KH,KK,KL,tx(i,3),LH,LK,LL
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
      nt(i,9) = (HH*KK-HK*KH)/IDET
   16 continue
      close(unit=1)
      px=-1.0
      px=2.0*dacos(px)
      pi2=px
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 14 it=1,3
      jit=MPX1
      if(it.eq.2) jit=MPY1
      if(it.eq.3) jit=MPZ1
      dp=px/(jit+1)
      do 20 i=0,jit
      arg=-i*dp
      if(it.eq.1) cx(i)=dcmplx(dcos(arg),dsin(arg))
      if(it.eq.2) cy(i)=dcmplx(dcos(arg),dsin(arg))
   20 if(it.eq.3) cz(i)=dcmplx(dcos(arg),dsin(arg))
   14 continue
C     set up FFT arrays for bit reversal on first pass
      do 12 it=1,2
      MPAR=MPY
      if(it.eq.2) MPAR=MPZ
      if(it.eq.1) nmy(0)=0
      if(it.eq.2) nmz(0)=0
      NX = alog(float(MPAR))/alog(2.0) - 0.9
      mm=1
      do 23 i=1,NX+1
      do 24 j=0,mm-1
      if(it.eq.1) nmy(j)=2*nmy(j)
      if(it.eq.2) nmz(j)=2*nmz(j)
      if(it.eq.1) nmy(j+mm) = nmy(j)+1
   24 if(it.eq.2) nmz(j+mm) = nmz(j)+1
      mm=mm+mm
   23 continue
   12 continue
      do 19 i=1,4
      do 19 j=1,4
   19 tq(i,j)=t(i,j)/42.
   21 continue
      do 26 i=1,NREF
      do 28 j=1,nsym
      jh = ih(i)*nt(j,1) + ik(i)*nt(j,2) + il(i)*nt(j,3)
      jk = ih(i)*nt(j,4) + ik(i)*nt(j,5) + il(i)*nt(j,6)
      jl = ih(i)*nt(j,7) + ik(i)*nt(j,8) + il(i)*nt(j,9)
      sh = ih(i)*tx(j,1) + ik(i)*tx(j,2) + il(i)*tx(j,3)
      sh = mod(sh,1.0)
      phi=ph(i) + sh*pi2
      pie2=pi2
      phi=amod(phi,pie2)
      call SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
      if(jh) 32,30,29
   30 if(jk) 32,31,29
   31 if(jl) 32,29,29
   32 jh=-jh
      jk=-jk
      jl=-jl
      phi=-phi
   29 nh(jk,jl)=1
      if(t(1,1).ne.18.) goto 28
      mh=max(jh,mh)
      mk=max(jk,mk)
      ml=max(jl,ml)
      mmk=min(jk,mmk)
      mml=min(jl,mml)
   28 f(jh,jk,jl) = e(i)*cmplx(cos(phi),sin(phi))
   26 continue
c      if(t(1,1).eq.18.) t(1,1)=t(1,1)/42.
c      write(*,*)mh,mk,ml,mmk,mml
c      write(*,*)mpax,mpz1,mpx1,mpy1,mpy,mpz
      rhomax=-9999.
C
C     calculate FOURIER
C
      do 1 ix=-1,MPAX+1
      ixx = ix+1
      ixx = 1+mod(ixx,3)
      do 3 iik = mmk,mk
      do 4 iz=0,MPZ1
    4 fx(iz)=0.
      do 5 iil = mml,ml
      if(nh(iik,iil).eq.0) goto 5
      jl=iand(iil,MPZ1)
      jl=nmz(jl)
      fs=0.
      do 2 iih=0,mh
      fhkl=f(iih,iik,iil)
      if(fhkl.eq.0.0) goto 2
      ihix=iih*ix
      ihix=iand(ihix,MPX1)
      fs=fs +fhkl*cx(ihix)
    2 continue
C     SUM the (IL+n*MPAR) components into the IL subspace
      fx(jl)=fx(jl)+fs
    5 continue
      call FFT(MPZ,3)
      jk=iand(iik,MPY1)
      do 6 iz=0,MPZ1
C     SUM the(IK+n*MPAR) results into the IK subspace
    6 p(jk,iz) = p(jk,iz) +fx(iz)
    3 continue
      do 7 iz=0,MPZ1
      do 8 jk=0,MPY1
    8 fx(jk)=0.
C     Scan subspace of p(jk,iz)
      do 9 jk=0,MPY1
      kj=nmy(jk)
      fx(kj)=p(jk,iz)
C     Zero P(jk,iz) array before calculating next section on x
    9 p(jk,iz)=0.
      call FFT(MPY,2)
      do 10 iy=0,MPY1
      amax=1000.*fx(iy)
      rhomax=max(rhomax,amax)
      if(amax.ne.rhomax) goto 10
      n11=ixx
      n12=iy
      n13=iz
   10 r(iy,iz,ixx)=fx(iy)
    7 continue
      if(ix.eq.(MPAX+1)) mid=-mid
      call cpu_time(pk0)
      if(ix.ge.1) call pkpick(ix-1,mid,npeaks,rlast)
      call cpu_time(pk1)
      pks=pks+pk1-pk0
      mid=ixx
    1 continue
      if(npeaks.lt.NATOM) write(*,13) npeaks
   13 format(' only ',i3,' peaks noted in map ')
      RETURN
      END
C     
      SUBROUTINE FFT(N,IT)
      PARAMETER (NPAR=256,NPA=NPAR-1)
      COMPLEX F(0:NPA),CX(0:NPA),CY(0:NPA),CZ(0:NPA),A
      COMMON /WORK/ F,CX,CY,CZ
      NX = ALOG(FLOAT(N))/ALOG(2.0) +0.1
      I1=1
      M2=N/2
      DO 1 I=1,NX
      I0=I1
      I2=I1-1
      I1=I1+I1
      IF(IT.EQ.3) GOTO 6
      DO 2 J=0,N-1,I1
      A=F(J+I0)
      F(J+I0)=F(J)-A
      F(J)=F(J)+A
      I4=0
      DO 3 K=J+1,J+I2
      I4=I4+M2
      A = F(K+I0)*CY(I4)
      F(K+I0)=F(K)-A
    3 F(K)=F(K)+A
    2 CONTINUE
      GOTO 1
    6 CONTINUE
      DO 4 J=0,N-1,I1
      A=F(J+I0)
      F(J+I0)=F(J)-A
      F(J)=F(J)+A
      I4=0
      DO 5 K=J+1,J+I2
      I4=I4+M2
      A = F(K+I0)*CZ(I4)
      F(K+I0)=F(K)-A
    5 F(K)=F(K)+A
    4 CONTINUE
    1 M2=M2/2
      RETURN
      END
C
      SUBROUTINE PKPICK(ixm1,MID,npeaks,rlast)
      parameter (NPAR=256)
      real*8 a,b,c,d,e,f,g,ah,ak,al,a1,a2,a3,q1,q2,q3,q4,det
      real*8 t11,t12,t13,t22,t23,t33
      integer*4 nh(27),nk(3)
      dimension r(-1:NPAR,-1:NPAR,3),t(4,4),r27(27),aname(2)
      dimension px(1500,4),rh(1500),ax(3)
      dimension cs(-0:65536),sn(0:65536),xxx(1500,3)
      COMMON /PKS/ r,t,peak,MPX,MPAX,MPY,MPZ,IHORZ,IVERT,ISECT,NATOM
      COMMON /SRT/ px,rh,rholast
      COMMON /SIGNS/ cs,sn,xxx,NLIST
      data nh/0,1,0,1,1,1,0,1,0,4*1,0,4*1,0,1,0,1,1,1,0,1,0/
      MPX1=MPX-1
      MPY1=MPY-1
      MPZ1=MPZ-1
      nstop=mid/iabs(mid)
      mid=iabs(mid)
      ntop=mid+1
      if(ntop.eq.4) ntop=1
      nbot=6-mid-ntop
C     fill out edges of the map
      do 1 ix=1,3
      do 2 iz=0,MPZ1
      r(-1,iz,ix) = r(MPY1,iz,ix)
    2 r(MPY,iz,ix) = r(0,iz,ix)
      do 3 iy=0,MPY1
      r(iy,-1,ix) = r(iy,MPZ1,ix)
    3 r(iy,MPZ,ix) = r(iy,0,ix)
      r(-1,-1,ix) = r(MPY1,MPZ1,ix)
      r(-1,MPZ,ix) = r(MPY1,0,ix)
      r(MPY,-1,ix) = r(0,MPZ1,ix)
      r(MPY,MPZ,ix) = r(0,0,ix)
    1 continue
      if(npeaks.gt.10) goto 12
      pick=0.
      do 4 iz=0,MPZ1
      do 4 iy=0,MPY1
      xx = r(iy,iz,mid)
      if(xx.gt.amax) amax=xx
    4 if(xx.gt.pick) pick=xx
      pick = 0.15*pick
      if(pick.gt.peak) peak=pick
   12 N = npeaks
      rholast=rlast
      if(N.eq.0)open(1,file='FRAG',status='unknown',form='formatted')
      PARX = 1.0/FLOAT(MPX)
      PARY = 1.0/FLOAT(MPY)
      PARZ = 1.0/FLOAT(MPZ)
C     find peaks surrounded by 19 positive points
      ix=mid
      ixm=nbot
      ixp=ntop
      nk(1)=nbot
      nk(2)=mid
      nk(3)=ntop
      do 5 iz=0,MPZ1
      izm=iz-1
      izp=iz+1
      do 6 iy=0,MPY1
      iym=iy-1
      iyp=iy+1
      xx=r(iy,iz,ix)
      zz = 0.00001
      if(xx.lt.peak) goto 6
C     IS THE PEAK A LOCAL MAXIMUM?
      i4=0
      do 50 i3=1,3
      imid=nk(i3)
      do 50 i2=-1,1
      do 50 i1=-1,1
      i4=i4+1
CC    add rho(peak) to all densities to make them all positive?
      r27(i4)=r(iy+i1,iz+i2,imid)
   50 if(nh(i4).eq.1.and.xx.lt.r27(i4)) goto 6
c      write(*,*)0,1,0,zz
C     are the surrounding points all positive densities?
      do 51 i1=1,27
c      if(nh(i1).eq.1) write(*,*)i1,r27(i1)
      if(nh(i1).eq.1.and.zz.ge.r27(i1)) goto 7
   51 CONTINUE
C     perform Rollet's ellipsoidal fit procedure
      oqq = log(r27(10))
      qqo = log(r27(4))
      oqo = log(r27(13))
      pqo = log(r27(22))
      oqp = log(r27(16))
      qoq = log(r27(2))
      ooq = log(r27(11))
      poq = log(r27(20))
      qoo = log(r27(5))
      ooo = log(xx)
      poo = log(r27(23))
      qop = log(r27(8))
      oop = log(r27(17))
      pop = log(r27(26))
      opq = log(r27(12))
      qpo = log(r27(6))
      opo = log(r27(15))
      ppo = log(r27(24))
      opp = log(r27(18))
      b = (poo+ppo+pqo+pop+poq -qoo-qqo-qpo-qoq-qop)/10.
      c = (opo+ppo+qpo+opp+opq -oqo-qqo-pqo-oqq-oqp)/10.
      d = (oop+opp+oqp+pop+qop -ooq-oqq-opq-qoq-poq)/10.
      ah = (opp+oqq-opq-oqp)/4.
      ak = (pop+qoq-poq-qop)/4.
      al = (ppo+qqo-pqo-qpo)/4.
      a1 = ppo+pqo+qpo+qqo
      a2 = pop+poq+qop+qoq
      a3 = opp+opq+oqp+oqq
      q1 = ooo+poo+qoo+opo+oqo+oop+ooq+a1+a2+a3
      q2 = poo+qoo+a1+a2
      q3 = opo+oqo+a1+a3
      q4 = oop+ooq+a2+a3
      a = t(1,1)*q1 + t(1,2)*q2 + t(1,3)*q3 + t(1,4)*q4
      e = t(2,1)*q1 + t(2,2)*q2 + t(2,3)*q3 + t(2,4)*q4
      f = t(3,1)*q1 + t(3,2)*q2 + t(3,3)*q3 + t(3,4)*q4
      g = t(4,1)*q1 + t(4,2)*q2 + t(4,3)*q3 + t(4,4)*q4
      det = -8*e*f*g -2*ah*ak*al + 2*f*ak*ak +2*ah*ah*e +2*g*al*al
      adet=det
      if(abs(adet).lt.0.0001) goto 7
      t11 = (4*f*g-ah*ah)/det
      t12 = (ah*ak-2*g*al)/det
      t13 = (ah*al-2*f*ak)/det
      t22 = (4*e*g-ak*ak)/det
      t23 = (ak*al-2*e*ah)/det
      t33 = (4*e*f-al*al)/det
      x = t11*b + t12*c + t13*d
      y = t12*b + t22*c + t23*d
      z = t13*b + t23*c + t33*d
      ax(1) = PARX*(float(ixm1)+x)
      if(ax(1).lt.-0.0001.or.ax(1).gt.1.0001) goto 6
      N=N+1
      if(N.gt.NATOM) N=NATOM+1
      rh(N) =dexp(a+b*x+c*y+d*z+e*x*x+f*y*y+g*z*z+ah*z*y+ak*x*z+al*x*y)
      ax(2) = PARY*(float(iy)+y)
      ax(3) = PARZ*(float(iz)+z)
      px(n,4) = 0.1
      goto 8
    7 N=N+1
      if(N.gt.NATOM) N=NATOM+1
      ax(1) = PARX*float(ixm1)
      ax(2) = PARY*float(iy)
      ax(3) = PARZ*float(iz)
      rh(N) = xx
      px(N,4) = 1.1
    8 continue
      px(N,1)=ax(IHORZ)
      px(N,2)=ax(IVERT)
      px(N,3)=ax(ISECT)
      if(N.eq.NATOM) call SORT(N)
      if(N.eq.(NATOM+1).and.rh(N).gt.rholast) call INSERT(NATOM)
      npeaks=N
      rlast=rholast
      peak=0.85*rholast
    6 continue
    5 continue
      if(nstop.eq.1) RETURN
      if(N.lt.NATOM.and.N.gt.0) call SORT(N)
      if(N.gt.NATOM) N = NATOM
      do 11 i=1,N
      xxx(i,1)=px(i,1)
      xxx(i,2)=px(i,2)
      xxx(i,3)=px(i,3)
      ichar=px(i,4)
      IF(NLIST.EQ.2)write(99,98)px(i,1),px(i,2),px(i,3),rh(i),i,ichar
   98 format(3f9.5,f10.2,2i5)
   11 CONTINUE
      RETURN
      END
C
      subroutine SORT(N)
C     sort the peaks in descending order of magnitude
      integer*4 num(1500)
      dimension px(1500,4),rh(1500),q(1500,4)
      COMMON /SRT/ px,rh,rholast
      do 1 i=1,N
    1 num(i)=i
      do 2 i=2,N
      j=i
      hold = rh(1)
      rh(1) = rh(j)
      do while (rh(j) .gt. rh(j-1))
         x=rh(j)
         rh(j)=rh(j-1)
         rh(j-1)=x
         ix=num(j)
         num(j)=num(j-1)
         num(j-1)=ix
         j = j-1
      end do
      rh(1) = hold
      if (rh(2) .gt. rh(1)) then
         x=rh(2)
         rh(2)=rh(1)
         rh(1)=x
         ix=num(2)
         num(2)=num(1)
         num(1)=ix
      end if
    2 continue
      do 3 i=1,N
      k=num(i)
      q(i,1)=px(k,1)
      q(i,2)=px(k,2)
      q(i,3)=px(k,3)
    3 q(i,4)=px(k,4)
      do 4 i=1,N
      px(i,1)=q(i,1)
      px(i,2)=q(i,2)
      px(i,3)=q(i,3)
    4 px(i,4)=q(i,4)
      rholast = rh(N)
      RETURN
      END
C
      subroutine INSERT(NATOM)
C     place newest peak within the sorted stack
      dimension px(1500,4),rh(1500)
      COMMON /SRT/ px,rh,rholast
      NATOM1=NATOM+1
      rho = rh(NATOM1)
      do 1 i=1,NATOM-1
      j=NATOM-i
      j1=j+1
      if(rh(j).gt.rho)  goto 2
      rh(j1)=rh(j)
      do 4 k=1,4
    4 px(j1,k)=px(j,k)
    1 continue
      rh(1)=rho
      do 5 i=1,4
    5 px(1,i)=px(NATOM1,i)
      goto 6
    2 rh(j1)=rh(NATOM1)
      do 3 j=1,4
    3 px(j1,j)=px(NATOM1,j)
    6 rholast = rh(NATOM)
      RETURN
      END
C
      SUBROUTINE SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
      integer*4 jh,jk,jl,nv(3)
      nv(IHORZ) = jh
      nv(IVERT) = jk
      nv(ISECT) = jl
      jh = nv(1)
      jk = nv(2)
      jl = nv(3)
      RETURN
      END
C
