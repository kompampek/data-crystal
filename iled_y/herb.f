C     PKPKing CC mode
      parameter mp=8000,mt=400000,mp2=16001,mt3=3*mt
      integer*2 IPR(9,3,3),narg,nar
      integer*2 ih(mp),ik(mp),il(mp),ns,ne,nr(mp),nq(6,1000)
      integer*2 na,m1,m2,n1,n2,n3,ix,iy,nn(3,mt),m(mp)
      integer*4 iseed,mm(0:mp),n0,nm(mt3)
      real*8 x,y,z,rand
      dimension p(-mp:mp),tt(mt),an(mt),e(mp),TR(9,3),pp(-mp:mp)
      dimension tnt(mt),ANAME(2),rr(4)
      dimension cs(-32768:32768),sn(-32768:32768),xx(900,3)
      COMMON /SIGNS/ cs,sn,xx
      COMMON /HKL/ MREF,ih,ik,il,e,p,NSYM,tr,ipr
      data m,mm/mp2*0/
      open(unit=2,file='fft.dat',status='old',form='formatted')
      read(2,*)NSYM
      read(2,*)NSYM
      do 4 i=1,NSYM
    4 read(2,19) (tr(i,j),(ipr(i,j,k),k=1,3),j=1,3)
   19 format(3(f10.5,3i5))
      close(unit=2)
      pi = acos(-1.0)
      pi2=2.*pi
      arg=pi/32768.
      do 43 i=-32768,32768
      ar=i*arg
      sn(i)=sin(ar)
   43 cs(i)=cos(ar)
      pi6=pi/6.
      pion2=pi/2.
      write(*,48)
   48 format(' Enter 4-letter registry name : ',$)
      read(*,49)ANAME(1)
   49 format(A4,x,2f8.3)
      ANAME(2)='.amb'
      open(unit=3,file=ANAME,status='unknown',form='formatted')
      write(*,11)
   11 format(' Odd-RAN#,NAT,#PS-Cycs,#REF-Cycs,#SETS,LIST,MAT :',$)
C     ix odd > 0 to scan random phase sets
      read(*,*)id,NATOM,iter,mcyc,msets,list,MATOM
    8 format(i7,9i5)
      if(id.eq.0) open(33,file='fort.333.old',form='unformatted')
      if(id.lt.0) open(33,file='IRAN-NOS')
      ix=abs(id)
      iran=ix
      ncyc=0
C     READ E-DATA FILE
      ANAME(2)='.eee'
      open(unit=1,file=ANAME,status='old',form='unformatted')
      do 1 i=1,mp
      read(1,end=47)ih(i),ik(i),il(i),ns,ne,fo,sf,fc,nar
      pp(i)=0.001*nar
      pp(-i)=-pp(i)
      nr(i)=13
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
   47 close(unit=1)
      ix=iabs(ix)
C     READ TRIPLES FILE
      ANAME(2)='.trp'
      open(unit=1,file=ANAME,status='old',form='unformatted')
      NREF=0
      bot=0.
      trip=0.
      ntrip=0
      qs=0
      do 33 i=1,mt
      read(1,end=2)na,m1,m2,n1,n2,n3
      a1=0.001*na
      it=-1+(m1/256).and.15
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
      a4=cos(pp(n1)+pp(n2)+pp(n3))
      if(it.eq.6) a4=-a4
      qs=qs+a1*(a4-T)**2
      tnt(ntrip)= it*pi6
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
      qs=qs/bot
      write(*,*)qs
      close(unit=1)
      MREF=NREF
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
      write(*,36)NREF,trip,ntrip,qest,mq
   36 format(i5,' phs, Cosav =',f6.3,' for ',i6,' trips, Nqest =',f7.3,'
     * for ',i4,' quarts'//'   Amb# Ncyc  Rmin  Cosav   Rcry '/)
      t0=secnds(0.0)
      t2=0.0
      t4=0.0
C
C     LOOP OVER THE NUMBER OF DIFFERENT TRIALS
C
      do 15 nsets=1,msets
      if(id.ne.0) goto 157
      read(33)iran
      write(*,*)iran
      read(33) (p(i),i=1,NREF)
      do 155 i=1,NREF
  155 p(-i)=-p(i)
  157 NLIST=0
      if(id.lt.0) read(33,*)iran
      if(nsets.eq.msets) NLIST=1
      if(id.eq.0) goto 156
      MAD=NATOM
C
C     GENERATE RANDOM ATOMS or read XYZ file
C
      IF(MATOM.eq.0) GOTO 46
      ANAME(2)='.xyz'
      open(unit=4,file=ANAME,status='old',form='formatted')
      NEND = MATOM*(NSETS-1)
      DO 28 I=1,NEND
   28 read(4,*)X
      DO 26 I=1,MATOM
   26 read(4,*)xx(i,1),xx(i,2),xx(i,3)
      CLOSE(4)
      MAD=MATOM
      GOTO 30
   46 ix=iran
      do 7 i=1,NATOM
      iy=899*ix
      if(iy.lt.0) iy=iy+32768
      if(ix.eq.0) goto 7
      iseed=iy
      call srand(iseed)
      x=rand()
      y=rand()
      z=rand()
      xx(i,1)=x
      xx(i,2)=y
      xx(i,3)=z
    7 ix=iy
   30 call SFCLC(MAD,RES)
      close(unit=1)
  156 CONTINUE
C
C     THIS IS THE MAIN CYCLES LOOP
C
      do 10 ncyc=1,mcyc
      do 21 nit=0,iter
      if(nit.eq.0) goto 20
C     PARMETER SHIFT - INNER LOOP 16
      do 16 i=1,NREF
      nrk=nr(i)
      do 166 j=1,4
  166 rr(j)=0.
      do 18 ij=mm(i-1)+1,mm(i)
      k=nm(ij)
      n1=nn(1,k)
      n2=nn(2,k)
      n3=nn(3,k)
      a1=an(k)
      tau=tt(k)
      narg = (p(n1)+p(n2)+p(n3)+tnt(k))/arg
      if(i.eq.-n2.or.i.eq.-n3) narg=-narg
      do 188 j=1,4
      narg=narg+16384
  188 rr(j)=rr(j)+a1*(cs(narg)-tau)**2
   18 CONTINUE
      rmi=99999.
      do 17 j=1,4
      narg=narg+16384
      if(nrk.eq.13) goto 133
      if(j.eq.1.or.j.eq.3) goto 17
  133 if(rr(j).gt.rmi) goto 17
      rmi=rr(j)
      ind=j
   17 continue
      p(i)=p(i)+ind*pion2
      p(-i)=-p(i)
   16 continue
C     CALCULATE RMIN
      goto 24
   20 continue
      Rold=0.
      Cosav=0.
      do 3 i=1,ntrip
      a1=an(i)
      narg=(p(nn(1,i))+p(nn(2,i))+p(nn(3,i))+tnt(i))/arg
      cost=cs(narg)
      Cosav=Cosav+a1*cost
      cost=tt(i)-cost
    3 Rold=Rold + a1*cost*cost
      Rold=Rold/bot
      Cv = Cosav/bot
      Cosa=Cosav
      kcyc=ncyc
   24 continue
   21 continue
      t1=secnds(0.0)
      call EDMAP(AME,NATOM,NLIST,ppktime)
      ppk=ppk+ppktime
      t1=secnds(t1)/60.
      if(t1.lt.0.0) t1=t1+1440.
      t4=t4+t1
      t1=secnds(0.0)
      call SFCLC(NATOM,RES)
      rms0=999.
      do 27 j=0,15
      i0=j
      i1=i0.and.1
      i0=i0/2
      i2=i0.and.1
      i0=i0/2
      i3=i0.and.1
      i0 = i0/2
      i4 = 2*i0 -1
      rms=0.
      do 25 i=1,NREF
      i0=abs(ih(i)*i1+ik(i)*i2+il(i)*i3).and.1
      dp=mod(abs(p(i)-i4*pp(i)+i0*PI),PI2)
      if(dp.gt.PI) dp=PI2-dp
   25 rms=rms+dp*dp
      rms=57.296*sqrt(rms/NREF)
   27 rms0=min(rms0,rms)
      if(list.ne.0) write(*,5)iran,ncyc,Rold,Cv,RES,rms0
    5 format(i7,i4,3f7.3,f6.1)
      t1=secnds(t1)/60.
      if(t1.lt.0.0) t1=t1+1440.
      t2=t2+t1
   10 continue
   23 jcyc=kcyc
      if(jcyc.le.0) jcyc=1
      write(3,44)iran,kcyc,Rold,Cv,RES,rms0
   44 format(i6,i4,3f7.3,f6.1)
      write(333)iran,kcyc,Rold,Cv,RES
      write(333) (p(i),i=1,NREF)
      iran=iran+2
   15 continue
      do 9 i=1,MREF
      ns=1000*p(i)
    9 write(7) ih(i),ik(i),il(i),ns,ne,e(i),e(i),e(i),ns
      t3=secnds(t0)/60.
      if(t3.lt.0.0) t3=t3+1440.
      ppk=ppk/60.
      if(ppk.lt.0.0) ppk=ppk+1440.
      t5=t3-t4-t2
      t4=t4-ppk
      if(list.ne.0) write(*,6)t3,t5,t2,t4,ppk
      write(3,6)t3,t5,t2,t4,ppk
    6 format(' Total time =',f9.3,' Minutes'/' PS =',f9.3,' SF ='f9.3,' 
     *FR =',f9.3,' PPK ='f9.3,' Minutes')
      close(unit=3)
      end
C
      SUBROUTINE SFCLC(NAT,RES)
      parameter mp=8000
      integer*2 jh(mp),jk(mp),jl(mp),IPR(9,3,3)
      integer*2 ih,ik,il,nx(8,499),ny(8,499),nz(8,499),narg
      dimension e(mp),p(-mp:mp),TR(9,3)
      dimension cs(-32768:32768),sn(-32768:32768),xxx(900,3)
      COMMON /SIGNS/ cs,sn,xxx
      COMMON /HKL/ NREF,jh,jk,jl,e,p,NSYM,tr,ipr
      do 1 i=1,NAT
      xx=xxx(i,1)
      yy=xxx(i,2)
      zz=xxx(i,3)
      do 2 j=1,NSYM
      x = ipr(j,1,1)*xx +ipr(j,1,2)*yy +ipr(j,1,3)*zz +tr(j,1)
      nx(j,i)=65536.*amod(x,1.0)
      y = ipr(j,2,1)*xx +ipr(j,2,2)*yy +ipr(j,2,3)*zz +tr(j,2)
      ny(j,i)=65536.*amod(y,1.0)
      z = ipr(j,3,1)*xx +ipr(j,3,2)*yy +ipr(j,3,3)*zz +tr(j,3)
    2 nz(j,i)=65536.*amod(z,1.0)
    1 continue
      top=0
      bot=0
      do 4 i=1,NREF
      ih=jh(i)
      ik=jk(i)
      il=jl(i)
      aa=0.
      bb=0.
      do 3 j=1,NAT
      do 6 k=1,NSYM
      narg = ih*nx(k,j) + ik*ny(k,j) + il*nz(k,j)
      aa = aa + cs(narg)
    6 bb = bb + sn(narg)
    3 continue
      sig=sqrt(float(NSYM*NAT))
      ec=sqrt(aa*aa+bb*bb)/sig
      top=top+abs(ec-e(i))
      bot=bot+e(i)
    5 format(i5,2f8.3)
      p(i)=atan2(bb,aa)
    4 p(-i)=-p(i)
      RES=top/bot
      RETURN
      end
C
C     TWO DIMENSIONAL FFT
      SUBROUTINE EDMAP(AME,MATOM,LIST,ppk)
      parameter (NPAR=128,NPA=NPAR-1,NP=84,mp=8000)
      byte nh(-NP:NP,-NP:NP),nt(16,9),HH,HK,HL
      real*8 arg,pi2,px
      integer*2 ih(mp),ik(mp),il(mp),jh,jk,jl,nm(0:NPA),IPR(9,3,3)
      dimension r(-1:NPAR,-1:NPAR,-1:NPAR),tx(16,3),p(-mp:mp),e(mp)
      dimension TR(9,3)
      complex f(0:NP,-NP:NP,-NP:NP),cx(0:NPA),fx(0:NPA)
      complex pp(0:NPA,0:NPA),fhkl,fs
      COMMON /WORK/ fx,cx,n1,NX
      COMMON /PKS/ r,peak,MPAR,MPA,MPAX,IHORZ,IVERT,ISECT,NATOM
      COMMON /HKL/ NREF,ih,ik,il,e,p,NSYM,tr,ipr
      COMMON /SYM/ nt,tx,msym
      data mh,mk,ml,mmk,mml/5*0/
      NATOM=MATOM
      px=-1.0
      open(unit=1,file='fft.dat',status='old',form='formatted')
      read(1,*)MPAR,MPAX,IHORZ,IVERT,ISECT
      read(1,*)nsym
      msym=nsym
      do 16 i=1,nsym
      read(1,*)tx(i,1),HH,HK,HL,tx(i,2),KH,KK,KL,tx(i,3),LH,LK,LL
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
      do 22 i=-1,MPAX
      do 22 j=-1,MPAR
      do 22 k=-1,MPAR
   22 r(k,j,i)=0.
      MPA=MPAR-1
      pi2=2.0*dacos(px)
      px=pi2/MPAR
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPA
      arg=-i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     set up FFT arrays for bit reversal on first pass
      n1=MPAR-1
      nm(0)=0
      nm(1)=1
      NX = alog(float(MPAR))/alog(2.0) - 0.9
      mm=1
      do 23 i=1,NX+1
      m2=mm-1
      do 24 j=0,m2
      nm(j)=2*nm(j)
   24 nm(j+mm) = nm(j)+1
      mm=mm+mm
   23 continue
C     read reflection data
      do 26 i=1,NREF
      ihh=ih(i)
      ikk=ik(i)
      ill=il(i)
      ff=e(i)
      phase=p(i)
      do 19 j=1,nsym
      jh = ihh*nt(j,1) + ikk*nt(j,2) + ill*nt(j,3)
      jk = ihh*nt(j,4) + ikk*nt(j,5) + ill*nt(j,6)
      jl = ihh*nt(j,7) + ikk*nt(j,8) + ill*nt(j,9)
      st = ihh*tx(j,1) + ikk*tx(j,2) + ill*tx(j,3)
      st = mod(st,1.0)
      ph=phase + st*pi2
      ph=mod(ph,6.283185)
      call SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
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
   19 f(jh,jk,jl) = ff*cmplx(cos(ph),sin(ph))
   26 continue
   27 continue
C     calculate FOURIER
      do 1 ix=-1,MPAX
      do 3 ikk = mmk,mk
      do 4 iz=0,MPA
    4 fx(iz)=0.
      do 5 ill = mml,ml
      if(nh(ikk,ill).eq.0) goto 5
      fs=0.
      do 2 ihh=0,mh
      fhkl=f(ihh,ikk,ill)
      if(fhkl.eq.0.0) goto 2
      fs=fs +fhkl*cx((ihh*ix).and.MPA)
    2 continue
      fx(nm(ill.and.MPA))=fs
    5 continue
      call FFT(MPAR)
      jk=ikk.and.MPA
      do 6 iz=0,MPA
    6 pp(jk,iz) = fx(iz)
    3 continue
      do 7 iz=0,MPA
      do 8 ikk=0,MPA
    8 fx(ikk)=0.
      do 9 ikk = mmk,mk
      jk=ikk.and.MPA
    9 fx(nm(jk))=pp(jk,iz)
      call FFT(MPAR)
      do 10 iy=0,MPA
   10 r(iy,iz,ix)=fx(iy)
    7 continue
    1 continue
C     end of FOURIER loops
      amax=0.
      amin=0.
      do 12 ix=-1,MPAX
      do 12 iz=0,MPA
      do 12 iy=0,MPA
      ama=r(iy,iz,ix)
      if(ama.lt.amin) amin=ama
   12 if(ama.gt.amax) amax=ama
      peak = 0.25*amax
      ppk=secnds(0.0)
      call PKPICK(AME,LIST)
      ppk=secnds(ppk)
      RETURN
      end
C     
      SUBROUTINE FFT(N)
      parameter (NPAR=128,NPA=NPAR-1)
      complex f(0:NPA),cx(0:NPA),a,b
      COMMON /WORK/ f,cx,n1,NX
      do 1 i1 = 0,n1,2
      a= f(i1)
      b= f(i1+1)
      f(i1)   = a + b
    1 f(i1+1) = a - b
      i1=2
      m2=N/4
      do 2 j0 = 1,NX
      i0=i1
      i2=i1-1
      i1=i1+i1
      do 3 i=0,n1,i1
      i4=-m2
      do 4 k=i,i+i2
      i4=i4+m2
      a = f(k)
      kk=k+i0
      b = f(kk)*cx(i4)
      f(k)=a+b
    4 f(kk)=a-b
    3 continue
      m2=m2/2
    2 continue
      RETURN
      end
C
      SUBROUTINE PKPICK(AME,LIST)
      parameter (NPAR=128,NPA=NPAR-1,NP=84)
      byte nt(16,9),jt(9)
      real*8 AZ,A0,A1,A2,A3,AA1,AA2,AA3
      dimension r(-1:NPAR,-1:NPAR,-1:NPAR),tx(16,3)
      dimension px(900,5),ax(3),aname(2),rr(5)
      dimension cs(-32768:32768),sn(-32768:32768),xxx(900,3)
      COMMON /PKS/ r,peak,MPAR,MPA,MPAX,IHORZ,IVERT,ISECT,NATOM
      COMMON /SRT/ px,rholast
      COMMON /SIGNS/ cs,sn,xxx
      COMMON /SYM/ nt,tx,nsym
      data good,bad,bname/'   G','   B','.pks'/
C     fill out edges of the map
      do 1 ix=-1,MPAX
      do 2 iz=0,MPA
      r(-1,iz,ix) = r(MPA,iz,ix)
    2 r(MPAR,iz,ix) = r(0,iz,ix)
      do 3 iy=0,MPA
      r(iy,-1,ix) = r(iy,MPA,ix)
    3 r(iy,MPAR,ix) = r(iy,0,ix)
      r(-1,-1,ix) = r(MPA,MPA,ix)
      r(-1,MPAR,ix) = r(MPA,0,ix)
      r(MPAR,-1,ix) = r(0,MPA,ix)
      r(MPAR,MPAR,ix) = r(0,0,ix)
    1 continue
      N = 0
      rholast=0.0
      PAR = 1.0/FLOAT(MPAR)
      MPA_2 = MPAR/2 -1
C     find peaks surrounded by 19 positive points
      do 4 ix=0,MPAX-1
      ixm=ix-1
      ixp=ix+1
      do 5 iz=0,MPA
      izm=iz-1
      izp=iz+1
      do 6 iy=0,MPA
      iym=iy-1
      iyp=iy+1
      x=r(iy,iz,ix)
      if(x.lt.peak) goto 6
C     is the peak a local maximum
      f5  = r(iy,izm,ix)
      f11 = r(iym,iz,ix)
      f13 = r(iy,iz,ixm)
      f15 = r(iy,iz,ixp)
      f17 = r(iyp,iz,ix)
      if(x.lt.f5.or.x.lt.f11.or.x.lt.f13.or.x.lt.f15.or.x.lt.f17) goto 6
      f2  = r(iym,izm,ix)
      f22 = r(iy,izp,ixm)
      f23 = r(iy,izp,ix)
      f24 = r(iy,izp,ixp)
      f26 = r(iyp,izp,ix)
      if(x.lt.f2.or.x.lt.f22.or.x.lt.f23.or.x.lt.f24.or.x.lt.f26) goto 6
      f4  = r(iy,izm,ixm)
      f6  = r(iy,izm,ixp)
      f8  = r(iyp,izm,ix)
      f10 = r(iym,iz,ixm)
      f12 = r(iym,iz,ixp)
      if(x.lt.f4.or.x.lt.f6.or.x.lt.f8.or.x.lt.f10.or.x.lt.f12) goto 6
      f16 = r(iyp,iz,ixm)
      f18 = r(iyp,iz,ixp)
      f20 = r(iym,izp,ix)
      if(x.lt.f16.or.x.lt.f18.or.x.lt.f20) goto 6
C     are the surrounding points all positive densities
      z = 0.0001
      if(z.ge.f2.or.z.ge.f22.or.z.ge.f23.or.z.ge.f24.or.z.ge.f26) goto 7
      if(z.ge.f4.or.z.ge.f16.or.z.ge.f17.or.z.ge.f18.or.z.ge.f20) goto 7
      if(z.ge.f6.or.z.ge.f11.or.z.ge.f12.or.z.ge.f13.or.z.ge.f15) goto 7
      if(z.ge.f5.or.z.ge.f8.or.z.ge.f10) goto 7
C     perform Pavelcik's ellipsoidal fit procedure
      f14 = x
      AZ = f2+f4+f6+f8+f10+f12+f16+f18+f20+f22+f24+f26
      A0 = (-AZ +4.*(f5+f11+f13+f15+f17+f23+3*f14))/24.
      A1 = (f6-f4+f12-f10+f18-f16+f24-f22 +4*(f15-f11))/16.
      A2 = (f8-f2+f16-f10+f18-f12+f26-f20 +4*(f17-f11))/16.
      A3 = (f22-f4+f20-f2+f26-f8+f24-f6 +4*(f23-f5))/16.
      AA1 = (AZ-f2-f8-f20-f26 +2*(f13+f15-f5-f11-f17-f23 -2*f14))/12.
      AA2 = (AZ-f4-f6-f22-f24 +2*(f11+f17-f5-f13-f15-f23 -2*f14))/12.
      AA3=(AZ-f10-f12-f16-f18 +2*(f5+f23-f11-f13-f15-f17 -2*f14))/12.
      x = -0.5*A1/AA1
      y = -0.5*A2/AA2
      z = -0.5*A3/AA3
      if(abs(x).gt.1.0.or.abs(y).gt.1.0.or.abs(z).gt.1.0) goto 7
      ax(1) = PAR*(float(ix)+x)
      if(ax(1).lt.0.0) goto 6
      N=N+1
      if(N.gt.NATOM) N=NATOM+1
      px(N,5) = A0 - 0.25*(A1*A1/AA1 + A2*A2/AA2 + A3*A3/AA3)
      ax(2) = PAR*(float(iy)+y)
      ax(3) = PAR*(float(iz)+z)
      char = good
      goto 8
    7 N=N+1
      if(N.gt.NATOM) N=NATOM+1
      ax(1) = PAR*float(ix)
      ax(2) = PAR*float(iy)
      ax(3) = PAR*float(iz)
      px(N,5) = x
      char = bad
    8 px(N,1)=ax(IHORZ)
      px(N,2)=ax(IVERT)
      px(N,3)=ax(ISECT)
      px(N,4)=char
      if(N.eq.NATOM) call SORT(N)
      if(N.eq.(NATOM+1).and.px(N,5).gt.rholast) call INSERT(N)
      peak=0.85*rholast
    6 continue
    5 continue
    4 continue
      if(N.lt.NATOM) call SORT(N)
      if(N.gt.NATOM) N = NATOM
      do 12 i=1,N
      xxx(i,1)=px(i,1)
      xxx(i,2)=px(i,2)
   12 xxx(i,3)=px(i,3)
C      if(LIST.eq.0) return
      scale=1000/px(1,5)
      do 9 i=1,N
      ix = (scale*px(i,5) + 0.5)
    9 write(99,10) (px(i,j),j=1,4),ix,i
   10 format(3f9.5,5x,a4,2i5)
      CLOSE(UNIT=99)
      RETURN
      END
C
      subroutine SORT(NATOM)
C     sort the peaks in descending order of magnitude
      integer*2 num(900)
      dimension px(900,4),rh(900),q(900,4)
      COMMON /SRT/ px,rh,rholast
      N=NATOM-1
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
      end
C
      subroutine INSERT(M)
C     place newest peak within the sorted stack, M = NATOM+1
      dimension px(900,4),rh(900)
      COMMON /SRT/ px,rh,rholast
      NATOM1=M
      NATOM =M-1
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
    2 rh(j1)=rho
      do 3 j=1,4
    3 px(j1,j)=px(NATOM1,j)
    6 rholast = rh(NATOM)
      RETURN
      end
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
      end
