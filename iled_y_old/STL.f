C     P212121 specific for GRAMA
      parameter mp=2000,mt=40000,mp2=4001,mt3=3*mt
      integer*2 IPR(9,3,3),narg
      integer*2 ih(mp),ik(mp),il(mp),ns,ne,nr(mp),nq(6,1000)
      integer*2 na,m1,m2,n1,n2,n3,ix,iy,nn(3,mt),m(mp)
      integer*4 iseed,mm(0:mp),n0,nm(mt3)
      real*8 x,y,z,rand
      dimension p(-mp:mp),tt(mt),an(mt),e(mp),TR(9,3)
      dimension tnt(mt),ANAME(2),BNAME(7),r(4),true(mp)
      dimension cs(-32768:32768),sn(-32768:32768),xx(1000,3)
      COMMON /SIGNS/ cs,sn,xx,NLIST
      COMMON /HKL/ NREF,ih,ik,il,e,p,NSYM,tr,ipr
      data m,mm/mp2*0/
      data BNAME/'.nqt','.stl','.trp','.pks','.xyz','.sfc','.amb'/
      open(unit=2,file='fft0.dat',status='old',form='formatted')
      read(2,*)NSYM
      read(2,*)NSYM,VOL
      do 4 i=1,NSYM
    4 read(2,19) (tr(i,j),(ipr(i,j,k),k=1,3),j=1,3)
   19 format(3(f10.5,3i5))
      close(unit=2)
      pi = acos(-1.0)
      pi2=2.*pi
      arg=pi/32768.
      do 32 i=-32768,32768
      ar=i*arg
      sn(i)=sin(ar)
   32 cs(i)=cos(ar)
      pi6=pi/6.
      pion2=pi/2.
      write(*,27)
   27 format(' Enter 4-letter registry name : ',$)
      read(*,25)AME,escape
   25 format(A4,x,f8.3)
      ANAME(1)=AME
      ANAME(2)=BNAME(7)
      open(unit=3,file=ANAME,status='unknown',form='formatted')
      write(*,11)
   11 format(' IRAN#,NAT,#PScycs,#REFcycs,#SETS,LIST,MAT,NSTART :',$)
C     ix odd > 0 to scan random phase sets
      read(*,*)id,NATOM,iter,mcyc,msets,list,MATOM,NSTART
      write(3,24)id,NATOM,iter,mcyc,msets,list,MATOM,NSTART
   24 format(8i5)
      if(id.lt.0) open(33,file='IRAN-NOS')
      ix=abs(id)
      iran=ix
C     READ NEGATIVE QUARTETS
      ANAME(2)=BNAME(1)
      open(unit=1,file=ANAME,status='old',form='unformatted')
      qsum=0.
      quest=0.
      do 35 i=1,1000
      read(1,end=37)nq(1,i),ns,ne,nq(2,i),nq(3,i),nq(4,i),nq(5,i)
      a1=0.002*nq(1,i)
      qsum=qsum+nq(1,i)
      it=(ns/256).and.15
      nq(6,i)=it-1
      t=0.
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
   35 quest=quest + t*nq(1,i)
   37 close(unit=1)
      qest=-quest/qsum
      mq=i-1
      ncyc=0
C     READ SYM FILE
      ANAME(2)=BNAME(2)
      open(unit=1,file=ANAME)
      do 1 i=1,mp
      read(1,*,end=29)ih(i),ik(i),il(i),e(i),true(i),stl
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
      esq=esq+e(i)*e(i)
    1 e(i)=2.*e(i)/VOL
   29 close(unit=1)
      var=NSYM*esq/(VOL*VOL)
C     READ TRIPLES FILE
      ANAME(2)=BNAME(3)
      open(unit=1,file=ANAME,status='old',form='unformatted')
      NREF=0
      bot=0.
      trip=0.
      ntrip=0
      do 33 i=1,mt
      read(1,end=2)na,m1,m2,n1,n2,n3
      it=m1/4096
      if(it.eq.7) goto 33
      ntrip=ntrip+1
      nn(1,ntrip)=n1
      nn(2,ntrip)=n2
      nn(3,ntrip)=n3
      a1=0.001*na
      t=0.
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      if(it.eq.1) T = tanh(0.5*a1)
      it=(m1/256).and.15
      tnt(ntrip)= (it-1)*pi6
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
     * for ',i4,' quarts'//'   Amb# Ncyc  Rmin  Cosav  Nqest '/)
      t0=secnds(0.0)
C
C     LOOP OVER THE NUMBER OF DIFFERENT TRIALS
C
      do 15 nsets=1,msets
      if(id.lt.0) read(33,*)iran
      NLIST=0
      if(nsets.eq.msets) NLIST=1
      MAD=NSTART
C
C     GENERATE RANDOM ATOMS or read XYZ file
C
      IF(MATOM.eq.0) GOTO 31
      ANAME(2)=BNAME(5)
      open(unit=4,file=ANAME,status='old',form='formatted')
      NEND = MATOM*(NSETS-1)
      DO 28 I=1,NEND
   28 read(4,*)X
      DO 26 I=1,MATOM
   26 read(4,*)xx(i,1),xx(i,2),xx(i,3)
      CLOSE(4)
      MAD=MATOM
      GOTO 30
   31 ix=iran
      do 7 i=1,NSTART
      iy=899*ix
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
   30 call SFCLC(MAD)
      close(unit=1)
C
C     THIS IS THE MAIN CYCLES LOOP
C
      do 10 ncyc=1,mcyc
      do 21 nit=0,iter
      if(nit.eq.0) goto 20
C     PARMETER SHIFT - INNER LOOP 16
      do 16 i=1,NREF
      do 17 j=1,4
   17 r(j)=0.
      do 18 ij=mm(i-1)+1,mm(i)
      k=nm(ij)
      ttk=tt(k)
      narg = (p(nn(1,k))+p(nn(2,k))+p(nn(3,k))+tnt(k))/arg
      if(i.eq.-n2.or.i.eq.-n3) narg=-narg
      ank=an(k)
      do 9 j=1,4
      narg=narg+16384
      cost=cs(narg)-ttk
    9 r(j)=r(j)+ank*cost*cost
   18 CONTINUE
      rmin=r(4)
      ind=0
      do 8 j=1,3
      if(r(j).gt.rmin) goto 8
      rmin=r(j)
      ind=j
    8 continue
      p(i)=p(i)+ind*pion2
      p(-i)=-p(i)
   16 continue
      goto 21
C     CALCULATE RMIN
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
      Cosav = Cosav/bot
      quest=0.
      do 34 i=1,mq
      narg=(p(nq(2,i))+p(nq(3,i))+p(nq(4,i))+p(nq(5,i))+nq(6,i)*pi6)/arg
   34 quest=quest + nq(1,i)*cs(narg)
      quest=quest/qsum
      kcyc=ncyc
      if(list.ne.0) write(*,5)iran,ncyc,Rold,Cosav,quest
    5 format(i7,i4,3f7.3)
   21 continue
      if(NLIST.eq.1.and.quest.lt.escape) NLIST=2
      if(NLIST.eq.1.and.ncyc.eq.mcyc) NLIST=2
      t1=secnds(0.0)
      tp=0.
      call EDMAP(tp,VAR)
      tp=tp/60.
      tF=secnds(t1)/60. -tp
      tFs=tFs+tF
      tPks=tPks+tp
      t1=secnds(0.0)
      call SFCLC(NATOM)
      tSF=tSF+secnds(t1)/60.
      if(quest.lt.escape) goto 23
   10 continue
   23 write(3,5)iran,kcyc,Rold,quest,Cosav
      iran=iran+2
   15 continue
      tsum=secnds(t0)/60.
      if(tsum.lt.0.0) tsum=tsum+1440.
      tPS=tsum-tFs-tPks-tSF
      write(*,6)tsum,tPS,tSF,tFs,tPks
      write(3,6)tsum,tPS,tSF,tFs,tPks
    6 format(' Total time =',f9.3,' Minutes'/' PS =',f9.3,' SF ='f9.3,' 
     * FR =',f9.3,' PPK =',f9.3,' Minutes')
      close(unit=3)
      end
C
      SUBROUTINE SFCLC(NAT)
      parameter mp=2000
      integer*2 jh(mp),jk(mp),jl(mp),IPR(9,3,3)
      integer*2 ih,ik,il,nx(3000,16),narg
      dimension e(mp),p(-mp:mp),TR(9,3)
      dimension cs(-32768:32768),sn(-32768:32768),xxx(1000,3)
      COMMON /SIGNS/ cs,sn,xxx,NLIST
      COMMON /HKL/ NREF,jh,jk,jl,e,p,NSYM,tr,ipr
      do 1 i=1,NAT
      xx=xxx(i,1)
      yy=xxx(i,2)
      zz=xxx(i,3)
      i3=3*i
      i2=i3-1
      i1=i3-2
      do 8 j=1,NSYM
      x = ipr(j,1,1)*xx +ipr(j,1,2)*yy +ipr(j,1,3)*zz +tr(j,1)
      nx(i1,j)=65536.*amod(x,1.0)
      y = ipr(j,2,1)*xx +ipr(j,2,2)*yy +ipr(j,2,3)*zz +tr(j,2)
      nx(i2,j)=65536.*amod(y,1.0)
      z = ipr(j,3,1)*xx +ipr(j,3,2)*yy +ipr(j,3,3)*zz +tr(j,3)
    8 nx(i3,j)=65536.*amod(z,1.0)
    1 continue
      NAT3=3*NAT
      do 10 i=1,NREF
      ih=jh(i)
      ik=jk(i)
      il=jl(i)
      aa=0.
      bb=0.
      do 3 k=1,NSYM
      do 7 j=1,NAT3,3
      narg = ih*nx(j,k) + ik*nx(j+1,k) + il*nx(j+2,k)
      aa = aa + cs(narg)
    7 bb = bb + sn(narg)
    3 continue
      p(i)=atan2(bb,aa)
   10 p(-i)=-p(i)
      RETURN
      end
C
C     3-D SECTIONAL FFT
      SUBROUTINE EDMAP(pks,VAR)
C     RADIX-2 FFT  -  PLUS PEAK LISTING SUBROUTINES
      parameter (NPAR=128,NPA=NPAR-1,NP=64,mp=2000)
      byte nt(24,9),nh(-NP:NP,-NP:NP),HH,HK,HL
      real*8 arg,pi2,px,dp
      integer*2 ih(mp),ik(mp),il(mp),jh,jk,jl,nmy(0:NPA),nmz(0:NPA)
      integer*2 IPR(9,3,3)
      dimension r(-1:NPAR,-1:NPAR,3),t(4,4),tx(24,3),TR(9,3)
      dimension e(mp),ph(-mp:mp)
      complex f(0:NP,-NP:NP,-NP:NP),cx(0:NPA),cy(0:NPA),cz(0:NPA)
      complex p(0:NPA,0:NPA),fx(0:NPA),fhkl,fs
      COMMON /WORK/ fx,cx,cy,cz
      COMMON /PKS/ r,t,peak,MPX,MPAX,MPY,MPZ,IHORZ,IVERT,ISECT,NATOM
      COMMON /HKL/ NREF,ih,ik,il,e,ph,NSYM,tr,ipr
      data mh,mk,ml,mmk,mml/5*0/
      data t/18.,4*-10.,11.,2*4.,-10.,4.,11.,4.,-10.,2*4.,11./
      npeaks=0
      rlast=1.5*var
      peak=1.5*var
      f(0,0,0)=0.1
      if(t(1,1).ne.18) goto 21
      open(unit=1,file='fft0.dat',status='old',form='formatted')
      READ(1,*)MPX,MPAX,MPY,MPZ,IHORZ,IVERT,ISECT,NATOM
      READ(1,*)NSYM
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
   19 if(i.ne.1.or.j.ne.1) t(i,j)=t(i,j)/42.
   21 continue
      do 26 i=1,NREF
      do 28 j=1,nsym
      jh = ih(i)*nt(j,1) + ik(i)*nt(j,2) + il(i)*nt(j,3)
      jk = ih(i)*nt(j,4) + ik(i)*nt(j,5) + il(i)*nt(j,6)
      jl = ih(i)*nt(j,7) + ik(i)*nt(j,8) + il(i)*nt(j,9)
      sh = ih(i)*tx(j,1) + ik(i)*tx(j,2) + il(i)*tx(j,3)
      sh = mod(sh,1.0)
      phi=ph(i) + sh*pi2
      phi=mod(phi,pi2)
      call SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
      if(jh) 32,30,29
   30 if(jk) 32,31,29
   31 if(jl) 32,29,29
   32 jh=-jh
      jk=-jk
      jl=-jl
      phi=-phi
   29 if(t(1,1).ne.18.) goto 28
      mh=max(jh,mh)
      mk=max(jk,mk)
      ml=max(jl,ml)
      mmk=min(jk,mmk)
      mml=min(jl,mml)
      nh(jk,jl)=1
   28 f(jh,jk,jl) = e(i)*cmplx(cos(phi),sin(phi))
   26 continue
      if(t(1,1).eq.18) t(1,1)=t(1,1)/42.
      pks=0.
C     calculate FOURIER
      do 1 ix=-1,MPAX+1
      ixx = ix+1
      ixx = 1+mod(ixx,3)
      do 3 iik = mmk,mk
      do 4 iz=0,MPZ1
    4 fx(iz)=0.
      do 5 iil = mml,ml
      if(nh(iik,iil).eq.0) goto 5
      jl=nmz(iil.and.MPZ1)
      fs=0.
      do 2 iih=0,mh
      fhkl=f(iih,iik,iil)
      if(fhkl.eq.0.0) goto 2
      fs=fs +fhkl*cx((iih*ix).and.MPX1)
    2 continue
C     SUM the (IL+n*MPAR) components into the IL subspace
      fx(jl)=fx(jl)+fs
    5 continue
      call FFT(MPZ,3)
      jk=iik.and.MPY1
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
   10 r(iy,iz,ixx)=fx(iy)
    7 continue
      if(ix.eq.(MPAX+1)) mid=-mid
      pk=secnds(0.0)
      if(ix.ge.1) call pkpick(ix-1,mid,npeaks,rlast)
      pks=pks+secnds(pk)
      mid=ixx
    1 continue
      if(npeaks.lt.NATOM) write(*,13) npeaks
   13 format(' only ',i3,' peaks noted in map ')
      end
C     
      SUBROUTINE FFT(N,IT)
      PARAMETER (NPAR=128,NPA=NPAR-1)
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
      parameter NPAR=128
      real*8 a,b,c,d,e,f,g,ah,ak,al,a1,a2,a3,q1,q2,q3,q4,det
      real*8 t11,t12,t13,t22,t23,t33
      dimension r(-1:NPAR,-1:NPAR,3),t(4,4),r27(27),aname(2)
      dimension px(1000,4),rh(1000),ax(3),nh(27),nk(3)
      dimension cs(-32768:32768),sn(-32768:32768),xxx(1000,3)
      COMMON /PKS/ r,t,peak,MPX,MPAX,MPY,MPZ,IHORZ,IVERT,ISECT,NATOM
      COMMON /SRT/ px,rh,rholast
      COMMON /SIGNS/ cs,sn,xxx,NLIST
      data good,bad,bname/'   G','   B','.pks'/
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
C     are the surrounding points all positive densities
      do 51 i1=1,27
   51 if(nh(i1).eq.1.and.zz.ge.r27(i1)) goto 7
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
      if(dabs(det).lt.0.0001) goto 7
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
      char = good
      goto 8
    7 N=N+1
      if(N.gt.NATOM) N=NATOM+1
      ax(1) = PARX*float(ixm1)
      ax(2) = PARY*float(iy)
      ax(3) = PARZ*float(iz)
      rh(N) = xx
      char = bad
    8 continue
      px(N,1)=ax(IHORZ)
      px(N,2)=ax(IVERT)
      px(N,3)=ax(ISECT)
      px(N,4)=char
      if(N.eq.NATOM) call SORT(N)
      if(N.eq.(NATOM+1).and.rh(N).gt.rholast) call INSERT(NATOM)
      npeaks=N
      rlast=rholast
      peak=0.85*rholast
    6 continue
    5 continue
      if(nstop.eq.1) return
      if(N.lt.NATOM.and.N.gt.0) call SORT(N)
      if(N.gt.NATOM) N = NATOM
      do 11 i=1,N
      xxx(i,1)=px(i,1)
      xxx(i,2)=px(i,2)
      xxx(i,3)=px(i,3)
      IF(NLIST.EQ.2)write(99,98)px(i,1),px(i,2),px(i,3),rh(i),i,px(i,4)
   98 format(3f9.5,f10.2,i5,5x,a4)
   11 CONTINUE
      RETURN
      end
C
      subroutine SORT(N)
C     sort the peaks in descending order of magnitude
      integer*2 num(1000)
      dimension px(1000,4),rh(1000),q(1000,4)
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
      return
      end
C
      subroutine INSERT(NATOM)
C     place newest peak within the sorted stack
      dimension px(1000,4),rh(1000)
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
      return
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
      return
      end
