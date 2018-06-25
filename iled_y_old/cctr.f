C     CORRELATION COEFFICIENT TRANSLATION FUNCTION
C     NX,NX,NY ARE THE RANGES OF HKL
      PARAMETER (NX=100,NY=150)
      REAL*8 PI,TR,DOUBL,EO2,EO4,AVE2,AVE4,VOBS,SOBS
      INTEGER*2 IH,IK,IL,NS,NE,NF(99),IPR(9,3,3),JH(8),JK(8),JL(8),NP(8)
      DIMENSION F(8),X(90),Y(90),Z(90),NFJ(90),EC(9999)
      DIMENSION AN(2),TR(12,3),CS(-32768:32768),SN(-32768:32768)        *73
      DIMENSION AM(0:NX,-NX:NX,-NY:NY),BM(0:NX,-NX:NX,-NY:NY)
      DIMENSION CM(0:NX,-NX:NX,-NY:NY),DM(0:NX,-NX:NX,-NY:NY)
      DIMENSION EM(0:NX,-NX:NX,-NY:NY),FM(0:NX,-NX:NX,-NY:NY)
      DATA MH,MK,ML,NK,NL,SUMF2,EO2,EO4/5*0,3*0./
C
      PI=3.141592654
      DOUBL=PI/32768.
      DO 5 I=-32768,32768
      ARG=I*DOUBL
      CS(I)=COS(ARG)
    5 SN(I)=SIN(ARG)
      OPEN(1,FILE='fft.dat',STATUS='OLD',FORM='FORMATTED')
      READ(1,*)NSYM
      READ(1,*)NSYM
      DO 1 I=1,NSYM
C     READ IN THE SYMMETRY MATRICES IN TRANSPOSED ORDER
      READ(1,2) (TR(I,J),(IPR(I,K,J),K=1,3),J=1,3)
    2 FORMAT(3(F10.5,3I5))
      DO 3 J=1,3
      SCALE=24.*TR(I,J)
      IX=NINT(SCALE)
    3 TR(I,J)=IX/24.
    1 CONTINUE
      CLOSE(1)
      WRITE(*,4)
    4 FORMAT(' DX, DY, DZ :',$)
      read(*,*)dx,dy,dz
   17 FORMAT(A4)
      OPEN(1,FILE='iled.MOL',STATUS='OLD',FORM='FORMATTED')
      DO 6 I=1,90
      READ(1,*,END=22)X(I),Y(I),Z(I),NFJ(I)
      sumf2=sumf2+nfj(i)**2
      x(i)=x(i)+dx
      y(i)=y(i)+dy
    6 z(i)=z(i)+dz
   22 NAT=I-1
      SCALE=1./SQRT(FLOAT(NSYM)*sumf2)
C
      CLOSE(1)
      OPEN(1,FILE='iled.E',STATUS='OLD',FORM='FORMATTED')
      DO 19 I=1,999999
      READ(1,*,END=18)IH,IK,IL,NS,NE
      ix=ih*(IK+IL) + ik*il
      E=0.001*NE
      if(ix.eq.0) e=1.4142*e
      a=0
      b=0
      do 102 j=1,nat
      xj=x(j)-dx
      yj=y(j)-dy
      zj=z(j)-dz
      do 101 k=1,NSYM
      xx = xj*ipr(k,1,1) + yj*ipr(k,2,1) + zj*ipr(k,3,1) + tr(k,1)
      yy = xj*ipr(k,1,2) + yj*ipr(k,2,2) + zj*ipr(k,3,2) + tr(k,2)
      zz = xj*ipr(k,1,3) + yj*ipr(k,2,3) + zj*ipr(k,3,3) + tr(k,3)
      NS=65536*(ih*xx +ik*yy +il*zz)
      a=a+nfj(j)*cs(NS)
  101 b=b+nfj(j)*sn(NS)
  102 continue
      ec(i)=scale*sqrt(a*a+b*b)
C
C      e=ec(i)
      EO2=EO2+E*E
   19 EO4=EO4+E**4
   18 NREF=I-1
      AVE2=EO2/NREF
      AVE4=EO4/NREF
      VOBS=AVE4-AVE2*AVE2
      write(*,*)ave2,ave4,vobs
C
      SOBS=SQRT(VOBS)
      REWIND(1)
      top=0.
      bot=0.
      T0 = SECNDS(0.0)
      DO 7 I=1,NREF
      READ(1,*,END=7)IH,IK,IL,NS,NE
      ix=ih*(IK+IL) + ik*il
      E=0.001*NE
      if(ix.eq.0) e=1.4142*e
C
C      e=ec(i)
      EO2=E*E
      DO 8 J=1,NSYM
      A=0
      B=0
      I1=IH*IPR(J,1,1) +IK*IPR(J,1,2) +IL*IPR(J,1,3)
      I2=IH*IPR(J,2,1) +IK*IPR(J,2,2) +IL*IPR(J,2,3)
      I3=IH*IPR(J,3,1) +IK*IPR(J,3,2) +IL*IPR(J,3,3)
      JH(J)=I1
      JK(J)=I2
      JL(J)=I3
      TS=TR(J,1)*IH +TR(J,2)*IK +TR(J,3)*IL
      DO 9 K=1,NAT
      NS = 65536*(I1*X(K) +I2*Y(K) +I3*Z(K)+TS)
      A=A+NFJ(K)*CS(NS)
    9 B=B+NFJ(K)*SN(NS)
      F(J)=SCALE*SQRT(A*A+B*B)
      NP(J)=ATAN2(B,A)/DOUBL
    8 CONTINUE
      top=top+abs(ec(i)-e)
      bot=bot+abs(e)
C
      AVE= (EO2-AVE2)
      DO 20 J1=1,NSYM
      DO 21 J2=1,NSYM
      IH=(JH(J1)-JH(J2))/2
      IK=(JK(J1)-JK(J2))/2
      IL=(JL(J1)-JL(J2))/2
      NS=NP(J1)-NP(J2)
      IF(IH) 13,14,15
   14 IF(IK) 13,16,15
   16 IF(IL) 13,15,15
   13 IH=-IH
      IK=-IK
      IL=-IL
      NS=-NS
   15 FF=F(J1)*F(J2)/NREF
C     SUM ECAL**2
      FG=SOBS*FF
      AM(IH,IK,IL)=AM(IH,IK,IL) +FG*CS(NS)
      BM(IH,IK,IL)=BM(IH,IK,IL) +FG*SN(NS)
      FG=FF*AVE
C     SUM (EO**2 -<EO**2>)*ECAL**2
      CM(IH,IK,IL)=CM(IH,IK,IL) +FG*CS(NS)
      DM(IH,IK,IL)=DM(IH,IK,IL) +FG*SN(NS)
   21 CONTINUE
   20 CONTINUE
C
      DO 120 J1=1,NSYM
      DO 121 J2=1,NSYM
      IH0=JH(J1)-JH(J2)
      IK0=JK(J1)-JK(J2)
      IL0=JL(J1)-JL(J2)
      MS=NP(J1)-NP(J2)
      FF=VOBS*F(J1)*F(J2)/NREF
      DO 122 J3=1,NSYM
      DO 123 J4=1,NSYM
      IH=(IH0+JH(J3)-JH(J4))/2
      IK=(IK0+JK(J3)-JK(J4))/2
      IL=(IL0+JL(J3)-JL(J4))/2
      NS=MS+NP(J3)-NP(J4)
      IF(IH) 113,114,115
  114 IF(IK) 113,116,115
  116 IF(IL) 113,115,115
  113 IH=-IH
      IK=-IK
      IL=-IL
      NS=-NS
  115 MH=MAX(MH,IH)
      MK=MAX(MK,IK)
      ML=MAX(ML,IL)
      NK=MIN(NK,IK)
      NL=MIN(NL,IL)
      FG=FF*F(J3)*F(J4)
      EM(IH,IK,IL)=EM(IH,IK,IL) +FG*CS(NS)
      FM(IH,IK,IL)=FM(IH,IK,IL) +FG*SN(NS)
  123 CONTINUE
  122 CONTINUE
  121 CONTINUE
  120 CONTINUE
    7 CONTINUE
      r = top/bot
      write(*,*)top,bot,r
      CLOSE(1)
C
C     WRITE FOURIER COEFFICIENT FILES FOUR-1,2,3 FOR MAPS-1,2,3
C     CC = MAP2/VOBS*SQRT(MAP3-MAP1**2)
C
      OPEN(1,FILE='FOUR1',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(2,FILE='FOUR2',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(3,FILE='FOUR3',STATUS='UNKNOWN',FORM='UNFORMATTED')
      DO 10 IL=NL,ML
      DO 11 IK=NK,MK
      DO 12 IH= 0,MH
      A=AM(IH,IK,IL)
      B=BM(IH,IK,IL)
      IF(A.EQ.0.0.AND.B.EQ.0.0) GOTO 112
      G=SQRT(A*A+B*B)
      NS=1000.*ATAN2(B,A)
      WRITE(1)IH,IK,IL,G,NS
  112 A=CM(IH,IK,IL)
      B=DM(IH,IK,IL)
      IF(A.EQ.0.0.AND.B.EQ.0.0) GOTO 111
      G=SQRT(A*A+B*B)
      NS=1000.*ATAN2(B,A)
      WRITE(2)IH,IK,IL,G,NS
  111 A=EM(IH,IK,IL)
      B=FM(IH,IK,IL)
      IF(A.EQ.0.0.AND.B.EQ.0.0) GOTO 12
      G=SQRT(A*A+B*B)
      NS=1000.*ATAN2(B,A)
      WRITE(3)IH,IK,IL,G,NS
   12 CONTINUE
   11 CONTINUE
   10 CONTINUE
      CALL EDMAP
      T1=SECNDS(T0)
      write(*,*)MH,MK,ML,NK,NL
      write(*,*)NAT,NREF,r,T1
      END
C      
C     RADIX-2 FFT  -  PLUS PEAK LISTING SUBROUTINES
      SUBROUTINE EDMAP
      parameter (NPAR=256,NPA=NPAR-1,NP=100)
      byte nt(24,9),nh(-NP:NP,-NP:NP,3),HH,HK,HL
      real*8 arg,pi2,px
      integer*2 ih,ik,il,jh,jk,jl,ns,ne,nph,nm(0:NPA)
      integer*2 mh(3),mk(3),ml(3),mmk(3),mml(3)
      dimension r(-1:NPAR,-1:NPAR,3),t(4,4),tx(24,3)
      dimension rr(-1:NPAR,-1:NPAR,3)
      complex f(0:NP,-NP:NP,-NP:NP,3),cx(0:NPA),fx(0:NPA),fhkl,cxl,fs
      complex p(0:NPA,0:NPA)
      COMMON /WORK/ fx,cx,n1,NX
      COMMON /PKS/ rr,t,peak,MPAR,MPA,MPAX,IHORZ,IVERT,ISECT,NATOM
      data npeaks,num,mh,mk,ml,mmk,mml,rlast,peak/17*0,2*0.0/
      data t/18.,4*-10.,11.,2*4.,-10.,4.,11.,4.,-10.,2*4.,11./
      do 19 i=1,4
      do 19 j=1,4
   19 t(i,j)=t(i,j)/42.
      px=-1.0
      OPEN(unit=4,file='fft.dat',status='old',form='formatted')
      READ(4,*)MPAR,MPAX,IHORZ,IVERT,ISECT,NATOM,LHALF
      close(unit=4)
      MPA=MPAR-1
      MP2=MPAR+MPA
      pi2=2.0*dacos(px)
      pie2=pi2
      px=pi2/MPAR
C     initialize exp(-twopi*x/MPAR), remember negative sign
      do 20 i=0,MPA
      arg=-i*px
   20 cx(i)=dcmplx(dcos(arg),dsin(arg))
C     zero data arrays
      do 21 iset=1,3
      do 21 il=-NP,NP
      do 21 ik=-NP,NP
      nh(ik,il,iset) = 0
      do 21 ih=  0,NP
   21 f(ih,ik,il,iset)=0.
C     set up FFT arrays for bit reversal on first pass
      n1=MPAR-1
      nm(0)=0
      NX = alog(float(MPAR))/alog(2.0) - 0.9
      mm=1
      do 23 i=1,NX+1
      do 24 j=0,mm-1
      nm(j)=2*nm(j)
   24 nm(j+mm) = nm(j)+1
      mm=mm+mm
   23 continue
C     read reflection data
      rewind(1)
      rewind(2)
      rewind(3)
      do 100 iset=1,3
      do 26 i=1,999999
      read(iset,end=27)jh,jk,jl,ff,nph
      ph=0.001*nph
      call SWAP(jh,jk,jl,IHORZ,IVERT,ISECT)
      if(jh) 32,30,29
   30 if(jk) 32,31,29
   31 if(jl) 32,29,29
   32 jh=-jh
      jk=-jk
      jl=-jl
      ph=-ph
   29 mh(iset)=max(jh,mh(iset))
      mk(iset)=max(jk,mk(iset))
      ml(iset)=max(jl,ml(iset))
      mmk(iset)=min(jk,mmk(iset))
      mml(iset)=min(jl,mml(iset))
      nh(jk,jl,iset)=1
      f(jh,jk,jl,iset) = ff*cmplx(cos(ph),sin(ph))
   26 continue
   27 close(unit=iset)
  100 continue
      write(*,765)mh,mk,ml,mmk,mml
  765 format(5i8)    
C     calculate FOURIER
      do 1 ix=-1,MPAX+1
      ixx = ix+1
      ixx = 1+mod(ixx,3)
C      
      do 200 iset=1,3
      do 3 ik = mmk(iset),mk(iset)
      do 4 iz=0,MPA
    4 fx(iz)=0.
      do 5 il = mml(iset),ml(iset)
      if(nh(ik,il,iset).eq.0) goto 5
      jl=iand(il,MPA)
      jl=nm(jl)
      fs=0.
      do 2 ih=0,mh(iset)
      fhkl=f(ih,ik,il,iset)
      if(fhkl.eq.0.0) goto 2
      ihix=ih*ix
      ihix=iand(ihix,MPA)
      cxl = cx(ihix)
      fs=fs +fhkl*cxl
    2 continue
      fx(jl)=fs
    5 continue
      call FFT(MPAR)
      jk=iand(ik,MPA)
      do 6 iz=0,MPA
    6 p(jk,iz) = fx(iz)
    3 continue
      do 7 iz=0,MPA
      do 8 ik=0,MPA
    8 fx(ik)=0.
      do 9 ik = mmk(iset),mk(iset)
      jk=iand(ik,MPA)
      kj=nm(jk)
    9 fx(kj)=p(jk,iz)
      call FFT(MPAR)
      do 10 iy=0,MPA
   10 r(iy,iz,iset)=fx(iy)
    7 continue
  200 continue
      do 300 iz=0,MPA
      do 300 iy=0,MPA
      cc = r(iy,iz,2)/SQRT(r(iy,iz,3)-r(iy,iz,1)**2)
  300 rr(iy,iz,ixx)=cc
      if(ix.eq.(MPAX+1)) mid=-mid
      if(ix.ge.1) call pkpick(ix-1,mid,npeaks,rlast,LHALF)
      mid=ixx
    1 continue
      if(npeaks.lt.NATOM) write(*,41) npeaks
   41 format(' only ',i3,' peaks noted in map ')
      RETURN
      END
C     
      SUBROUTINE FFT(N)
      parameter (NPAR=256,NPA=NPAR-1)
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
      a=f(i)
      ii=i+i0
      b=f(ii)
      f(i)=a+b
      f(ii)=a-b
      i4=0
      do 4 k=i+1,i+i2
      i4=i4+m2
      a = f(k)
      kk=k+i0
      b = f(kk)*cx(i4)
      f(k)=a+b
    4 f(kk)=a-b
    3 continue
      m2=m2/2
    2 continue
      return
      end
C
      SUBROUTINE PKPICK(ixm1,MID,npeaks,rlast,LHALF)
      parameter (NPAR=256,NPA=NPAR-1,NP=100)
      dimension r(-1:NPAR,-1:NPAR,3),t(4,4),r27(27)
      dimension px(1000,4),rh(1000),ax(3)
      INTEGER*2 nh(27),nk(3),NGOOD(1000)
      COMMON /PKS/ r,t,peak,MPAR,MPA,MPAX,IHORZ,IVERT,ISECT,NATOM
      COMMON /SRT/ px,rh,rholast
      data nh/0,1,0,1,1,1,0,1,0,4*1,0,4*1,0,1,0,1,1,1,0,1,0/
      nstop=mid/iabs(mid)
      mid=iabs(mid)
      ntop=mid+1
      if(ntop.eq.4) ntop=1
      nbot=6-mid-ntop
C     fill out edges of the map
      do 1 ix=1,3
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
      if(npeaks.gt.10) goto 12
      pick=0.
      MPA2=MPA
      IF(LHALF.NE.0) MPA2=MPA/2
      do 4 iz=0,MPA2
      do 4 iy=0,MPA2
      xx = r(iy,iz,mid)
    4 if(xx.gt.pick) pick=xx
      pick = 0.15*pick
      if(pick.gt.peak) peak=pick
   12 N = npeaks
      rholast=rlast
      if(N.eq.0) OPEN(1,file='PEAKS',status='unknown',form='formatted')
      PAR = 1.0/FLOAT(MPAR)
C     find peaks surrounded by 19 positive points
      ix=mid
      ixm=nbot
      ixp=ntop
      nk(1)=nbot
      nk(2)=mid
      nk(3)=ntop
      do 5 iz=0,MPA2
      izm=iz-1
      izp=iz+1
      do 6 iy=0,MPA2
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
      if(abs(det).lt.0.0001) goto 7
      t11 = (4*f*g-ah*ah)/det
      t12 = (ah*ak-2*g*al)/det
      t13 = (ah*al-2*f*ak)/det
      t22 = (4*e*g-ak*ak)/det
      t23 = (ak*al-2*e*ah)/det
      t33 = (4*e*f-al*al)/det
      x = t11*b + t12*c + t13*d
      y = t12*b + t22*c + t23*d
      z = t13*b + t23*c + t33*d
      ax(1) = PAR*(float(ixm1)+x)
      if(ax(1).lt.-0.0001.or.ax(1).gt.1.0001) goto 6
      N=N+1
      if(N.gt.NATOM) N=NATOM+1
      rh(N) = exp(a+b*x+c*y+d*z+e*x*x+f*y*y+g*z*z+ah*z*y+ak*x*z+al*x*y)
      ax(2) = PAR*(float(iy)+y)
      ax(3) = PAR*(float(iz)+z)
      MGOOD = 0
      goto 8
    7 N=N+1
      if(N.gt.NATOM) N=NATOM+1
      ax(1) = PAR*float(ixm1)
      ax(2) = PAR*float(iy)
      ax(3) = PAR*float(iz)
      rh(N) = xx
      MGOOD = -1
    8 px(N,1)=ax(IHORZ)
      px(N,2)=ax(IVERT)
      px(N,3)=ax(ISECT)
      NGOOD(N)=MGOOD
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
      do 9 i=1,N
    9 write(1,10) (px(i,j),j=1,3),rh(i),NGOOD(I),i
   10 format(4f8.4,2i5)
      return
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
