C     PROGRAM = P10 formula
      PARAMETER (IE=60000)
      CHARACTER an*12,bn*12
      byte JH(IE,24),JK(IE,24),JL(IE,24),NR(IE),IR(24,9),NC(0:IE)
      dimension e(0:IE),p(-IE:IE),fj(8)
      integer*4 LOC(99999),H1,H2,H3,NA(8)
      COMMON /HKL/ LAUE,MH,MK,ML,NK,NL,IHR,IKR,LOC
      COMMON /MADD/W,WT,LLX,LLL,LSX,LSL,LIM1,LIM2,NOP,NC
      E(0)=0.
      mum=0
      nss=0
      mss=0
      open(1,file='sym.dat',form='formatted')
      read(1,*)NSYM,LAUE,ICENT,LAT
      read(1,18)sct
      read(1,18)atms
      DO 30 I = 1,16
      Z2 = Z2 + (2-ICENT)*NEQ*ATMS(I)*SCT(I)**2
   30 Z3 = Z3 + (2-ICENT)*NEQ*ATMS(I)*SCT(I)**3
      SCALE = 2.0*Z3/(Z2**1.5)
      do 1 i=1,NSYM
      READ(1,2)T1,IH,IK,IL,T2,KH,KK,KL,T3,LH,LK,LL
    2 FORMAT(3(F10.5,3I5),i5)
      IX=ABS(IH)+ABS(KK)+ABS(LL)
      IY=ABS(IK)+ABS(IL)+ABS(KH)+ABS(KL)+ABS(LH)+ABS(LK)
      IF(IY.EQ.0.AND.IX.NE.3) GOTO 3
      IDET = IH*KK*LL+IK*KL*LH+KH*LK*IL-LH*KK*IL-KH*IK*LL-IH*LK*KL
      IF(ABS(IDET).NE.1) GOTO 3
      IR(I,1) = (KK*LL-LK*KL)/IDET
      IR(I,2) = -(KH*LL-LH*KL)/IDET
      IR(I,3) = (KH*LK-KK*LH)/IDET
      IR(I,4) = -(IK*LL-IL*LK)/IDET
      IR(I,5) = (IH*LL-LH*IL)/IDET
      IR(I,6) = -(IH*LK-LH*IK)/IDET
      IR(I,7) = (IK*KL-KK*IL)/IDET
      IR(I,8) = -(IH*KL-IL*KH)/IDET
    1 IR(I,9) = (IH*KK-IK*KH)/IDET
      CLOSE(1)
      write(*,21)
   21 format(' E-data file : name.E : ',$)
      read(*,22)an
   22 format(a)
      write(*,27)
   27 format(' triple file : name.T : ',$)
      read(*,22)bn
      write(*,26)
   26 format(' Ebig,NTRIPS,NOP (0 or 1 for wts): ',$)
      read(*,*)EBIG,NTRIP,NOP
      open(1,file=an,form='formatted')
      do 19 I=1,IE
      read(1,*,end=6)ih,ik,il,ns,ne
      en=0.001*ne
      E(I)=en*en-1.0
   19 if(en.ge.EBIG) LIM1=I
    6 NREF=I-1
      rewind(1)
      LIM2=NREF-1.25*LIM1
      do 5 I=1,NREF
      read(1,*)ih,ik,il,ns,ne
      NC(I)=2
      if(I.le.LIM1) NC(I)=1
      if(I.EQ.LIM2) ESMALL=0.001*Ne
      if(I.ge.LIM2) NC(I)=4
      mh=max(ih,mh)
      mk=max(ik,mk)
      ml=max(il,ml)
      nk=min(ik,nk)
    5 nl=min(il,nl)
      LIM3=NREF-LIM2
      REWIND(1)
      V0=FLOAT(LIM1)/FLOAT(NREF)
      write(*,*)mh,mk,ml,nk,nl
      write(*,23)NREF,LIM1,EBIG,LIM3,ESMALL,V0
   23 format(i5,' Es,',i5,' Es >',f5.2,',',i5,' Es <',f5.2,', V =',f6.3)  
      IHR=MH+1
      IKR=MK-NK+1
      nref2=nref/2
      do 7 I=1,NREF
      read(1,*)ih,ik,il,ns,ne,a,b,c,np
      JH(I,1)=ih
      JK(I,1)=ik
      JL(I,1)=il
      p(i)=0.001*np
      p(-i)=-p(i)
      NCT=1
      DO 9 J=2,NSYM
      KH=IH*IR(J,1) +IK*IR(J,2) + IL*IR(J,3)
      KK=IH*IR(J,4) +IK*IR(J,5) + IL*IR(J,6)
      KL=IH*IR(J,7) +IK*IR(J,8) + IL*IR(J,9)
      NX=0
      DO 10 K=1,J-1
      IX=IABS(JH(I,K)-KH) + IABS(JK(I,K)-KK) + IABS(JL(I,K)-KL)
      JX=IABS(JH(I,K)+KH) + IABS(JK(I,K)+KK) + IABS(JL(I,K)+KL)
   10 IF(IX.EQ.0.OR.JX.EQ.0) NX=1
      IF(NX.EQ.1) GOTO 9
      NCT=NCT+1
      JH(I,NCT)=KH
      JK(I,NCT)=KK
      JL(I,NCT)=KL
    9 CONTINUE
      NR(I)=NCT
    8 NLOC = 1 + IH +IHR*(IK-NK +IKR*(IL-NL))
    7 LOC(NLOC)=I
      close(1)
      open(1,file=bn,form='formatted')
C      
      do 11 i=1,NTRIP
      read(1,*,end=12)aa,ires,it,n1,n2,n3
      if(ires.eq.7) goto 11
      arg=cos(p(n1)+p(n2)+p(n3))
      if(it.eq.6) arg=-arg
      H1=jh(n1,1)
      H2=jk(n1,1)
      H3=jl(n1,1)
      eh=e(n1)
      n4=abs(n2)
      ek=e(n4)
      j2=n4/n2
      n5=abs(n3)
      el=e(n5)
      j3=n5/n3
      n6=nc(N1)
      n7=nc(n4)
      n8=nc(n5)
      nc(N1)=-9
      nc(N4)=-9
      nc(n5)=-9
      do 13 j=1,NR(n4)
      K1=j2*jh(n4,j)
      K2=j2*jk(n4,j)
      K3=j2*jl(n4,j)
      do 14 k=1,NR(n5)
      L1=j3*jh(n5,k)
      L2=j3*jk(n5,k)
      L3=j3*jl(n5,k)
      ix=abs(h1+k1+l1)+abs(h2+k2+l2)+abs(h3+k3+l3)
   14 if(ix.eq.0) goto 15
   13 continue
      ix=1
   15 if(ix.ne.0) write(*,*)99   
C
      AS=0.
      BS=0.
C      
      do 16 J=1,LIM1
      em=e(j)
      if(j.eq.n1.or.j.eq.n4.or.j.eq.n5) goto 16
      JNC=NC(J)
      NC(J)=-9
      do 17 k=1,NR(J)
      M1=jh(j,k)
      M2=jk(j,k)
      M3=jl(j,k)
C     H+M, IS1      
      IH=H1+M1
      IK=H2+M2
      IL=H3+M3
      CALL STD(IH,IK,IL,I1)
C     K+M, IS2   
      IH=K1+M1
      IK=K2+M2
      IL=K3+M3
      CALL STD(IH,IK,IL,I2)
C     L+M, IS3
      IH=L1+M1
      IK=L2+M2
      IL=L3+M3
      CALL STD(IH,IK,IL,I3)
C     H-M, IS4
      IH=H1-M1
      IK=H2-M2
      IL=H3-M3
      CALL STD(IH,IK,IL,I4)
C     K-M, IS5
      IH=K1-M1
      IK=K2-M2
      IL=K3-M3
      CALLSTD(IH,IK,IL,I5)
C     L-M, IS6
      IH=L1-M1
      IK=L2-M2
      IL=L3-M3
      CALL STD(IH,IK,IL,I6)
      IQ=1
      if(i1.le.lim1.and.i1.ne.0) iq=iq+1
      if(i2.le.lim1.and.i2.ne.0) iq=iq+1
      if(i3.le.lim1.and.i3.ne.0) iq=iq+1
      if(i4.le.lim1.and.i4.ne.0) iq=iq+1
      if(i5.le.lim1.and.i5.ne.0) iq=iq+1
      if(i6.le.lim1.and.i6.ne.0) iq=iq+1
      e1=e(i1)
      e2=e(i2)
      e3=e(i3)
      e4=e(i4)
      e5=e(i5)
      e6=e(i6)
      AM = Em*(e1*(e5+e6)+e2*(e4+e6)+e3*(e4+e5))
      BM =      eh*(em*(e1+e4)+e2*e6+e5*e3)
      bm = bm + ek*(em*(e2+e5)+e1*e6+e4*e3)
      bm = bm + el*(em*(e3+e6)+e1*e5+e4*e2)
      as=as+am
      bs=bs+bm
      if(iq.gt.4) write(42,29)IQ,J,I1,I2,I3,I4,I5,I6
   29 format(8i6)       
   17 CONTINUE
      NC(J)=JNC
   16 CONTINUE
   20 CONTINUE
      NN=84*4
      v1=0.
      v2=0.
      if(ires.eq.1) as=as/2.
      if(ires.eq.1) bs=bs/2.
      as = as/NN
      ehkl = eh*ek*el
      cs=0.5*(ehkl+bs)/nn    
      if(cs.lt.0.0) cs=0.
      q=as/(1+cs)
      G = AA*(1.+q)
      write(2,28)aa,ires,it,n1,n2,n3,v1,v2,g,arg
   28 format(f8.3,3i5,2i6,3f8.4,f6.2,2f5.1)
   24 format(f6.2,3(i5,2i4),7i5,3f8.3,f6.2,2f5.1)
   11 continue
   12 continue
      nss = NINT(NSS/FLOAT(NTRIP))
      mss = NINT(MSS/FLOAT(NTRIP))
      write(*,*)mum,ntrip,nss,mss
      goto 25
    3 write(*,4)
    4 format(' SYMMETERY CARD ERROR')
   25 end
      SUBROUTINE STD(IH,IK,IL,ISER)
C     ******************************************************************
C     SUBROUTINE TO DETERMINE STANDARD PARENT FORM
C     ******************************************************************
      INTEGER*4 LOC(99999)
      COMMON /HKL/ LAU,MH,MK,ML,NK,NL,IHR,IKR,LOC
      LAUE=IABS(LAU)
      IF(LAUE.GT.3) GOTO 1
      IF(LAUE-2) 2,3,4
C     TRICLINIC SYMMETRY
    2 IF(IH) 5,6,7
    6 IF(IK) 5,8,7
    8 IF(IL) 5,7,7
    5 IH = -IH
      IK = -IK
      IL = -IL
    7 ISER=0
      IF(IH.GT.MH.OR.IK.GT.MK.OR.IL.GT.ML) RETURN
      IF(IK.LT.NK.OR.IL.LT.NL) RETURN
      NLOC = 1 + IH +IHR*(IK-NK +IKR*(IL-NL))
      ISER=LOC(NLOC)
      RETURN
C     MONOCLINIC SYMMETRY -  B-AXIS UNIQUE
    3 IK = IABS(IK)
      IF(IH.GT.0) GOTO 7
      IH = -IH
      IL = -IL
      IF(IH.EQ.0) IL=IABS(IL)
      GOTO 7
C     ORTHORHOMBIC SYMMETRY
    4 IH = IABS(IH)
      IK = IABS(IK)
      IL = IABS(IL)
      GOTO 7
    1 IF(LAUE-5) 9,9,10
C     4/M SYMMETRY
    9 IL = IABS(IL)
      IF(IH*IK) 12,13,14
   12 IX = IABS(IH)
      IH = IABS(IK)
      IK = IX
      GOTO 11
   13 IK = IABS(IH) + IABS(IK)
      IH = 0
      GOTO 11
   14 IH = IABS(IH)
      IK = IABS(IK)
   11 IF(LAUE.EQ.4) GOTO 7
      IX = MIN(IH,IK)
      IK = MAX(IH,IK)
      IH = IX
      GOTO 7
   10 IF(LAUE-10) 15,16,16
C     BAR-3 SYMMETRY
   15 IF(IH*IK) 17,18,19
   17 IX = IH
      IH = -IK
      IK = IX+IK
      IL = -IL
      GOTO 15
   18 IF(IH.EQ.0) GOTO 20
      IF(IH.GT.0) IL=-IL
      IK=IABS(IH)
      IH=0
      GOTO 21
   20 IF(IK.LT.0) IL=-IL
      IK=IABS(IK)
      GOTO 21
   19 IF(IH.GT.0) GOTO 21
      IH=-IH
      IK=-IK
      IL=-IL
   21 IF(LAUE.EQ.6) GOTO 7
      IF(LAUE.GE.8) IL=IABS(IL)
      IF(LAUE.EQ.8) GOTO 7
      IF(IH.LT.IK) GOTO 7
      IX=IH
      IH=IK
      IK=IX
      IF(LAU.EQ.-7) IL=-IL
      IF(LAU.EQ.-7.AND.IH.EQ.IK) IL=IABS(IL)
      GOTO 7
C     CUBIC SYMMETRY
   16 IH=IABS(IH)
      IK=IABS(IK)
      IL=IABS(IL)
      IX=MIN(IH,IK)
      MINH=MIN(IX,IL)
   23 IX=IH
      IH=IK
      IK=IL
      IL=IX
      IF(IH.NE.MINH) GOTO 23
      IF(IH.EQ.MINH.AND.IK.NE.MINH.AND.IL.EQ.MINH) GOTO 23
      IF(LAUE.EQ.10) GOTO 7
      IF(IK.LT.IL) GOTO 7
      IX=IK
      IK=IL
      IL=IX
      GOTO 7
      END
