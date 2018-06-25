      PARAMETER (IE=40000)
      CHARACTER an*6,bn*6
      byte JH(IE,24),JK(IE,24),JL(IE,24),NR(IE),IR(24,9),NN(0:IE)
      integer*4 LOC(99999),H1,H2,H3,SL,SX
      dimension p(-IE:IE)
      COMMON /HKL/ LAUE,MH,MK,ML,NK,NL,IHR,IKR,LOC
      NN(0)=-9
      mum=0
      nss=0
      mss=0
      open(1,file='sym.dat',form='formatted')
      read(1,*)NSYM,LAUE
      read(1,*)a
      read(1,*)n
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
   26 format(' Ebig,NTRIPS,delv : ',$)
      read(*,*)EBIG,NTRIP,dvv
      open(1,file=an,form='formatted')
C     READ THE SORTED E LIST      
      do 29 I=1,IE
      read(1,*,end=6)ih,ik,il,ns,ne
      mh=max(ih,mh)
      mk=max(ik,mk)
      ml=max(il,ml)
      nk=min(ik,nk)
      nl=min(il,nl)
      en=0.001*ne
      NN(I)=2
      if(en.ge.EBIG) NN(I)=1
      if(en.ge.EBIG) LIM1=I
   29 CONTINUE
    6 NREF=I-1
      rewind(1)
      LIM2=NREF-1.25*LIM1
      LIM3=NREF-LIM2
      v0=FLOAT(LIM1)/float(NREF)
      write(*,*)mh,mk,ml,nk,nl
      write(*,23)NREF,LIM1,EBIG,LIM3,ESMALL,v0
   23 format(i5,'=NREF,',i5,' E >',f6.2,',',i5,' E <',f6.2,', V0=',f6.3)
      REWIND(1)
      IHR=MH+1
      IKR=MK-NK+1
C     RE-READ THE SORTED E LIST     
      do 7 I=1,NREF
      read(1,*)ih,ik,il,ns,ne,a,b,c,np
      if(I.ge.LIM2) nn(i)=4
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
C     READ THE TRIPLES LIST      
      do 11 i=1,NTRIP
      read(1,*,end=12)a,ires,it,n1,n2,n3
      arg=cos(p(n1)+p(n2)+p(n3))
      if(it.eq.6) arg=-arg
      H1=jh(n1,1)
      H2=jk(n1,1)
      H3=jl(n1,1)
      n22=abs(n2)
      j2=n22/n2
      n33=abs(n3)
      j3=n33/n3
C     set marker for TRIPLES MAIN TERMS N1,N2,N3
      m11= nn(n1)
      m22=nn(n22)
      m33=nn(n33)  
      NN(n1) =-9
      NN(n22)=-9
      NN(n33)=-9
      do 13 j=1,NR(N22)
      K1=j2*jh(n22,j)
      K2=j2*jk(n22,j)
      K3=j2*jl(n22,j)
      do 14 k=1,NR(N33)
      L1=j3*jh(n33,k)
      L2=j3*jk(n33,k)
      L3=j3*jl(n33,k)
      ix=abs(h1+k1+l1)+abs(h2+k2+l2)+abs(h3+k3+l3)
   14 if(ix.eq.0) goto 15
   13 continue
      ix=1
   15 if(ix.ne.0) write(*,5)99,h1,h2,h3,k1,k2,k3,l1,l2,l3
    5 format(i5,3(i7,2i4))
      LL=0
      LX=0
      SL=0
      SX=0
      do 16 J=1,LIM1
      if(NN(j).eq.-9) goto 16
      do 17 k=1,NR(J)
      do 18 l=-1,1,2
      M1=l*jh(j,k)
      M2=l*jk(j,k)
      M3=l*jl(j,k)
C     M,H-M,K+M, M,ISER1,ISER2      
      IH=H1-M1
      IK=H2-M2
      IL=H3-M3
      CALL STD(IH,IK,IL,ISER1)
      IH=K1+M1
      IK=K2+M2
      IL=K3+M3
      CALL STD(IH,IK,IL,ISER2)
      IX=NN(ISER1)+NN(ISER2)
      IF(IX.LT.0.or.ix.eq.4) goto 30
      IF(IX.EQ.2) LL=LL+1
      IF(IX.LE.5) LX=LX+1
      IF(IX.EQ.5) SL=SL+1
      IF(IX.GE.5) SX=SX+1
C     H-M,M,L+M, ISER1,M,ISER3                  
   30 IH=L1+M1
      IK=L2+M2
      IL=L3+M3 
      CALL STD(IH,IK,IL,ISER3)
      IY=NN(ISER1)+NN(ISER3)
      IF(IY.LE.0.or.iy.eq.4) GOTO 31
      IF(IY.EQ.2) LL=LL+1
      IF(IY.LE.5) LX=LX+1
      IF(IY.EQ.5) SL=SL+1
      IF(IY.GE.5) SX=SX+1
C     K+M,M-l,M, ISER2,ISER4,M
   31 IH=M1-L1
      IK=M2-L2
      IL=M3-L3
      CALL STD(IH,IK,IL,ISER4)
      IZ=NN(ISER2)+NN(ISER4)
      IF(IZ.LE.0.or.iz.eq.4) goto 18
      IF(IZ.EQ.2) LL=LL+1
      IF(IZ.LE.5) LX=LX+1
      IF(IZ.EQ.5) SL=SL+1
      IF(IZ.GE.5) SX=SX+1
   18 continue
   17 continue
   16 continue
      NN(N1) =m11
      NN(N22)=m22
      NN(N33)=m33
      LL=LL/3
      LX=(LX-LL)/2
      SL=SL/2.
      V1= FLOAT(LL)/FLOAT(LX)
      V2= FLOAT(SL)/FLOAT(SX)
      dv=v1-v2
      write(99,24)A,H1,H2,H3,K1,K2,K3,L1,L2,L3,N1,N2,
     *N3,LL,LX,SL,SX,V1,V2,dv,arg
      nss=nss+LX
      mss=mss+SX
      if(dv.gt.dvv) mum=mum+1
      write(2,28)a,ires,it,n1,n2,n3,v1,v2,dv,arg
   28 format(f8.3,3i5,2i6,3f8.4,2x,f6.2)
   24 format(f6.2,3(i5,2i4),7i5,3f8.4,F6.2)
   11 continue
   12 continue
      write(*,*)i-1,mum,ntrip,nss,mss
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

