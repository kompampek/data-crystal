      PARAMETER (IE=10000)
      CHARACTER an*6,bn*6
      byte JH(IE,24),JK(IE,24),JL(IE,24),NR(IE),IR(24,9)
      integer*4 LOC(99999)
      dimension p(-IE:IE)
      COMMON /HKL/ LAUE,MH,MK,ML,NK,NL,IHR,IKR,LOC
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
   26 format(' Ebig,Esmall,NTRIPS,delv : ',$)
      read(*,*)EBIG,ESMALL,NTRIP,dvv
      Ebig = Ebig -0.005
      Esmall = Esmall + 0.005
      open(1,file=an,form='formatted')
      do 5 I=1,99999
      read(1,*,end=6)ih,ik,il,ns,ne
      en = 0.001*ne
      if(en.ge.EBIG) LIM1=I
      if(en.ge.ESMALL) LIM2=I
      mh=max(ih,mh)
      mk=max(ik,mk)
      ml=max(il,ml)
      nk=min(ik,nk)
    5 nl=min(il,nl)
    6 rewind(1)
      NUM=I-1
      LIM3=NUM-LIM2
      write(*,23)NUM,LIM1,EBIG,LIM3,ESMALL,mh,mk,ml,nk,nl
   23 format(i5,' total data,',i5,' Es >',f6.2,',',i6,' Es <',f6.2,5i5)  
      IHR=MH+1
      IKR=MK-NK+1
      do 7 I=1,NUM
      read(1,*)ih,ik,il,ns,ne,a,b,c,np
      if(i.gt.IE) goto 8
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
      do 11 i=1,NTRIP
      read(1,*,end=12)a,ires,it,n1,n2,n3
      arg=cos(p(n1)+p(n2)+p(n3))
      if(it.eq.6) arg=-arg
      I1=jh(n1,1)
      I2=jk(n1,1)
      I3=jl(n1,1)
      m2=abs(n2)
      j2=m2/n2
      m3=abs(n3)
      j3=m3/n3
      do 13 j=1,NR(M2)
      JJ1=j2*jh(m2,j)
      JJ2=j2*jk(m2,j)
      JJ3=j2*jl(m2,j)
      do 14 k=1,NR(M3)
      KK1=j3*jh(m3,k)
      KK2=j3*jk(m3,k)
      KK3=j3*jl(m3,k)
      ix=abs(i1+jj1+kk1)+abs(i2+jj2+kk2)+abs(i3+jj3+kk3)
   14 if(ix.eq.0) goto 15
   13 continue
      ix=1
   15 if(ix.ne.0) write(*,*)99

      LLL=0
      LLX=0
      LLS=0
      LSX=0
      do 16 J=1,LIM1
      do 17 k=1,NR(J)
      do 18 l=-1,1,2
      L1=l*jh(j,k)
      L2=l*jk(j,k)
      L3=l*jl(j,k)
      IH=I1+L1
      IK=I2+L2
      IL=I3+L3
      CALL STD(IH,IK,IL,ISER1)
      IH=JJ1-L1
      IK=JJ2-L2
      IL=JJ3-L3
      CALL STD(IH,IK,IL,ISER2)
      IF(ISER1.LE.0.OR.ISER2.LE.0) GOTO 19
      IF(ISER1.LE.LIM1.AND.ISER2.LE.LIM1) LLL=LLL+1
      IF(ISER1.GE.LIM2.AND.ISER2.LE.LIM1) LLS=LLS+1
      IF(ISER1.LE.LIM1.AND.ISER2.GE.LIM2) LLS=LLS+1
      IF(ISER1.LE.LIM1.OR.ISER2.LE.LIM1) LLX=LLX+1
      IF(ISER1.GE.LIM2.OR.ISER2.GE.LIM2) LSX=LSX+1
   19 IH=KK1+L1
      IK=KK2+L2
      IL=KK3+L3
      CALL STD(IH,IK,IL,ISER3)
      IH=JJ1-L1
      IK=JJ2-L2
      IL=JJ3-L3
      CALL STD(IH,IK,IL,ISER4)
      IF(ISER3.LE.0.OR.ISER4.LE.0) GOTO 20
      IF(ISER3.LE.LIM1.AND.ISER4.LE.LIM1) LLL=LLL+1
      IF(ISER3.GE.LIM2.AND.ISER4.LE.LIM1) LLS=LLS+1
      IF(ISER3.LE.LIM1.AND.ISER4.GE.LIM2) LLS=LLS+1
      IF(ISER3.LE.LIM1.OR.ISER4.LE.LIM1) LLX=LLX+1
      IF(ISER3.GE.LIM2.OR.ISER4.GE.LIM2) LSX=LSX+1
   20 IH=I1+L1
      IK=I2+L2
      IL=I3+L3
      CALL STD(IH,IK,IL,ISER5)
      IH=KK1-L1
      IK=KK2-L2
      IL=KK3-L3
      CALL STD(IH,IK,IL,ISER6)
      IF(ISER5.LE.0.OR.ISER6.LE.0) GOTO 18
      IF(ISER5.LE.LIM1.AND.ISER6.LE.LIM1) LLL=LLL+1
      IF(ISER5.GE.LIM2.AND.ISER6.LE.LIM1) LLS=LLS+1
      IF(ISER5.LE.LIM1.AND.ISER6.GE.LIM2) LLS=LLS+1
      IF(ISER5.LE.LIM1.OR.ISER6.LE.LIM1) LLX=LLX+1
      IF(ISER5.GE.LIM2.OR.ISER6.GE.LIM2) LSX=LSX+1
   18 continue
   17 continue
   16 continue
      LLL=LLL/3
      LLS=LLS/2
      LLX=(LLX-LLL)/2
      LSX=LSX-LLS
      V1 = FLOAT(LLL)/FLOAT(LLX)
      V2 = FLOAT(LLS)/FLOAT(LSX)
      d1=v1*lll
      d2=v2*lls
      dv=v1-v2
      write(99,24)A,I1,I2,I3,JJ1,JJ2,JJ3,KK1,KK2,KK3,N1,N2,
     *N3,LLL,LLX,LLS,LSX,V1,V2,dv,arg,d1,d2
      nss=nss+LLX
      mss=mss+LSX
      if(dv.gt.dvv) mum=mum+1
      write(2,28)a,ires,it,n1,n2,n3,v1,v2,dv,arg,d1,d2
   28 format(f8.3,3i5,2i6,3f8.4,f6.2,2f5.1)
   24 format(f6.2,3(i5,2i4),7i5,3f8.3,f6.2,2f5.1)
   11 continue
   12 continue
      write(*,*)i-1,mum,nss,mss
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

