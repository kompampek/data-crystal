C     SIR data
      CHARACTER Car,Nit,OX,SU,FE,ALP(7),n1,n2   !!YV
      byte NB(6),nj(5200)                        !!YV
C !!YV     byte Car,Nit,OX,SU,FE,ALP(7),NB(6),n1,n2,nj(5200)
      integer*2 ih,ik,il,ns,ne,np
      dimension x(5200),y(5200),z(5200),bi(5200),g(5200),cell(6)
      dimension ih1(4108),ik1(4108),il1(4108),FN1(4108),PN1(4108) !!YV
      dimension ff(9,8),fj(8),CNT(0:8),fp(8)
      data Car,Nit,OX,SU,FE/'C','N','O','S','E'/
      data ALP,NB/'B','G','D','E','H','Z',' ',5,8,13,20,20,20/
      write(*,14)
   14 format(' RESOLUTION OF DATA IN ANG = ',$)
      read(*,*)RES
      t0=secnds(0.0)
      e2=0.0   !!YV
      nume=0   !!YV
      f2=0.0   !!YV
      OPEN(1,file='tox2_external')          !!YV
      read(1,88,end=2) ih1,ik1,il1,FN1,PN1  !!YV
CCCCCCCCC      write (*,88) ih1,ik1,il1,FN1,PN1      !!YV  
   88 format(3x,I2,4x,I2,4x,I2,5x,F11.6,3x,F6.2)          !!YV
      close(1)                              !!YV 
      OPEN(1,file='sas.dat')
      read(1,*)cell
      write(*,*)cell   !!YV
      read(1,*)ntypes
      do 6 i=1,ntypes
    6 read(1,*) (ff(j,i),j=1,9)
      read(1,*) (fp(i),i=1,ntypes)
      do 10 i=1,ntypes
      a=0
CCC      write(*,*) 'fp==', fp(i)  !!YV
      ff(9,i)=ff(9,i)+fp(i)
      do 11 j=1,9,2 
   11 a=a+ ff(j,i)
   10 write(*,*)a
      CLOSE(1)
      PI=acos(-1.0)
      P2=2.0*PI
      RAD=P2/360.
      CA=COS(RAD*cell(4))
      CB=COS(RAD*cell(5))
      CG=COS(RAD*cell(6))
      SA=SIN(RAD*cell(4))
      SB=SIN(RAD*cell(5))
      SG=SIN(RAD*cell(6))
      V=cell(1)*cell(2)*cell(3)*SQRT(1.-CA*CA-CB*CB-CG*CG+2*CA*CB*CG)
      AS=cell(2)*cell(3)*SA/V
      BS=cell(1)*cell(3)*SB/V
      CS=cell(1)*cell(2)*SG/V
      CAS=(CB*CG-CA)/(SB*SG)
      CBS=(CA*CG-CB)/(SA*SG)
      CGS=(CA*CB-CG)/(SA*SB)
      aa1=2.0*as*bs*cgs
      aa2=2.0*as*cs*cbs     
      OPEN(1,file='X.pdb')
      do 1 i=1,1017
      read(1,8,end=2)atom,n1,n2,xx,yy,zz,goc,bis
    8 format(A4,9x,2a1,15x,3f8.3,2f6.2)
      ix=0
      if(n1.eq.Car) ix=5
      if(n1.eq.FE) ix=1
      if(n1.eq.Nit) ix=4
      if(n1.eq.OX) ix=3
      if(n1.eq.SU) ix=2
      if(ix.eq.0) write(*,8)atom,n1,i
      CNT(ix)=cnt(ix)+goc
      if(ix.eq.0) ix=8
      nj(i)=ix
      g(i)=goc
      bis=15.
      if(n1.eq.OX.and.n2.eq.ALP(7)) bis=bis+5
      if(n1.eq.OX.and.n2.ne.ALP(7)) bis=bis+10
      do 12 j=1,6
      if(n1.eq.Car.and.n2.eq.alp(j)) bis=bis+nb(j)
   12 if(n1.eq.Nit.and.n2.eq.alp(j)) bis=bis+nb(j)
      bi(i)= bis
      x(i)=xx/cell(1)
      y(i)=yy/cell(2)
      z(i)=zz/cell(3)
    1 write(11,111)i,x(i),y(i),z(i),bi(i),nj(i)  
  111 format(i6,3f8.4,f6.1,i3)
    2 natoms=i-1
      RE=0.5/RES
      write(*,*)(cnt(i),i=0,ntypes),natoms
      CLOSE(1)
      d1=-1  
      do 33   k=1,4108
         write(*,*) k,ih1,ik1,il1 
33     continue	   
       do 31 k=1,4108                                        !!YV    
       ih=ih1(k)                                             !!YV
       ik=ik1(k)                                             !!YV
       il=il1(k)                                             !!YV
       write(*,*) k,ih1,ik1,il1          !!!YV
C!!YV      do 33 ih=0,26
      ah1=(ih/cell(1))**2
CCCCC!!YV     do 32 ik=0,51
      ah2=ah1 +(ik/cell(2))**2
CC!!!YV      do 31 il=0,20
      if(ih.eq.0.and.ik.eq.0.and.il.eq.0) goto 31
C     filter out axial extinctions
      ix=ih*(ik+il) + ik*il
      iy=mod((ih+ik+il),2)
      if(ix.eq.0.and.iy.eq.1) goto 31
      ah3=ah2 +(il/cell(3))**2
      s = 0.5*sqrt(ah3)
      if(s.gt.RE) goto 31
      ss=s   !!YV
      ns=10000.*s
      s=-s*s
      do 5 i=1,ntypes
      f2=ff(1,i)*exp(ff(2,i)*s) + ff(3,i)*exp(ff(4,i)*s) + ff(9,i)
    5 fj(i)= f2 + ff(5,i)*exp(ff(6,i)*s)+ ff(7,i)*exp(ff(8,i)*s)
      ah=(ih-ik)/4.
      ak=(ik-il)/4.
      al=(il-ih)/4.
      a=0
      b=0  
      sz2=0
C!!YV     do 7 i=2,natoms
      do 7 i=1,natoms             !!YV  
      if (nj(i).eq.8) go to 7     !!YV
      fh =fj(nj(i))*g(i)*exp(bi(i)*s)
      sz2=sz2+fh*fh
      a1 = p2*(ih*x(i)-ah)
      a2 = p2*(ik*y(i)-ak)
      a3 = p2*(il*z(i)-al)
      a = a + fh*cos(a1)*cos(a2)*cos(a3)
    7 b = b - fh*sin(a1)*sin(a2)*sin(a3)
      PN = atan2(b,a)
      FN = 4.0*sqrt(a*a+b*b)
      esq=FN*FN/(4*sz2)
      U=sqrt(esq)
C     write(*,*)'s,fh.FN.E,sz2',s,fh,FN,U,sz2  !!YV
      NE=1000*U              !!YV
      NS=10000*ss            !!YV
      e2=e2+esq
      fh=fj(nj(1))*g(1)*exp(bi(1)*s)
      a1 = p2*(ih*x(1)-ah)
      a2 = p2*(ik*y(1)-ak)
      a3 = p2*(il*z(1)-al)
      a = fh*cos(a1)*cos(a2)*cos(a3)
      b =-fh*sin(a1)*sin(a2)*sin(a3)
      fh = 4.0*sqrt(a*a+b*b)            !!YV
      PH=atan2(b,a)
      PM=mod((2*PH-PN),P2)
      if(PM.gt.PI) PM=PM-P2
      if(PM.lt.-PI) PM=PM+P2
      d1=-d1
      ix=ih*ik*il
      if(ix.eq.0) d1=1
      nume=nume+1
      FA=4.*sqrt(a*a+b*b)
      ARG=0.5*FN*FA/(sz2)
      ar=mod((PN-PM+P2+P2),P2)
      if(ar.gt.PI) ar=ar-P2
      if(ar.gt.PI) ar=ar-P2
      amp=FN*abs(sin(0.5*ar))
C     FN=native amplitude, PN=native phase
C     PH=HEAVY phase, PM=false phase, d1&d2 +/-1
      if(amp.lt.1.) amp=1.
      if((ih*ik*il).eq.0) amp=0.
      WT=1.-exp(-ESQ)
      if(WT.lt.0.001) WT=0.001
      a=cos(PN)+cos(PM)
      b=sin(PN)+sin(PM)
      ar=atan2(b,a)
C !YV    write(3,9)ih,ik,il,FN,PN,PM,ar,AMP,WT
C !YV  9 format(3i4,f8.2,3f7.3,f8.2,f6.3)
CCCCCCCCC      write(3,9)ih,ik,il,NS,NE,fh,FN,PH,PN    !!YV
      write(3,9)ih,ik,il,NS,NE,FN1(k),FN,PN1(k),PN    !!YV
      write(*,9)ih,ik,il,NS,NE,FN1(k),FN,PN1(k),PN    !!YV
     
    9 format(3i4,2x,i5,2x,i4,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2)            !!YV
   31 continue
C!!V   32 continue
C !!YV  33 continue
      e2=e2/nume
      write(*,*)nume,e2
      end
