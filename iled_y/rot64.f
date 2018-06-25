      BYTE EL,CAR,NIT,OXY
      integer*2 ih(7000),ik(7000),il(7000),nf(82)
      integer*2 nh,nk,nl,N1,N2,N3,NS,nx(82),ny(82),nz(82)
      dimension e(7000),cs(-32768:32768),sn(-32768:32768)
      dimension x(82),y(82),z(82),xx(82),yy(82),zz(82)
      open(1,file='2221.E')
      ntot=0
      neg=0
      read(*,*)thres,INC
      thres2=thres/2.
      PI=acos(-1.0)
      PI2=PI+PI
      rad=PI/32768.
      do 11 i=-32768,32768
      arg=i*rad
      cs(i)=cos(arg)
   11 sn(i)=sin(arg)   
      CC0=0.
      f11=0.
      f22=0.
      n0=0
      n1=0
      f1=0
      f2=0
      do 8 i=1,6685
      read(1,*)ih(i),ik(i),il(i),ns,ne
      e(i)=0.001*ne      
      ix=ih(i)*ik(i)*il(i)
      if(ix.eq.0) goto 8
      j=mod(i,10)
      if(j.ne.1) goto 81
      f11=f11+e(i)**2
      f22=f22+e(i)**4
      n0=n0+1
   81 f1=f1+e(i)**2
      f2=f2+e(i)**4
      n1=n1+1
    8 CONTINUE
      write(*,*)n1,n0
      f1=f1/n1
      f2=f2/n1
      f11=f11/n0
      f22=f22/n0
      close(1)
      open(1,file='64.pdb')
      read(1,111)CAR,NIT,OXY
  111 format(3a1)    
      sx=0.
      sy=0.
      sz=0.
      do 1 i=1,62
      read(1,2)xx(i),yy(i),zz(i),EL
      if(EL.EQ.CAR) it=6
      if(EL.eq.NIT) it=7
      if(EL.eq.OXY) it=8
      nf(i)=it
    2 format(30x,3f8.3,23x,a1)
      sx=sx+xx(i)
      sy=sy+yy(i)
    1 sz=sz+zz(i)
      CLOSE(1)
      sx=sx/62.
      sy=sy/62.
      sz=sz/62.
      do 3 i=1,62
      x(i)=xx(i)-sx
      y(i)=yy(i)-sy
    3 z(i)=zz(i)-sz
      SUM=99999.
      rad=PI/180.
      CALL CPU_TIME(T0)
C
C     NON-EULER ROTATIONS?
C
      do 4 i=0,179,INC
      ang=i*rad
      c1=cos(ang)
      s1=sin(ang)
      do 5 j=0,179,INC
      ang=j*rad
      c2=cos(ang)
      s2=sin(ang)
      do 6 k=0,179,INC
      ntot=ntot+1
      ang=k*rad
      c3=cos(ang)
      s3=sin(ang)
      do 7 ii=1,62
C     rotate on z     
      x1 = x(ii)*c1 -y(ii)*s1
      y1 = x(ii)*s1 +y(ii)*c1
      z1 = z(ii)
C     rotate on x      
      y2 = y1*c2 -z1*s2
      z2 = y1*s2 +z1*c2
      x2 = x1
C     rotate on y      
      zz(ii) = 65536*((z2*c3 -x2*s3)/38.303)
      xx(ii) = 65536*((z2*s3 +x2*c3)/25.79)
    7 yy(ii) = 65536*(y2/35.479)
C    
      do 99 npass=1,2
      mpass=10
      if(npass.eq.2) mpass=1
      e1=0.
      e2=0.
      ef=0.
      do 9 ii=1,6685,mpass
      ix=ih(ii)*ik(ii)*il(ii)
      if(ix.eq.0) goto 9
      ah=ih(ii)
      ak=ik(ii)
      al=il(ii)
      a1=0.
      b1=0.
      a2=0
      b2=0
      a3=0
      b3=0
      a4=0
      b4=0
      do 10 jj=1,62
      nj=nf(jj)
      x1=ah*xx(jj)
      x2=ak*yy(jj)
      x3=al*zz(jj)
      NS= x1+x2+x3
      a1=a1+cs(NS)
      b1=b1+sn(NS)
      NS=-x1-x2+x3
      a2=a2+cs(NS)
      b2=b2+sn(NS)
      NS=-x1+x2-x3
      a3=a3+cs(NS)
      b3=b3+sn(NS)
      NS= x1-x2-x3
      a4=a4+cs(NS)
   10 b4=b4+sn(NS)
      a1=a1*a1 +b1*b1
      a2=a2*a2 +b2*b2
      a3=a3*a3 +b3*b3
      a4=a4*a4 +b4*b4
      av=a1+a2+a3+a4
      e1=e1+av
      e2=e2+av**2
      ef=ef+av*e(ii)**2
    9 CONTINUE
      m0=n0
      if(npass.eq.2) m0=n1
      e1=e1/m0
      e2=e2/m0
      ef=ef/m0
      a1=ef-e1*f11
      a2=ef-e1*f1
      a3=f22-f11*f11
      a4=f2-f1*f1
      a5=e2-e1*e1
      if(npass.eq.1) cold=(a1)/sqrt((a3)*(a5))
      if(npass.eq.1.and.cold.le.0.0) neg=neg+1
      if(npass.eq.2) cc=(a2)/sqrt((a4)*(a5))
      if(npass.eq.1.and.cold.lt.thres2) goto 6
      if(npass.eq.1) goto 99
      if(npass.eq.2.and.cc.lt.thres) goto 6
      write(*,33) i,j,k,cold,cc
      write(3,33)i,j,k,cold,cc
   33 format(3i5,2f8.4)   
   99 CONTINUE    
    6 CONTINUE
    5 CONTINUE
    4 CONTINUE
      CALL CPU_TIME(T1)
      time = (t1-t0)/60.
      write(*,*)neg,ntot,time
      END
