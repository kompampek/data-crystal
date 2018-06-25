      Character*1 EL,CA,NI,OX
      Character*12 FILN
      integer*2 ih(25000),ik(25000),il(25000),nf(99)
      dimension x(399),y(399),z(399),xx(399),yy(399)
      dimension zz(399),e(25000)
      DATA CC0,f11,f22,f1,f2,sx,sy,sz,n0,ii/8*0.,2*0/
c      read(*,13)FILN
      open(1,file='iled.P1')
      read(*,*)INC,MS,thres
      thres3=thres/3.
      PI2=2.*acos(-1.0)
      do 8 i=1,24999
      read(1,*,END=7)jh,jk,jl,ns,ne
      if(ns.gt.ms) goto 8
      ii=ii+1
      ee=0.001*ne
      ix=jh*jk*jl
      if(ix.ne.0) goto 12
      ix=abs(jh*jk) +abs(jl*(abs(jh)+abs(jk)))
      if(ix.eq.0) ee=1.414*ee
   12 j=mod(ii,10)
      if(j.ne.1) goto 11
      f11=f11+ee
      f22=f22+ee*ee
      n0=n0+1
   11 e(ii)=ee
      f1=f1+ee
      f2=f2+ee*ee
      ih(ii)=jh
      ik(ii)=jk
      il(ii)=jl
    8 continue
    7 NREF=ii
      write(*,*)NREF
      CLOSE(1)
      f1=f1/NREF
      f2=f2/NREF
      f11=f11/n0
      f22=f22/n0
      f3=f2-f1*f1
      f33=f22-f11*f11
      write(*,*)f1,f2,f3
      close(1)
      read(*,13)FILN
   13 format(A)    
      open(1,file=FILN)
      READ(1,14)CA,NI,OX
   14 format(3A1)   
      do 1 i=1,399
      read(1,2,END=15)EL,xx(i),yy(i),zz(i)
    2 format(13X,A1,16X,3f8.3)
      nf(i)=6
      if(EL.eq.NI) nf(i)=7
      if(EL.eq.OX) nf(i)=8
      sx=sx+xx(i)
      sy=sy+yy(i)
    1 sz=sz+zz(i)
   15 NAT=i-2  
      CLOSE(1)
      sx=sx/NAT
      sy=sy/NAT
      sz=sz/NAT
      write(*,*)NAT,sx,sy,sz
C     center molecule at origin      
      do 3 i=1,NAT
      x(i)=xx(i)-sx
      y(i)=yy(i)-sy
    3 z(i)=zz(i)-sz
      PI=acos(-1.0)
      rad=PI/180.
      write(*,*),INC,NAT,NREF,N0
      CALL CPU_TIME(T0)
C
C     Spherical Polar rotations (theta,phi)
C
      do 4 i=0,359,INC
      ang=i*rad
      c1=cos(ang)
      s1=sin(ang)
      do 5 j=0,359,INC
      if(j.eq.0.and.i.ne.0) goto 5
      ang=j*rad
      c2=cos(ang)
      s2=sin(ang)
      do 6 k=1,NAT   
C     rotate +theta around z-axis
      x1 = x(k)*c1 -y(k)*s1
      y1 = x(k)*s1 +y(k)*c1
      z1 = z(k)
C     then +phi around x-axis
      y2 = y1*c2 -z1*s2
      z2 = y1*s2 +z1*c2
      x2 = x1
C     then backrotate -theta around z-axis
      xx(k) = ( x2*c1 +y2*s1)/11.615
      yy(k) = (-x2*s1 +y2*c1)/15.705  
    6 zz(k) = z2/39.1
C    
      do 16 npass=1,2
      mpass=10
      if(npass.eq.2) mpass=1
      e1=0.
      e2=0.
      ef=0.
      do 9 ii=1,NREF,mpass
      ah=PI2*ih(ii)
      ak=PI2*ik(ii)
      al=PI2*il(ii)
      a=0.
      b=0.
      do 10 jj=1,NAT
      arg=ah*xx(jj)+ak*yy(jj)+al*zz(jj)
      a=a+nf(jj)*cos(arg)
   10 b=b+nf(jj)*sin(arg)
      ej=sqrt((a*a+b*b)/NAT)
      e1=e1+ej
      e2=e2+ej*ej
    9 ef=ef+ej*e(ii)
      m0=n0
      if(npass.eq.2) m0=NREF
      e1=e1/m0
      ef=ef/m0
      sig=e2/m0 -e1*e1
      if(npass.eq.1) cx=(ef-e1*f11)/sqrt(f33*sig)
      if(npass.eq.2) cc=(ef-e1*f1)/sqrt(f3*sig)
      if(npass.eq.1.and.cx.lt.thres3) goto 5
      if(npass.eq.1) goto 16
      if(npass.eq.2.and.cc.lt.thres) goto 5
      write(*,17) i,j,cc,cx
      write(3,17)i,j,cc,cx
   17 format(2i4,2f9.5)    
      cc0=max(cc0,cc)
   16 CONTINUE    
    5 CONTINUE
    4 CONTINUE
      CALL CPU_TIME(T1)
      time = (t1-t0)/60.
      write(*,*)time,cc0
      END
      
      
 
