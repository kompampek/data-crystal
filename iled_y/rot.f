      Character*1 EL,CA,NI,OX
      Character*12 FILN
      integer*2 ih(25000),ik(25000),il(25000),nf(99)
      dimension x(399),y(399),z(399),xx(399),yy(399)
      dimension zz(399),e(25000)
      open(1,file='cell')
      read(1,*)ax,bx,cx
      close(1)
      open(1,file='P1.E')
      read(*,*)INC,MS,thres
      thres3=thres/3.
      CC0=0.
      f11=0.
      f22=0.
      n0=0
      f1=0
      f2=0
      PI2=2.*acos(-1.0)
      ii=0
      do 8 i=1,25000
      read(1,*,END=109)jh,jk,jl,ns,ne
      if(ns.gt.ms) goto 8
      ii=ii+1
      ee=0.001*ne
      ix=jh*jk*jl
      if(ix.ne.0) goto 82
      ix=abs(jh*jk) +abs(jl*(abs(jh)+abs(jk)))
      if(ix.eq.0) ee=1.414*ee
   82 j=mod(ii,10)
      if(j.ne.1) goto 81
      f11=f11+ee
      f22=f22+ee*ee
      n0=n0+1
   81 e(ii)=ee
      f1=f1+ee
      f2=f2+ee*ee
      ih(ii)=jh
      ik(ii)=jk
      il(ii)=jl
    8 continue
  109 NREF=ii
      f1=f1/NREF
      f2=f2/NREF
      f11=f11/n0
      f22=f22/n0
      f3=f2-f1*f1
      f33=f22-f11*f11
      write(*,*)f1,f2,f3
      close(1)
      read(*,110)FILN
  110 format(A)    
      open(1,file=FILN)
      READ(1,111)CA,NI,OX
  111 format(3A1)   
      sx=0.
      sy=0.
      sz=0.
      do 1 i=1,399
      read(1,2,END=112)EL,xx(i),yy(i),zz(i)
    2 format(13X,A1,16X,3f8.3)
      nf(i)=6
      if(EL.eq.NI) nf(i)=7
      if(EL.eq.OX) nf(i)=8
      sx=sx+xx(i)
      sy=sy+yy(i)
    1 sz=sz+zz(i)
  112 NAT=i-1  
      CLOSE(1)
      sx=sx/NAT
      sy=sy/NAT
      sz=sz/NAT
      do 3 i=1,NAT
      x(i)=xx(i)-sx
      y(i)=yy(i)-sy
    3 z(i)=zz(i)-sz
      SUM=99999.
      PI=acos(-1.0)
      rad=PI/180.
      write(*,*),INC,NAT,NREF,N0
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
      ang=k*rad
      c3=cos(ang)
      s3=sin(ang)
      do 7 ii=1,NAT
C     rotate on z     
      x1 = x(ii)*c1 -y(ii)*s1
      y1 = x(ii)*s1 +y(ii)*c1
      z1 = z(ii)
C     rotate on x      
      y2 = y1*c2 -z1*s2
      z2 = y1*s2 +z1*c2
      x2 = x1
C     rotate on y      
      zz(ii) = (z2*c3 -x2*s3)/cx
      xx(ii) = (z2*s3 +x2*c3)/ax
    7 yy(ii) =  y2/bx
C    
      do 99 npass=1,2
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
      e2=e2/m0
      ix=abs(1)+abs(j)+abs(k)
      ef=ef/m0
      sig=e2-e1*e1
      if(npass.eq.1) cc=(ef-e1*f11)/sqrt(f33*sig)
      if(npass.eq.2) cc=(ef-e1*f1)/sqrt(f3*sig)
      if(npass.eq.1.and.cc.lt.thres3) goto 6
      if(npass.eq.1) goto 99
      if(npass.eq.2.and.cc.lt.thres) goto 6
      write(*,113) i,j,k,cc
      write(3,113)i,j,k,cc
  113 format(3i4,f9.6)    
      cc0=max(cc0,cc)
   99 CONTINUE    
    6 CONTINUE
    5 CONTINUE
    4 CONTINUE
      CALL CPU_TIME(T1)
      time = (t1-t0)/60.
      write(*,*)time,cc0
      END
      
      
 
