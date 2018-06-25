      Character*1 EL,CA,NI,OX
      Character*12 FILN
      integer*2 ih(25000),ik(25000),il(25000),nf(3)
      integer*2 nx(399),ny(399),nz(399),NARG,n1(3),n2(3)
      dimension x(399),y(399),z(399),e(25000)
      dimension cs(-32768:32768),sn(-32768:32768)
      DATA CC0,f11,f22,f1,f2,sx,sy,sz,n0,ii/8*0.,2*0/
      data n1,n2,nf,mh,mk,ml/1,61,67,60,66,84,6,7,8,3*0/
      ccmax=0.
      ax=11.615
      bx=15.705
      cx=39.31
      rd=acos(-1.0)/32768.
      do 18 i=-32768,32768
      arg=i*rd
      cs(i)=cos(arg)
   18 sn(i)=sin(arg)
      open(1,file='iled.P1')
      write(*,22)
   22 format(' INC,MS,thres,NAT  ',$)   
      read(*,*)INC,MS,thres,NAT
      thres3=thres/3.
      do 8 i=1,24999
      read(1,*,END=7)jh,jk,jl,ns,ne
      mh=max(mh,jh)
      mk=max(mk,jk)
      ml=max(ml,jl)
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
      do 1 i=1,NAT
      read(1,2,END=15)EL,x(i),y(i),z(i)
    2 format(13X,A1,16X,3f8.3)
      sx=sx+x(i)
      sy=sy+y(i)
    1 sz=sz+z(i)
   15 CLOSE(1)
      sx=sx/NAT
      sy=sy/NAT
      sz=sz/NAT
      write(*,*)NAT,sx,sy,sz,mh,mk,ml
C     center molecule at origin      
      do 3 i=1,NAT
      x(i)=x(i)-sx
      y(i)=y(i)-sy
    3 z(i)=z(i)-sz
      rad=acos(-1.0)/180.
      write(*,*),INC,NAT,NREF,N0
      CALL CPU_TIME(T0)
C
C     Spherical Polar rotations (theta,phi)
C
      do 4 i=0,179,INC
      ang=i*rad
      c1=cos(ang)
      s1=sin(ang)
      do 5 j=0,179,INC
      ang=j*rad
      c2=cos(ang)
      s2=sin(ang)
      do 20 ij=0,179,INC
      ang=ij*rad
      c3=cos(arg)
      s3=sin(arg)
      do 6 k=1,NAT   
C     rotate +theta around z-axis
      x1 = x(k)*c1 -y(k)*s1
      y1 = x(k)*s1 +y(k)*c1
      z1 = z(k)
C     then +phi around x-axis
      y2 = y1*c2 -z1*s2
      z2 = y1*s2 +z1*c2
      x2 = x1
C     then rotate around y-axis
      zz = (z2*c3 -x2*s3)/cx
      xx = (z2*s3 +x2*c3)/ax
      yy = y2/bx
      nx(k)=65536*xx
      ny(k)=65536*yy
    6 nz(k)=65536*zz  
C    
      do 16 npass=1,2
      mpass=10
      if(npass.eq.2) mpass=1
      e1=0.
      e2=0.
      ef=0.
      do 9 ii=1,NREF,mpass
      jh=ih(ii)
      jk=ik(ii)
      jl=il(ii)
      aa=0.
      bb=0.
      do 19 kk=1,3
      mf=nf(kk)
      m1=n1(kk)
      m2=n2(kk)
      a=0.
      b=0.
      do 10 jj=m1,m2
      NARG=jh*nx(jj)+jk*ny(jj)+jl*nz(jj)
      a=a+cs(narg)
   10 b=b+sn(narg)
      aa=aa+a*mf
   19 bb=bb+b*mf
      ej=sqrt((aa*aa+bb*bb)/NAT)
      e1=e1+ej
      e2=e2+ej*ej
    9 ef=ef+ej*e(ii)
      m0=n0
      if(npass.eq.2) m0=NREF
      e1=e1/m0
      ef=ef/m0
      sig=e2/m0 -e1*e1
      if(npass.eq.1) cc0=(ef-e1*f11)/sqrt(f33*sig)
      if(npass.eq.2) cc=(ef-e1*f1)/sqrt(f3*sig)
      if(npass.eq.1.and.cc0.lt.thres3) goto 5
      if(npass.eq.1) goto 16
      if(npass.eq.2.and.cc.lt.thres) goto 5
      write(*,17) i,j,ij,cc,cc0
      write(3,17)i,j,ij,cc,cc0
   17 format(3i4,2f9.6)    
      ccmax=max(ccmax,cc)
   16 CONTINUE
   20 continue   
    5 CONTINUE
    4 CONTINUE
      CALL CPU_TIME(T1)
      time = (t1-t0)/60.
      write(*,*)time,ccmax,n1,n2
      END
      
      
 
