      Character*1 EL,CA,NI,OX
      Character*12 FILN
      integer*2 ih(25000),ik(25000),il(25000),nf(99)
      integer*2 nx(99),ny(99),nz(99),NARG
      dimension x(99),y(99),z(99),xx(99),yy(99),zz(99)
      dimension e(25000),cs(-32768:32768),sn(-32768:32768)
      rd=acos(-1.0)/32768.
      do 18 i=-32768,32768
      arg=i*rd
      cs(i)=cos(arg)
   18 sn(i)=sin(arg)   
      open(1,file='cell')
      read(1,*)ax,bx,cx
      close(1)
      open(1,file='iled.P1')
      write(*,22)
   22 format('inc,ms,thres,NAT ')   
      read(*,*)INC,MS,thres,NAT
      thres3=thres/3.
      CC0=0.
      f11=0.
      f22=0.
      n0=0
      f1=0
      f2=0
      ii=0
      do 8 i=1,25000
      read(1,*,END=109)jh,jk,jl,ns,ne
      if(ns.gt.ms) goto 8
      ii=ii+1
      ee=0.001*ne
c      ee=ee*ee
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
      write(*,199)n0,ii,f1,f2,f3,f11,f22,f33
  199 format(2i6,6f8.4)    
      close(1)
      read(*,110)FILN
  110 format(A)    
      open(1,file=FILN)
      READ(1,111)CA,NI,OX
  111 format(3A1)   
      sx=0.
      sy=0.
      sz=0.
      do 1 i=1,NAT
      read(1,2,END=112)EL,xx(i),yy(i),zz(i)
    2 format(13X,A1,16X,3f8.3)
      nf(i)=6
      if(EL.eq.NI) nf(i)=7
      if(EL.eq.OX) nf(i)=8
      sx=sx+xx(i)
      sy=sy+yy(i)
    1 sz=sz+zz(i)
  112 CLOSE(1)
      sx=sx/NAT
      sy=sy/NAT
      sz=sz/NAT
      do 3 i=1,NAT
      x(i)=xx(i)-sx
      y(i)=yy(i)-sy
    3 z(i)=zz(i)-sz
      write(*,*),INC,NAT,NREF,N0,sx,sy,sz
      CALL CPU_TIME(T0)
C
C     NON-EULER ROTATIONS?
C
      rad=acos(-1.0)/180.
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
      z3=(z2*c3 -x2*s3)/cx
      x3=(z2*s3 +x2*c3)/ax
      y3=y2/bx
      nz(ii) = 65536*z3
      nx(ii) = 65536*x3
    7 ny(ii) = 65536*y3
C    
      do 99 npass=1,2
      mpass=10
      if(npass.eq.2) mpass=1
      e1=0.
      e2=0.
      ef=0.
      do 9 ii=1,NREF,mpass
      jh=ih(ii)
      jk=ik(ii)
      jl=il(ii)
      a=0.
      b=0.
      do 10 jj=1,NAT
      narg=jh*nx(jj)+jk*ny(jj)+jl*nz(jj)
      a=a+nf(jj)*cs(narg)
   10 b=b+nf(jj)*sn(narg)
      ej=sqrt((a*a+b*b)/NAT)
c      ej=ej*ej
      e1=e1+ej
      e2=e2+ej*ej
    9 ef=ef+ej*e(ii)
      m0=n0
      if(npass.eq.2) m0=NREF
      e1=e1/m0
      ef=ef/m0
      sig=e2/m0 -e1*e1
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
      
      
 
