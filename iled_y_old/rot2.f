      Character*1 EL,CA,NI,OX
      real*8 a,b,e1,e2,ef,ej,f1,f2,f11,f22,PI2,ah,ak,al,arg,ang
      real*8 c1,s1,c2,s2,c3,s3,x1,x2,y1,y2,z1,z2,rad
      integer*2 ih(25000),ik(25000),il(25000),nf(99)
      dimension x(99),y(99),z(99),xx(99),yy(99),zz(99),e(25000)
      open(1,file='P1.E')
      read(*,*)I1,I2,I3
      sig2=0.
      CC0=0.
      f1=0
      f2=0
      PI2=2.*acos(-1.0)
      ii=0
      do 8 i=1,17904
      read(1,*)jh,jk,jl,ns,ne
      ii=ii+1
      ee=0.001*ne
      e(ii)=ee
      f1=f1+ee
      f2=f2+ee*ee
      ih(ii)=jh
      ik(ii)=jk
      il(ii)=jl
    8 continue
      NREF=ii
      f1=f1/NREF
      f2=f2/NREF
      f3=f2-f1*f1
      write(*,*)f1,f2,f3
      close(1)
      open(1,file='iled.pdb')
      READ(1,111)CA,NI,OX
  111 format(3A1)   
      sx=0.
      sy=0.
      sz=0.
      do 1 i=1,99
      read(1,2,END=112)EL,xx(i),yy(i),zz(i)
    2 format(13X,A1,16X,3f8.3)
      nf(i)=6
      if(EL.eq.NI) nf(i)=7
      if(EL.eq.OX) nf(i)=8
      sig2=sig2+nf(i)*nf(i)
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
      rad=PI2/360.
      write(*,*)NAT,NREF
      CALL CPU_TIME(T0)
C
C     NON-EULER ROTATIONS?
C
      do 4 i=I1-2,I1+2,2
      ang=i*rad
      c1=cos(ang)
      s1=sin(ang)
      do 5 j=I2-2,I2+2,2
      ang=j*rad
      c2=cos(ang)
      s2=sin(ang)
      do 6 k=I3-2,I3+2,2
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
      zz(ii) = (z2*c3 -x2*s3)/39.31
      xx(ii) = (z2*s3 +x2*c3)/11.615
    7 yy(ii) =  y2/15.705
C
      e1=0.
      e2=0.
      ef=0.
      do 9 ii=1,NREF
      ah=PI2*ih(ii)
      ak=PI2*ik(ii)
      al=PI2*il(ii)
      a=0.
      b=0.
      do 10 jj=1,NAT
      arg=ah*xx(jj)+ak*yy(jj)+al*zz(jj)
      a=a+nf(jj)*cos(arg)
   10 b=b+nf(jj)*sin(arg)
      ej=sqrt((a*a+b*b)/sig2)
      e1=e1+ej
      e2=e2+ej*ej
    9 ef=ef+ej*e(ii)
      e1=e1/NREF
      e2=e2/NREF
      ef=ef/NREF
      sig=e2-e1*e1
      if(i.eq.0.and.j.eq.0.and.k.eq.0) write(*,*)e1,e2,sig
      cc=(ef-e1*f1)/sqrt(f3*sig)
      write(*,*) i,j,k,cc
      write(3,33)i,j,k,cc
   33 format(3i5,f8.4)    
      cc0=max(cc0,cc)   
    6 CONTINUE
    5 CONTINUE
    4 CONTINUE
      CALL CPU_TIME(T1)
      time = (t1-t0)/60.
      write(*,*)time,cc0
      END
      
      
 
