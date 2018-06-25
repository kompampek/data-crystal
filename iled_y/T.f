      Character*1 CA,NI,OX
      BYTE A(30)
      dimension x(99),y(99),z(99),xx(88),yy(88),zz(88)
      read(*,*)NTHE,NPHI,NCHI
      open(1,file='ILED.ent')
      READ(1,2)CA,NI,OX
    2 format(30A1,3f8.3) 
      sx=0.
      sy=0.
      sz=0.
      do 4 i=1,84
      read(1,2)A,x(i),y(i),z(i)
      sx=sx+x(i)
      sy=sy+y(i)
      sz=sz+z(i)
    4 CONTINUE
C     center molecule at origin 
      do 3 i=1,84
      x(i)=x(i)-sx/84.
      y(i)=y(i)-sy/84.
    3 z(i)=z(i)-sz/84.
      rad=acos(-1.0)/180.
      read(*,*)inc,t
      t3=3.*t
      CALL CPU_TIME(T0)
C
C     Euler rotations (theta,phi,chi)
C
      ang=NTHE*rad
      c1=cos(ang)
      s1=sin(ang)
      ang=NPHI*rad
      c2=cos(ang)
      s2=sin(ang)
      ang=NCHI*rad
      c3=cos(ang)
      s3=sin(ang)
      a11 =-s1*c2*s3 +c1*c3
      a12 = c1*c2*s3 +s1*c3
      a13 = s2*s3
      a21 =-s1*c2*c3 -c1*s3
      a22 = c1*c2*c3 -s1*s3
      a23 = s2*c3
      a31 = s1*s2
      a32 =-c1*s2
      a33 = c2
      do 1 k=1,84
C      read(1,2)A
C      
      xx(k) = x(k)*a11 +y(k)*a21 + z(k)*a31
      yy(k) = x(k)*a12 +y(k)*a22 + z(k)*a32
      zz(k) = x(k)*a13 +y(k)*a23 + z(k)*a33
c      
      xx(k) = x(k)*a11 +y(k)*a12 + z(k)*a13
      yy(k) = x(k)*a21 +y(k)*a22 + z(k)*a23
      zz(k) = x(k)*a31 +y(k)*a32 + z(k)*a33
    1 CONTINUE
      rmax=0.
      rmin=99.
      nu=0
      do 5 i=1,83
      do 11 j=i+1,84
      r=sqrt((xx(i)-xx(j))**2 +(yy(i)-yy(j))**2 +(zz(i)-zz(j))**2 )
      if(r.gt.1.7) goto 11
      nu=nu+1
      rmin=min(rmin,r)
      rmax=max(rmax,r)
   11 continue
    5 continue
      write(*,*)nu,rmin,rmax
      read(*,*)junk   
C    
      rad=acos(-1.0)/180.
      do 6 i=0,359,INC
      ang= i*rad
      c1=cos(ang)
      s1=sin(ang)
      do 7 j=0,359,INC
      ang= j*rad
      c2=cos(ang)
      s2=sin(ang)
      do 8 k=0,359,INC
      ang= k*rad
      c3=cos(ang)
      s3=sin(ang)
c      
      a11 =-s1*c2*s3 +c1*c3
      a12 = c1*c2*s3 +s1*c3
      a13 = s2*s3
      a21 =-s1*c2*c3 -c1*s3
      a22 = c1*c2*c3 -s1*s3
      a23 = s2*c3
      a31 = s1*s2
      a32 =-c1*s2
      a33 = c2
c      
      sm1=0.
      sm2=0.
      sm3=0.
      sm4=0.
      do 9 l=1,84
      x1 = xx(l)*a11 +yy(l)*a21 + zz(l)*a31
      y1 = xx(l)*a12 +yy(l)*a22 + zz(l)*a32
      z1 = xx(l)*a13 +yy(l)*a23 + zz(l)*a33
      dx=x(l)-x1
      dy=y(l)-y1
      dz=z(l)-z1
      r1=sqrt(dx*dx+dy*dy+dz*dz)
      dx=x(l)-x1
      dy=y(l)+y1
      dz=z(l)+z1
      r2=sqrt(dx*dx+dy*dy+dz*dz)
      dx=x(l)+x1
      dy=y(l)-y1
      dz=z(l)+z1
      r3=sqrt(dx*dx+dy*dy+dz*dz)
      dx=x(l)+x1
      dy=y(l)+y1
      dz=z(l)-z1
      r4=sqrt(dx*dx+dy*dy+dz*dz)
      if(r1.gt.t3.and.r2.gt.t3.and.r3.gt.t3.and.r4.gt.t3) goto 8
      sm1=sm1+r1
      sm2=sm2+r2
      sm3=sm3+r3
      sm4=sm4+r4
    9 CONTINUE 
      sm1=sm1/84.
      sm2=sm2/84.
      sm3=sm3/84.
      sm4=sm4/84
      if(sm1.lt.T) write(2,10)i,j,k,sm1
      if(sm2.lt.T) write(2,10)i,j,k,sm2
      if(sm3.lt.T) write(2,10)i,j,k,sm3
      if(sm4.lt.T) write(2,10)i,j,k,sm4
   10 format(3i5,f8.3)   
    8 continue
    7 continue
    6 continue  
      END
      
      
 
