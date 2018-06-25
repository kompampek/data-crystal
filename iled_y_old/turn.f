      Character*1 CA,NI,OX
      BYTE A(30)
      integer*2 ih(25000),ik(25000),il(25000),nf(99)
      dimension x(399),y(399),z(399)
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
      do 5 i=1,84
      x(i)=x(i)-sx/84.
      y(i)=y(i)-sy/84.
    5 z(i)=z(i)-sz/84.
      rewind(1)
      read(1,2)CA,NI,OX
      open(2,file='turn.ent')
      WRITE(2,2)CA,NI,OX   
      rad=acos(-1.0)/180.
      CALL CPU_TIME(T0)
C
C     Spherical Polar rotations (theta,phi)
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
      do 6 k=1,84
      read(1,2)A
      x1 = x(k)*a11 +y(k)*a21 + z(k)*a31
      y1 = x(k)*a12 +y(k)*a22 + z(k)*a32
      z1 = x(k)*a13 +y(k)*a23 + z(k)*a33
C      
      x1 = x(k)*a11 +y(k)*a12 + z(k)*a13
      y1 = x(k)*a21 +y(k)*a22 + z(k)*a23
      z1 = x(k)*a31 +y(k)*a32 + z(k)*a33
    6 write(2,2)a,x1,y1,z1
      END
      
      
 
