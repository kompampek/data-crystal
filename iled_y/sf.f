      byte NF(15)
      dimension x(85),y(85),z(85)
      COMMON/xyz/x,y,z,P4,sum,ns,ne,np,nsum,nf
      open(1,file='frag.xyz')
      PI2=2.0*acos(-1.0)
      P4=PI2/4.
      sum=0
      nsum=0
      do 1 i=1,84
      read(1,*)x(i),y(i),z(i)
      x(i)=PI2*x(i)
      y(i)=PI2*y(i)
      z(i)=PI2*z(i)
    1 write(*,*)x(i),y(i),z(i)
      close(1)
      open(1,file='P1.E')
      do 2 i=1,17904
      read(1,*)ih,ik,il,ns,ne,a,b,c,np
      CALL FCAL(ih,ik,il)
    2 CONTINUE
      e2=sum/nsum
      write(*,*)nsum,e2
      END  
       
    
      SUBROUTINE FCAL(ih,ik,il)
      byte NF(15)
      dimension x(85),y(85),z(85)
      COMMON/xyz/x,y,z,P4,sum,ns,ne,np,nsum,nf
      s=84
      sig=2.*sqrt(s)
      a=0
      b=0
      b1=(ik-ih)*p4
      b2=(il-ik)*p4
      b3=(ih-il)*p4
      do 1 i=1,84
      a1=ih*x(i)+b1
      a2=ik*y(i)+b2
      a3=il*z(i)+b3
      a=a +cos(a1)*cos(a2)*cos(a3)
    1 b=b -sin(a1)*sin(a2)*sin(a3)
      f = 4.*sqrt(a*a+b*b)/sig
      jh=abs(ih)
      jk=abs(ik)
      jl=abs(il)
      iy=jh*jk*jl
      ix=jh*jk + jl*(jk+jh)
      if(ix.eq.0) f=f/1.4142
      mp=1000*atan2(b,a)
      ne = 1000*f
      write(2,2)ih,ik,il,ns,ne,f,0.,0.,mp
    2 format(3i4,2i5,3f6.2,2i6)
      nsum=nsum+1
      sum=sum+f*f
      RETURN
      END  
