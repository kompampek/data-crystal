      integer*2 ih(900),ik(900),il(900)
      dimension x(84),y(84),z(84),p(900)
      open(1,file='FRAG.xyz')
      do 1 i=1,84
    1 read(1,*)x(i),y(i),z(i)
      close(1)
      open(1,file='iled.E')
      do 2 i=1,840
      read(1,*)ih(i),ik(i),il(i),ns,ne,a,b,c,np
    2 p(i)=0.001*np
      close(1)
      PI2=2.*acos(-1.0)
      P2=PI2/4.
      del=0.
      do 3 i=1,840
      b1=P2*(ih(i)-ik(i))
      b2=P2*(ik(i)-il(i))
      b3=P2*(il(i)-ih(i))
      a=0
      b=0
      do 4 j=1,84
      a1=PI2*ih(i)*x(j)-b1
      a2=PI2*ik(i)*y(j)-b2
      a3=PI2*il(i)*z(j)-b3
      a=a+cos(a1)*cos(a2)*cos(a3)
    4 b=b+sin(a1)*sin(a2)*sin(a3)
      ph=atan2(-b,a)
      if(i.lt.10) write(*,*)ph,p(i)
    3 del = del + acos(cos(p(i)-ph))
      del = 57.296*del/840.
      write(*,*)del
      END
      
      
      
