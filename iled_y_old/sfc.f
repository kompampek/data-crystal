      dimension x(84),y(84),z(84)
      open(1,file='frag.xyz')
      do 1 i=1,84
    1 read(1,*)x(i),y(i),z(i)
      P2 = 2.*acos(-1.0)
      P4 = P2/4.
      E2=0.
      do 3 ii=1,50
      read(45,*)j,ih,ik,il
      b1=P4*(ih-ik)
      b2=P4*(ik-il)
      b3=P4*(il-ih)
      a=0
      b=0
      do 2 i=1,84
      a1=P2*ih*x(i) +b1
      a2=P2*ik*y(i) -b2
      a3=P2*il*z(i) -b3
      a=a+cos(a1)*cos(a2)*cos(a3)
    2 b=b+sin(a1)*sin(a2)*sin(a3)
      f = 2*sqrt(a*a+b*b)/sqrt(84.)
      E2=e2+f*f
    3 write(*,*)ih,ik,il,f
      E2=e2/50
      write(*,*)50,e2
      END
      
