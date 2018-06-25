      dimension x(88),y(88),z(88)\
      integer*2 nf(88)
      open(1,file='LL')
      PI2=2.*acos(-1.0)
      P4=PI2/4.
      read(*,*)ntop
      s = 30.283
      do 1 i=1,84
      read(1,*)x(i),y(i),z(i),ii,dis,nf(i)
      x(i)=x(i)*pi2
      y(i)=y(i)*pi2
      z(i)=z(i)*pi2
    1 write(*,*)x(i),y(i),z(i),nf(i)
      close(1)
      open(1,file='iled.E')
      x1=0.
      x2=0.
      y1=0.
      y2=0.
      xy=0.
      top=0.
      bot=0.
      do 2 i=1,NTOP
      read(1,*)ih,ik,il,ns,ne
      e=0.001*ne
      x1=x1+e
      x2=x2+e*e
      a=0.
      b=0.
      do 3 j=1,84
      a1=ih*x(j)+(ih-ik)*P4
      a2=ik*y(j)+(ik-il)*p4
      a3=il*z(j)+(il-ih)*P4
      a=a+nf(j)*cos(a1)*cos(a2)*cos(a3)
    3 b=b+nf(j)*sin(a1)*sin(a2)*sin(a3)
      ec = sqrt(a*a+b*b)/s
      top=top+abs(ec-e)
      bot=bot+e
      y1=y1+ec
      y2=y2+ec*ec
    2 xy=xy+e*ec
      x1=x1/ntop
      y1=y1/ntop
      ss = y1/x1
      R=top/bot
      cc = (xy/ntop -x1*y1)/sqrt((x2/ntop -x1*x1)*(y2/ntop -y1*y1))
      write(*,*)x1,y1,cc,R
      
      END
      
