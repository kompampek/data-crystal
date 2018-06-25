      dimension n(80),m(80)
      do 1 i=1,80
      n(i)=0
    1 m(i)=0
      do 2 i=1,9999
      read(23,*,end=4)j,k,l,d1,d2,ang
      if(ang.gt.145.0.or.ang.lt.95.) goto 3
      n(j)=n(j)+1
      n(k)=n(k)+1
      n(l)=n(l)+1
      goto 2
    3 m(j)=m(j)+1
      m(k)=m(k)+1
      m(l)=m(l)+1
    2 continue
    4 continue
      do 5 i=1,80
    5 write(98,6)i,n(i),m(i)
    6 format(3i5)
      END  
          
      
