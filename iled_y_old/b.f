      integer N(257)
      dimension cc(257),ang(257)
      open(1,file='ppt')
      read(1,*)nn
      write(*,*)nn
      kk=0
      do 1 i=1,257
      read(1,*,end=4)n(i),cc(i),an,ang(i)
    1 write(*,*)n(i),cc(i),an,ang(i)
    4 LAST=i
      lost=last/2
      write(*,*)last,lost
      do 2 i=1,lost
      j=last-i
      jj = n(i)+n(j)
      if(jj.ne.kk) write(*,*)jj
      kk=jj
      d=cc(i)-cc(j)
      if(d.ge.0.0) write(2,3)n(i),n(j),cc(i),d,ang(i)
    3 format(2i8,2f8.5,f5.1,i4)
      if(d.lt.0.0) write(2,3)n(j),n(i),cc(j),-d,ang(j),-1
    2 continue
      end
        
