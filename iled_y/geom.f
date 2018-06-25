      character CH*8
      dimension tx(4),ty(4),tz(4),x(284),y(284),z(284)
      dimension xx(10),yy(10),zz(10),rr(10)
      integer*2 ix(4),iy(4),iz(4),ll(10),IND(284)
      open(1,file='cell')
      read(1,*)a,b,c
      close(1)
      nang=0
      kkmax=0
      num=0
      write(*,12)
   12 format(' XYZ file name ',$)
      read(*,13)CH
   13 format(A)
      open(98,file=CH,status='old',form='formatted')  
      do 1 i=1,80
      ind(i)=0
      read(98,*,END=7)x(i),y(i),z(i)
    1 num=num+1
    7 open(1,file='fft0.dat')
      read(1,*)nx
      read(1,*)nsym
      do 2 i=1,4
      read(1,3) tx(i),ix(i),ty(i),iy(i),tz(i),iz(i)
    2 write(*,3)tx(i),ix(i),ty(i),iy(i),tz(i),iz(i)
    3 format(f10.5,i5,10x,f10.5,5x,i5,5x,f10.5,10x,i5)
      close(1)
      nall=0
      nbad=0
      do 4 i=1,num
      kk=0
      do 5 j=1,num
      if(i.eq.j) goto 5
      do 6 k=1,nsym
      x1=x(j)*ix(k)+tx(k)
      dx=mod(abs(x(i)-x1),1.0)
      if(dx.gt.0.5) dx=1-dx
      y1=y(j)*iy(k)+ty(k)
      dy=mod(abs(y(i)-y1),1.0)
      if(dy.gt.0.5) dy=1-dy
      z1=z(j)*iz(k)+tz(k)
      dz=mod(abs(z(i)-z1),1.0)
      if(dz.gt.0.5) dz=1-dz
      r = sqrt((a*dx)**2 +(b*dy)**2 +(c*dz)**2)
      if(r.gt.2.0) goto 6
      kk=kk+1
      kkmax=max(kkmax,kk)
      xx(kk)=x1
      yy(kk)=y1
      zz(kk)=z1
      rr(kk)=r
      ll(kk)=j
    6 CONTINUE 
    5 CONTINUE
      if(kk.lt.2) goto 4
      do 8 l=1,kk-1
      do 8 m=l+1,kk
      dx=mod(abs(xx(l)-xx(m)),1.0)
      if(dx.gt.0.5) dx=1.-dx
      dy=mod(abs(yy(l)-yy(m)),1.0)
      if(dy.gt.0.5) dy=1.-dy
      dz=mod(abs(zz(l)-zz(m)),1.0)
      if(dz.gt.0.5) dz=1-dz
      r2=(a*dx)**2 +(b*dy)**2 +(c*dz)**2
      r=sqrt(r2)
      d=-(r2-rr(l)**2 -rr(m)**2)/(2.*rr(l)*rr(m))
      ang=57.296*acos(d)
      nall=nall+1
      if(ang.ge.95.0.and.ang.le.135.0) goto 14
      nbad=nbad+1
      ind(i)=ind(i)+1
      ind(ll(l))=ind(ll(l))+1
      ind(ll(m))=ind(ll(m))+1
   14 write(23,11)i,ll(l),ll(m),rr(l),rr(m),ang
   11 format(3i4,2f6.2,f6.1)   
    8 CONTINUE  
      nang=nang+1
    4 CONTINUE
      do 9 i=1,num
    9 write(24,11)i,ind(i)
      fig=float(nbad)/float(nall)
      write(*,*)num,fig
      rewind(98)
      do 10 i=1,num
      if(ind(i).gt.6) goto 10
      write(98,*)x(i),y(i),z(i)
   10 continue   
      close(98)
      END
      
