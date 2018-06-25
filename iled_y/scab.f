      integer*2 ih(900),ik(900),il(900)
      dimension p(-840:840),as(900),bs(900),t(840)
      data as,bs/1800*0./
      open(1,file='iled.e')
      do 1 i=1,840
      read(1,*)ih(i),ik(i),il(i),ns,ne,a,b,c,np
      p(i)=0.001*np
      t(i)=p(i)
    1 p(-i)=-p(i)
      close(1)
      open(1,file='iled.TRIP')
      do 2 i=1,9999
      read(1,*,end=3)a,it,is,n1,n2,n3
      arg=p(n1)+p(n2)+p(n3)
      it=1
      if(is.eq.6) it=-1
      cs=it*cos(arg)
      sn=it*sin(arg)
      n2=abs(n2)
      n3=abs(n3)
      as(n1)=as(n1)+a*cs
      as(n2)=as(n2)+a*cs
      as(n3)=as(n3)+a*cs
      bs(n1)=bs(n1)+a*sn
      bs(n2)=bs(n2)+a*sn
    2 bs(n3)=bs(n3)+a*sn
    3 continue
      do 4 i=1,840
      if(i.le.10) write(*,*)as(i),bs(i)
    4 as(i)=SQRT(as(i)**2 + bs(i)**2)
      do 6 j=1,5
      read(*,*)thres
      n=0
      do 5 i=1,840
    5 if(as(i).gt.thres) n=n+1
      M=n
      f=float(n)/840.
      write(*,*)f,M
    6 continue
      read(13,*)ii
      read(13,*)(p(i),i=1,840)
      write(*,11)ii,(p(i),i=831,840)
   11 format(i5,10f7.2)   
      close(13)
      PI=acos(-1.0)
C      
      do 7 i=0,15
      i0=i
      i1=iand(i0,1)
      i0=i0/2
      i2=iand(i0,1)
      i0=i0/2
      i3=iand(i0,1)
      i0=i0/2
      i4=2*i0-1
      s=0
      do 8 j=1,840
      if(as(j).le.thres) goto 8
      ix=ih(j)*i1 +ik(j)*i2 +il(j)*i3
      ix=iand(ix,1)
      del=t(j)-i4*p(j)+ix*PI
      del=acos(cos(del))
      s=s+del
    8 continue
      del=57.296*s/M
      write(*,9)i1,i2,i3,i4,del
    9 format(4i4,f6.1)  
    7 continue    
      end  
     
      
