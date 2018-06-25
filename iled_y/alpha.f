      parameter (mp=1200,mt=12000)
      integer*2 n1(mt),n2(mt),n3(mt),nt(mt),ih(mp),ik(mp),il(mp)
      dimension p(-mp:mp),ahkl(mt),a(mp),b(mp),c(mp),q(mp),cs(mt)
      dimension r(0:7,mp)
      p6=acos(-1.0)/6.
      p4=acos(-1.0)/4.
      pmin=99.
      NREF=840
      do 14 i=1,mp
      a(i)=0.
      b(i)=0.
   14 c(i)=0.
      open(1,file='F13')
      do 1 i=1,99
      read(1,*,end=2)nam,ncyc,rmin,cosav,pend,pset
      read(1,3) (q(j),j=1,NREF)
    3 format(20f6.2)
      if(pset.gt.pmin) goto 1
      pmin=pset
      mam=nam
      do 4 j=1,NREF
      p(j)=q(j)
    4 p(-j)=-q(j)
    1 continue
    2 continue
      write(*,*)mam,pmin
      close(1)
      open(1,file='iled.T')
      do 5 i=1,10000
    5 read(1,*,end=6)ahkl(i),it,is,n1(i),n2(i),n3(i),nt(i),cs(i)
    6 ntrip=i-1
      close(1)
      open(1,file='iled.E')
      do 9 i=1,NREF
      read(1,*)ih(i),ik(i),il(i),ns,ne,ax,bx,cx,np
    9 q(i)=0.001*np
      close(1)
      pmin=99.
      do 10 i=0,15
      i1=iand(i,1)
      i0=i/2
      i2=iand(i0,1)
      i0=i0/2
      i3=iand(i0,1)
      i4=i0/2
      i4=i4+i4-1
      s=0.
      do 11 j=1,NREF
      ix=abs(ih(j)*i1+ik(j)*i2+il(j)*i3)
      ix=iand(ix,1)
      del=abs(q(j)-i4*p(j)+ix*3.141592)
      del=acos(cos(del))
   11 s=s+del
      del=57.296*s/NREF
      write(*,*)i1,i2,i3,i4,del
      if(del.gt.pmin) goto 10
      j1=i1
      j2=i2
      j3=i3
      j4=i4
      pmin=del
   10 CONTINUE
      write(*,*)ntrip,j1,j2,j3,j4
      do 12 i=1,NREF
      ix=abs(ih(i)*j1+ik(i)*j2+il(i)*j3)
      ix=iand(ix,1)
      pp=j4*p(i) +ix*3.141592
      p(i)=mod(pp,6.283185)
      q(i)=57.296*acos(cos(p(i)-q(i)))
   12 p(-i)=-p(i)
      do 19 i=1,8
      do 19 j=1,nref
   19 r(i,j)=0.   
      do 7 i=1,ntrip
      m1=n1(i)
      m2=abs(n2(i))
      m3=abs(n3(i))
      arg=p(m1)+p(n2(i))+p(n3(i))+nt(i)*P6
      ah=ahkl(i)
      a2=ah*ah
      a3=aH*a2
      t=0.
      if(ah.gt.4.0) T = 1.0 -0.5/ah - 0.16/a2
      if(ah.lt.1.0) T =  0.50287*ah -0.15708*a2 -0.041078*a3
      if(t.eq.0.0) T = 0.57745*ah -0.13869*a2 -0.012091*a3
      ra=ah*cos(arg)
      rb=ah*sin(arg)
      a(m1)=a(m1)+ra
      a(m2)=a(m2)+ra
      a(m3)=a(m3)+ra
      c(m1)=c(m1)+ah
      c(m2)=c(m2)+ah
      c(m3)=c(m3)+ah
      b(m1)=b(m1)+rb
      do 15 k=0,7
   15 r(k,m1)=r(k,m1) +ah*(T-cos(arg+k*P4))**2
      rbb=rb
      ar=arg
      if(n2(i).lt.0) rbb=-rb
      if(n2(i).lt.0) ar=-ar
      b(m2)=b(m2)+rbb
      do 16 k=0,7
   16 r(k,m2)=r(k,m2) +ah*(T-cos(ar+k*p4))**2   
      rbb=rb
      ar=arg
      if(n3(i).lt.0) rbb=-rb
      if(n3(i).lt.0) ar=-ar
      do 17 k=0,7
   17 r(k,m3)=r(k,m3) +ah*(T-cos(ar+k*P4))**2   
      b(m3)=b(m3)+rbb
    7 CONTINUE
      open(1,file='ALFA')
      do 8 i=1,NREF
      alpha=sqrt(a(i)*a(i)+b(i)*b(i))
      a1=0.
      a2=0.
      do 18 j=0,7
      a1=a1 +r(j,i)
   18 a2=a2 +r(j,i)**2
      a2=a2/8.
      a1=a1/8.
      sig=sqrt(a2-a1*a1)  
    8 write(1,13)i,alpha,q(i),sig
   13 format(i5,3f8.2)
      END
         
      
