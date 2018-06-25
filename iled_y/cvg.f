      dimension A(840),B(2200)
      integer*2 N(4,2200),ID(840),JD(840)
      data ID/840*0/
      NAX=0
      open(1,file='KK')
      do 1 i=1,2147
      read(1,*)b(i),it,is,n1,n2,n3,n4
      n2=abs(n2)
      n3=abs(n3)
      n4=abs(n4)
      n(1,i)=n1
      n(2,i)=n2
      n(3,i)=n3
      n(4,i)=n4
      NAX=MAX(NAX,N4)
      a(n1)=a(n1)+b(i)
      a(n2)=a(n2)+b(i)
      a(n3)=a(n3)+b(i)
    1 a(n4)=a(n4)+b(i)
      CLOSE(1)
      do 2 jit = 1,NAX
      amin=999.
      do 3 i=1,400
      if(a(i).lt.0.5) goto 3
      if(a(i).gt.amin) goto 3
      amin=a(i)
      next=i
    3 CONTINUE
      if(amin.eq.999) goto 6
      id(next)=1
      nqt=0
      do 4 i=1,2147
      n1=n(1,i)
      n2=n(2,i)
      n3=n(3,i)
      n4=n(4,i)
      ix=(n1-NEXT)*(n2-next)*(n3-next)*(n4-next)
      if(ix.ne.0) goto 4
      nqt=nqt+1
      a(n1)=a(n1)-b(i)
      a(n2)=a(n2)-b(i)
      a(n3)=a(n3)-b(i)
      a(n4)=a(n4)-b(i)
    4 CONTINUE
      write(2,5)next,nqt,amin
    5 format(2i5,f8.3)  
    2 CONTINUE
    6 CONTINUE
      j=0
      DO 7 i=1,NAX
      if(id(i).ne.0) goto 7
      j=j+1
      jd(j)=i
    7 ID(I)=0
      NJ=J
      ID(jd(1))=1
      ID(jd(2))=1
      ID(jd(3))=1
      write(*,*)jd(1),jd(2),jd(3)
      mj=4
      do 9 j=1,50
      it=0
      do 8 i=2147
      i1=ID(n(1,i))
      i2=id(n(2,i))
      i3=id(n(3,i))
      i4=id(n(4,i))
      ix=i1+i2+i3+i4
      if(ix.ne.3) goto 8
      ix=(1-i1)+2*(1-i2)++3*(1-i3)+4*(1-i4)
      write(*,*)ix,99
      ix=n(ix,i)
      write(*,*)ix
      read(*,*)jj
    8 CONTINUE
      if(it.ne.0) goto 9
      write(*,*)jd(mj),88
      ID(jd(mj))=1
      mj=mj+1
    9 CONTINUE   
      END
      
      
