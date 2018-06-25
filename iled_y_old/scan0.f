      byte it(5000)
      dimension ac(5000),cost(5000)
      data NUM,tip,top,tup,bip/100,4*0.0/
      open(1,file='five')
      NTRIP=0
      do 1 i=1,20
      top=0.
      bot=0.
      bs=0.
      MUM=NUM
      do 2 j=1,NUM
      read(1,*,END=5)a1,it(j),is,n1,n2,n3,v,w,ac(j),cost(j)
      ntrip=ntrip+1
      bot=bot+a1
      t=0.
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      if(it(j).eq.1) T = tanh(0.5*a1)
      top=top+a1*(cost(j)-t)**2
    2 bs=bs+ac(j)
    5 NUM=J-1
      tip=tip+top
      bip=bip+bot
      s = bot/bs
      Rmin=top/bot
      cs=0.
      do 3 j=1,num
      a1=s*ac(j)
      if(a1.lt.0.0) a1=0.
      t=0.
      a2=a1*a1
      a3=a1*a2
      IF(A1.GT.4.0) T = 1.0 - 0.5/A1 - 0.16/A2
      IF(A1.LT.1.0) T = 0.50287*A1 -.015708*A2 -0.041078*A3
      IF(T.EQ.0.0) T = 0.57745*A1 -0.13869*A2 +0.012091*A3
      if(it(j).eq.1) T = tanh(0.5*a1)
      if(a0.lt.0.0) t=-t
    3 cs=cs+a1*(cost(j)-t)**2
      tup=tup+cs
      rm = cs/bot
      bot=bot/num
      write(*,4)NUM,bot,rmin,rm
      if(NUM.NE.MUM) goto 6 
    4 format(i5,3f7.4)
      NUM=NUM+NUM
      if(NUM.gt.3200) NUM=3200
    1 CONTINUE
    6 rmin=tip/bip
      rm=tup/bip
      a1=bip/NTRIP
      write(*,4)ntrip,a1,rmin,rm
      END  
