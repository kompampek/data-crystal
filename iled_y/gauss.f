      integer*2 ih,ik,il,ns,ne,np,IX
      open(unit=1,file='iled.eee',form='unformatted')
      open(unit=2,file='iled.sfc',form='unformatted')
      read(*,*)RMS,IX
      vs=0.
      amine=0.
      rad = acos(-1.0)/180.
      rms = rms*rad
      do 1 i=1,1000
      read(1)ih,ik,il,ns,ne,a,b,c,np
      CALL GAUSS(IX,RMS,V)
      np = np + 1000*v
      if(np.gt.3141) np=np-6283
      if(np.lt.-3141) np=np+6283
      write(2)ih,ik,il,ns,ne,a,b,c,np
      amine = amine + abs(v)
    1 vs=vs+v*v
      vs = vs/i
      amine = amine/(i*rad)
      vs=sqrt(vs)
      vs = vs/rad
      write(*,*)vs,amine
      end
      SUBROUTINE GAUSS(IX,S,V)
      INTEGER *2 IXX,IY,IX
      IXX=IX
      IY=IXX*899
      IF(IY)1,2,2
    1 IY=IY+32768
    2 YFL=IY
      DEL=YFL/32767.
      DEL = -ALOG(DEL)
      DEL=SQRT(DEL)
      IX=IY.AND.128
      IF(IX.EQ.0) DEL=-DEL
      IX=IY.AND.16
      IF(IX.EQ.0) DEL=-DEL
      IX=IY
      V=DEL*S
      RETURN
      END

