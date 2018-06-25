      open(1,file='iled.e')
      do 1 i=1,2000
      read(1,*)ih,ik,il,ns,ne,a,b,c,np
      IX=IH*IK*IL
      IF(IX.NE.0) GOTO 1
      e=0.001*NE
      P=0.001*NP
      write(2,2)ih,ik,il,e,p
    2 FORMAT(3I4,3f7.2,i6)
    1 CONTINUE  
      END
      
