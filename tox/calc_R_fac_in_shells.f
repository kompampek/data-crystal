C      calc_R_fac_in_shells
       dimension FO1(30998),FN1(30998),FN(30998)
       OPEN(1,file='s_sort_output')
       do 2 k=1,30998
       read(1,9,end=2)ih,ik,il,NS,NE,FO1(k),FN1(k),FN(k),PN1,PN    !!YV 
       write(*,9)ih,ik,il,NS,NE,FO1(k),FN1(k),FN(k),PN1,PN    !!YV
    2 continue      
      close(1)
    9 format(3i4,2x,i5,2x,i4,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2)  
      iflag=1
      R_fac_fo_fc_ccp4=0
      R_fac_fo_fc=0
      R_fac_fc_ccp4_fc=0
      R_fo=0
      do 3 i=1,30800
      if (i.lt.200*iflag) then
      R_fac_fo_fc_ccp4=R_fac_fo_fc_ccp4+abs(FO1(i)-FN1(i))
      R_fac_fo_fc=R_fac_fo_fc+abs(FO1(i)-FN(i))
      R_fac_fc_ccp4_fc= R_fac_fc_ccp4_fc+abs(FN1(i)-FN(i))
      R_fo=R_fo+FO1(i)
      else
      R_fac_fo_fc_ccp4=R_fac_fo_fc_ccp4/R_Fo
      R_fac_fo_fc=R_fac_fo_fc/R_Fo
      R_fac_fc_ccp4_fc=R_fac_fc_ccp4_fc/R_Fo
      write (3,10) R_fac_fo_fc_ccp4,R_fac_fo_fc,R_fac_fc_ccp4_fc
      R_fo=0
      R_fac_fo_fc_ccp4=0
      R_fac_fo_fc=0
      R_fac_fc_ccp4_fc=0 
      iflag=iflag+1
      end if
    3 continue  
   10 format (F10.4,2x,F10.4,2x,F10.4)
       end