      real*8 x,y,z,rand
      read(*,*)iseed
      call srand(iseed)
      do 1 i=1,1000
      x=rand()
    1 write(2,2)x
    2 format(f8.4)
      END
