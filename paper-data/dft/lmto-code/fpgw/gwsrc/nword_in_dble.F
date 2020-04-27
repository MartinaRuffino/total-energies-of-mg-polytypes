      integer(4) function nword_in_dble()
C- Returns size of real in units of integer(4)
      real(8):: r
      integer(4):: len
c      nword_in_dble = 4 ! in alpha
c      nword_in_dble =1 ! in AIX
c      nword_in_dble=NWORD_RECORDSIZE
      inquire(iolength=len) r
      nword_in_dble = 8/len
c      write(6,*)' nword=',nword_in_dble
      end
