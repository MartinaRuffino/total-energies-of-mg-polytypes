      subroutine inline(w,t,leof,ifi)
c  read line from file ifi, separate off first word
      character buf*72,w*10,x*1,y*1,t*72
      leof=0
      read(ifi,100,end=99) buf
  100 format(a72)
C|    write(6,*) 'buf=<',buf,'>'
c ------ take out comments -------
      lok=1
      do 10 i=1,72
        x=buf(i:i)
        y=x
        if(x == '{'.or.x == '}'.or.lok == 0) y=' '
        buf(i:i)=y
        if(x == '{') lok=0
  10    if(x == '}') lok=1
C|    write(6,*) 'buf=<',buf,'>'
c ------ translate < and > to ' ------------
      do 11 i=1,72
  11  if(buf(i:i) == '<'.or.buf(i:i) == '>') buf(i:i)=''''
C|    write(6,*) 'buf=<',buf,'>'
c --------- find first word -----------------
      do 1 k=1,72
      k1=k
      if(buf(k1:k1) /= ' ') goto 2
  1   continue
  2   if(k1 == 72) k1=0
      do 3 k=k1+1,72
      k2=k
      if(buf(k2:k2) == ' ') goto 4
  3   continue
  4   k2=k2-1
c ------- split into first word and rest --------------
      w=' '
      t=' '
      if(k1 <= k2) w=buf(k1:k2)
      if(k2 <= 72) t(1:k2)=' '
      if(k2 < 72) t(k2+1:72)=buf(k2+1:72)
      t=buf(k2+1:72)
C|    write(6,*) 'w=<',w,'>'
C|    write(6,*) 't=<',t,'>'
      return
  99  leof=1
      return
      end
