      subroutine rdfa(spid,z,r,konf,q,lmax,nxi,lxi,exi,
     .   nspec,nspmx,n0,ifi)
c  read free-atom data from control file
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(1),z(1),konf(n0,1),q(n0,1),lmax(1),nxi(1),
     .   lxi(n0,1),exi(n0,1)
      character*8 spid(1),sid,sl,dup
      character t*72,w*10,wl*10,ttt*20
      character dum*8
      nt=72
      nspec=0
      wl=' '
      sid='    '
      sl='    '
      dup='"   '
      rewind ifi
  90  call inline(w,t,leof,ifi)
      if(leof == 1 .or. w == 'exit') return
      if(w == '"') w=wl
c ------- spec ---------------
      if(w == 'spec') then
        nspec=nspec+1
        if(nspec > nspmx) call rx('readmc: too many species''')
        read(t,*,end=9) spid(nspec),z(nspec)
        sl=spid(nspec)
        endif
c ------- rad ----------------
      if(w == 'rad') then
        read(t,*) sid
C|        write(6,*) w,'  sid=',sid,'   sl=',sl,'    dup=',dup
        if(sid == dup) sid=sl
C|        write(6,*) w,'  sid=',sid,'   sl=',sl,'    dup=',dup
        call spcfnd(sid,spid,nspec,j)
        if(j == 0) goto 95
        read(t,*,end=9) dum,r(j)
        endif
c ------- chd ----------------
      if(w == 'chd') then
        read(t,*,end=9) sid,ii,ttt
C|        write(6,*) w,'  sid=',sid,'   sl=',sl,'    dup=',dup
        if(sid == dup) sid=sl
C|        write(6,*) w,'  sid=',sid,'   sl=',sl,'    dup=',dup
        call spcfnd(sid,spid,nspec,j)
        if(j == 0) goto 95
        call iunpac(ttt,lxi(1,j),nxi(j))
        read(t,*,end=9) sid,ii,ttt,(exi(i,j),i=1,nxi(j))
        endif
c ------- fa -----------------
      if(w == 'fa') then
        read(t,*,end=9) sid,ttt
C|        write(6,*) w,'  sid=',sid,'   sl=',sl,'    dup=',dup
        if(sid == dup) sid=sl
C|        write(6,*) w,'  sid=',sid,'   sl=',sl,'    dup=',dup
        call spcfnd(sid,spid,nspec,j)
        if(j == 0) goto 95
        call iunpac(ttt,konf(1,j),nl)
        lmax(j)=nl-1
        read(t,*,end=9) sid,ttt,(q(i,j),i=1,nl)
        endif
c -----------------------------------------
      wl=w
      goto 90
  95  write(6,*) '----  ',w,t
      call rx('rdfa: species not yet defined''')
  9   call rx('rdfa: eol reached, 1st word='//w)
      end
