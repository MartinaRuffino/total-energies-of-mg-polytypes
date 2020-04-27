      subroutine ioatom(atid,pos,pol,nbas,nbx,ifi)
c  read in atomic positions from control file
      implicit real*8 (a-h,p-z), integer (o)
      dimension pos(3,1),pol(3,1)
      character*8 atid(1)
      character t*73,w*10,wl*10,ttt*20
      nbas=0
      wl=' '
      rewind ifi
  90  call inline(w,t,leof,ifi)
      t(73:73)='/'
      if(leof == 1 .or. w == 'exit') return
      if(w == '"') w=wl
c ------- pos ----------------
      if(w == 'atom') then
        nbas=nbas+1
        if(nbas > nbx) call rx('ioatom: too many atoms')
        do 1 m=1,3
        pos(m,nbas)=0d0
  1     pol(m,nbas)=0d0
        read(t,*,err=99) atid(nbas),
     .     (pos(m,nbas),m=1,3),(pol(m,nbas),m=1,3)
        endif
c -----------------------------------------
      wl=w
      goto 90
   99 call rx('ioatom: eol reached, 1st word='//w)

      end
