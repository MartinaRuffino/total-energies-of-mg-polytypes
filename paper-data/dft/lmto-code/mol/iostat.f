      subroutine iostat(nhs,evl,ewt,t,nstate,nread,ifi)
c-  i/o for states.
      implicit real*8 (a-h,p-z), integer (o)
      dimension evl(nhs),ewt(nhs),t(nhs,nhs)
      jfi=iabs(ifi)
c ---- write nstate states; all negative evl's if nstate=0 ---
      if(ifi < 0) then
      rewind jfi
      nst=nstate
      if(nstate == 0) then
        do 1 i=1,nhs
  1     if(evl(i) <= 0d0) nst=i
        endif
      write(jfi) nst
      do 10 i=1,nst
      write(jfi) evl(i),ewt(i)
  10  write(jfi) (t(j,i),j=1,nhs)
      if(iprint() >= 10) write(6,100) nst,jfi
  100 format(' iostat: write',i5,'  states on file',i3)
      endif
c ---- read states; return number in nread ---------
      if(ifi > 0) then
      rewind jfi
      read(jfi) nst
      do 20 i=1,nst
      read(jfi) evl(i),ewt(i)
  20  read(jfi) (t(j,i),j=1,nhs)
      if(iprint() >= 10) write(6,200) nst,jfi
  200 format(' iostat: read',i5,'  states from file',i3)
      nread=nst
      endif
      end
