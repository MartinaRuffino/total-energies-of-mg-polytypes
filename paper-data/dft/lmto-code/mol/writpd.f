      subroutine writpd(time,ekin,ehtot,diff,psav,f,nbas,
     .     lrx,ldyn,lswtch)
c-  This subroutine decides what goes into the 'plot data' file
      implicit real*8 (a-h,p-z), integer (o)
      dimension psav(3,nbas),f(3,nbas),lswtch(20),vec(50)
      if (lrx == 0 .and .ldyn == 0) return
c   spare positions in lswtch say for which atoms...
      n=0
      do i=11,14
        ib=lswtch(i)
        if (ib >= 1 .and. ib <= nbas) then
           vec(n+1)=psav(1,ib)
           vec(n+2)=psav(2,ib)
           vec(n+3)=psav(3,ib)
           n=n+3
        endif
      enddo
      if (n > 0) then
         write(82,810) time,tomry(ehtot),tomry(ekin+ehtot),diff,
     .        (vec(i),i=1,n)
      else
         nn=min0(4,nbas)
         write(82,810) time,tomry(ehtot),tomry(ekin+ehtot),diff,
     .        ((psav(m,ib),m=1,3),ib=1,nn)
      endif
  810 format(f9.3, 2f10.3, f9.6, 12f10.5)
      end
      double precision function tomry(e)
c- tomry: returns mry-part only
      implicit real*8 (a-h,p-z)
      ii=dabs(e)
      e1=dabs(e)-ii
      if (e < 0d0) e1=-e1
      tomry=1000d0*e1
      end
