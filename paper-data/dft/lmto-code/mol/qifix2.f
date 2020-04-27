      subroutine qifix2(alfa,n,nstate,si,z,wk,evl,wgt,dsev)
      implicit none
      integer n,nstate,i,j,iv,iprint
      double precision si(n,n),z(n,n),wk(n,n),evl(n),alfa,corr,dsev,wgt

      dsev = 0
      call dgemm('N','N',n,nstate,n,1d0,si,n,z,n,0d0,wk,n)
      iv = 1
      if (iprint() >= 50) print 345, alfa
  345 format(' qifix2:  alfa=',f10.5,/
     .  ' state  Old eval    Change      New eval')
      do  10  i = 1, nstate
        corr = 0
        do  12  j = 1, n
   12   corr = corr + z(j,i)*wk(j,i)

        if (iprint() > 10 .and. alfa*corr < 1d0 .and.
     .      alfa*corr > .3d0)
     .    print 333, i, iv, alfa*corr > 1, alfa*corr > .1d0,
     .    evl(i), corr, evl(i)+alfa*corr, evl(i)-alfa*corr
  333   format(' qifix2:',2i4,2l3,4f12.5)
        if (dabs(alfa*corr) >= 1d0) then
          iv = iv-1
        else
          dsev = dsev - wgt*alfa*corr
          if (iprint() >= 50)
     .      print 346, i, evl(i), alfa*corr, evl(i) - alfa*corr
  346     format(i4,3f12.6)
          evl(i) = evl(i) - alfa*corr
        endif
        if (dabs(alfa*corr) < 1d0 .and. i /= iv) then
          evl(iv) = evl(i)
          call dpcopy(z(1,i),z(1,iv),1,n,1d0)
        endif
        iv = iv+1
   10 continue
      end
