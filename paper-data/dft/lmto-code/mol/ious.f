      subroutine ious(is,oul,osl,a,rmt,nr,lmxa,ifi)
C- Assembles u,s for a one species by averaging over correponding basis
      implicit none
      integer is,ifi,oul,osl,lmxa
      integer ouj,osj,ib,np,nr,nbas,ips,nis,l
      double precision rmt,a,pnu(10),v0
      real w(1)
      common /w/ w

      rewind (ifi)
      read(ifi) nbas

      nis = 0
      do  10  ib = 1, nbas
        read(ifi) nr,a,rmt,lmxa,np,ips
        read(ifi) v0, (pnu(l+1), l=0,lmxa)
C       print *, v0, (pnu(l+1), l=0,lmxa)
        if (ib == 1) then
          call defrr(oul,     np)
          call dpzero(w(oul), np)
          call defrr(osl, np)
          call dpzero(w(osl), np)
          call defrr(ouj,     np)
          call defrr(osj,     np)
        endif
        call dpdump(w(ouj),np,ifi)
        call dpdump(w(osj),np,ifi)
        if (ips == is) then
          nis = nis+1
          call dpadd(w(oul),w(ouj),1,np,1d0)
          call dpadd(w(osl),w(osj),1,np,1d0)
        endif
   10 continue
      call rlse(ouj)
      call dscal(np,1/dble(nis),w(oul),1)
      call dscal(np,1/dble(nis),w(osl),1)
      end
