      subroutine trysop(ib,jb,g,nbas,ipclas,tau,rb,qb,lok,trans)
C- Tests whether crystal as seen from site ib same as that seen from jb
C  after point operation g.  If so, returns trans and lok=1.
      implicit none
      integer ipclas(*)
      integer lok,ib,jb,nbas
      double precision rb(9),qb(9),g(3,3),tau(3,1),d(3),dif(3),trans(3)
      integer kb,m,k,kbp,kk,iprint,i1mach,iclas1,iclas2

      lok = 0
      do  kb = 1, nbas
        do  m = 1, 3
          d(m) = 0
          do  k = 1, 3
            d(m) = d(m) + g(m,k)*(tau(k,kb)-tau(k,ib))
          enddo
        enddo
        kbp = 0
        do  kk = 1, nbas
          do  m = 1, 3
            dif(m) = d(m) - (tau(m,kk)-tau(m,jb))
          enddo
          call shorbz(dif,dif,rb,qb)
          if (dif(1)**2+dif(2)**2+dif(3)**2 < 1d-10) then
            kbp = kk
            goto 5
          endif
        enddo
        if (iprint() >= 60) then
          call awrit3(' trysop ib=%i:  no site analogous to'//
     .      ' site %i from site %i',' ',80,i1mach(2),ib,kb,jb)
          if (iprint() >= 80) then
         call awrit1(' missing atom at: %3:1,6;6d',' ',80,i1mach(2),dif)
          do  m = 1, 3
            dif(m) = tau(m,kb)-tau(m,ib)
          enddo
         call awrit1('    vector kb-ib: %3:1,6;6d',' ',80,i1mach(2),dif)
         call awrit1('  after rotation: %3:1,6;6d',' ',80,i1mach(2),d)
         endif
       endif
       return
    5  continue
       iclas1 = ipclas(kb)
       iclas2 = ipclas(kbp)
       if (iprint() >= 60.and.iclas1 /= iclas2 .or. iprint() >= 70) then
         call awrit6(' trysop ib=%i jb=%i:  site %i class %i'//
     .     ' maps to site %i class %i',' ',80,i1mach(2),ib,jb,
     .     kb,iclas1,kbp,iclas2)
         if (iprint() >= 80) then
         do  m = 1, 3
           dif(m) = tau(m,kb)-tau(m,ib)
         enddo
         call awrit1('    vector kb-ib: %3:1,6;6d',' ',80,i1mach(2),dif)
         call awrit1('  after rotation: %3:1,6;6d',' ',80,i1mach(2),d)
         endif
       endif

       if (iclas1 /= iclas2) return
      enddo

      lok = 1
      do  m = 1, 3
      trans(m) = tau(m,jb)
        do  k = 1, 3
          trans(m) = trans(m)-g(m,k)*tau(k,ib)
        enddo
      enddo
      call shorbz(trans,trans,rb,qb)
      end
