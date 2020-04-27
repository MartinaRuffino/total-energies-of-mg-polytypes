      subroutine pairp6(nbas,npadl,npadr,iax,ntab)
C- Find effective matching pairs for padded sites
      implicit none

      integer nbas,npadl,npadr,ntab(nbas+npadl+npadr+1)
      integer, parameter :: niax=10
      integer iax(niax,*)
      integer ip,it,jp,ib,npad,jb,jt,ipr,nbasp

      call getpr(ipr)
      nbasp = nbas+npadl+npadr
      npad = nbasp-nbas
      call info0(80,0,0,' pairp6: ibp   jbpp   ib   jbp  iax(6)')

C --- For each pair, find matching pairs for padded sites ---
      do  ip = nbas+1, nbasp
      do  it = ntab(ip)+1, ntab(ip+1)
        jp = iax(2,it)
        if (jp > nbasp) then
C   ...   ip derived from ib by shift as follows; see pgfpad
          if (ip <= nbas+npadl) then
            ib = ip - nbas
          else
            ib = ip - npad
          endif
C   ...   jb (a doubly padded site) came from jp by the following shift:
          jb = jp - npad
C   ...   Sanity check ... comment out because left padded and 1st PL are not identical
C          if (ntab(ip+1)-ntab(ip) /= ntab(ib+1)-ntab(ib))
C     .      call rx1('bug in pairp6: sites %s,%2i have different number of pairs',[ip,ib])

C   ...   Equivalent jp-ip pair is the jb-ib one:
          do  jt = ntab(jb)+1, ntab(jb+1)
            if (iax(2,jt) == ib .and.
     .        iax(3,it) == -iax(3,jt) .and.
     .        iax(4,it) == -iax(4,jt) .and.
     .        iax(5,it) == -iax(5,jt))  then
              iax(6,it) = jt
              call info5(80,0,0,'%,12i%,6i%,6i%,6i%,7i',iax(1,it),iax(2,it),ib,jb,iax(6,it))
              goto 5
            endif
          enddo
          call rx('bug in pairp6')
        endif
    5 enddo
      enddo

      end
