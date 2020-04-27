      subroutine prrhat(nbas,s_site,s_spec,s_rhat)
C- Print info about the atomic charge densities
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl a nr rmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Cu Updates
Cu   05 Dec 12 Complete migration to f90 structures
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(*)
C ... Local parameters
      integer ib,intopt,is,lmxl,nglob,nlml,nr,nsp,stdo
      double precision a,rmt
      real(8), allocatable :: rofi(:),rwgt(:)

      nsp  = nglob('nsp')
      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')

      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        if (lmxl == -1) cycle
        allocate(rofi(nr),rwgt(nr))
        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)

        nlml = (lmxl+1)**2
        write(stdo,200) ib,rmt,nr,nlml
  200   format(/' Density at site',i3,'   rmt=',f8.4,
     .     '   nr=',i5,'   nlml=',i3)
        call prlrho('true density',nr,nlml,nsp,rofi,rwgt,
     .    s_rhat(ib)%rho1)
        call prlrho('smooth density',nr,nlml,nsp,rofi,rwgt,
     .    s_rhat(ib)%rho2)

        deallocate(rofi,rwgt)
      enddo

      end

      subroutine prlrho(str,nr,nlm,nsp,rofi,rwgt,rho)
C- Print info about charge for a single site density
      implicit none
C ... Passed parameters
      integer nr,nlm,nsp
      double precision rho(nr,nlm,nsp),rofi(nr),rwgt(nr)
      character*(*) str
C ... Local parameters
      integer stdo,nglob,ilm,l,ll,itop,ibot,i
      double precision xx,xx2,pi,srfpi,top,bot,sum,sum2

      stdo = nglob('stdo')
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)

      write(stdo,601) str
      do  ilm = 1, nlm
        l = ll(ilm)
        top = -1d10; bot = 1d10
        sum = 0d0; sum2 = 0d0
        itop = 0; ibot = 0
        do  i = 1, nr
          xx = (rho(i,ilm,1)+rho(i,ilm,nsp))/(3-nsp)
          xx2 = (rho(i,ilm,1)-rho(i,ilm,nsp))
          if (xx > top) then
            itop = i
            top = xx
          endif
          if (xx < bot) then
            ibot = i
            bot = xx
          endif
          sum = sum + rwgt(i)*xx*(rofi(i)+1d-32)**l
          sum2 = sum2 + rwgt(i)*xx2
        enddo
        xx = dmax1(dabs(top),dabs(bot))
        if (xx > 1d-6) then
          if (ilm == 1 .and. nsp == 2) then
            write(stdo,600) ilm,bot,ibot,top,itop,sum,srfpi*sum,sum2*srfpi
          elseif (ilm == 1) then
            write(stdo,600) ilm,bot,ibot,top,itop,sum,srfpi*sum
          else
            write(stdo,600) ilm,bot,ibot,top,itop,sum
          endif
        endif
  600   format(i4,f15.8,i5,f15.8,i5,f15.8,f12.6:f12.6)
  601   format(/' prlrho: ',a/'  ilm',6x,'rhomin    pnt',
     .     7x,'rhomax    pnt       moment',7x,'charge',6x,'spin')
      enddo
      end
