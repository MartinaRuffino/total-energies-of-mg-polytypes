      subroutine gfpdos(ib1,ib2,ndos,dosw,nl,nlo,nsp,nz,zp,offpd,gd,
     .  pdos)
C- Partial dos by Pade approximation
Ci Inputs
Ci   ndos,dosw:  number of dos, and dos window
Co Outputs
Co   pdos
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ib1,ib2,nz,nsp,nl,nlo,ndos,offpd
      double precision zp(2,nz),gd(2,nz,nsp,0:1)
      double precision dosw(2),pdos(ndos,nsp,0:1)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: cof(:)
C ... Local parameters
      logical lall
      integer isp,ib,lm,l,m,ie,offd
      integer scrwid,lgunit,iprint,ld,ksp
      double precision gpad(2),z(2),de,pi
      parameter (scrwid=80)

      if (iprint() >= 20) then
        call awrit2('%N GFPDOS: partial dos (Pade approximation)'//
     .    ' ndos=%i dosw =%2:1d',' ',scrwid,lgunit(1),ndos,dosw)
      endif
      pi = 4*datan(1d0)
      lall = nlo == nl*nl
      de = (dosw(2)-dosw(1))/(ndos-1)
      allocate(cof(nz*(nz+2)))
      do  10  isp = 1, nsp
      ksp = isp
      offd = 0
      offpd = 0
      do  20  ib = ib1, ib2
        lm = -1
        do  21  l = 0, nl-1
        do  21  m = -l, l
          lm = lm + 1
          ld = lm
          if (.not. lall) ld = l

C         call zprm('g',2,gd(1,1,isp,offpd+lm),nz,nz,1)
          call padcof(nz,zp,gd(1,1,isp,offpd+lm),nz,cof)

          do  25  ie = 1, ndos
            z(1) = dosw(1) + (ie-1)*de
            z(2) = 0
            call pade(1,z,nz,zp,nz,cof,gpad)
            pdos(ie,ksp,offd+ld) = pdos(ie,ksp,offd+ld) + gpad(2)/pi
   25     continue
C         call prmx('pdos',pdos(1,ksp,offd+ld),ndos,ndos,1)

   21   continue
        offd  = offd+ld+1
        offpd = offpd+lm+1
   20 continue
   10 continue

      deallocate(cof)
      end
