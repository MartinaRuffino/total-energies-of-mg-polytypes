      subroutine pgpdos(ndos,dosw,pgplp,glist,npl,nl,nlo,nsp,nspc,nz,zp,
     .  gd,pdos)
C- Partial dos by Pade approximation
Ci Inputs
Ci   ndos,dosw:  number of dos, and dos window
Co Outputs
Co   pdos
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nz,nsp,nspc,npl,nl,nlo,ndos
      integer glist(-1:*),pgplp(6,-1:npl)
      double precision zp(2,nz),gd(2,nz,nsp,nspc,0:*)
      double precision dosw(2),pdos(ndos,nsp,0:*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: cof(:)
C ... Local parameters
      logical lall
      integer isp,ib,lm,l,m,ib1,ib2,ipl,kpl,ie,offpd,offd
      integer scrwid,lgunit,iprint,ld,jsp
      double precision gpad(2),z(2),de,pi,wt
      parameter (scrwid=80)
      if (iprint() >= 20) then
        call awrit2('%N PGPDOS: partial dos (Pade approximation)'//
     .    ' ndos=%i dosw =%2:1d',' ',scrwid,lgunit(1),ndos,dosw)
      endif
      pi = 4*datan(1d0)
      wt = -(2d0/nsp)/pi
      lall = nlo == nl*nl
      de = (dosw(2)-dosw(1))/(ndos-1)
      allocate(cof(nz*(nz+2)))
C --- For each spin, do ---
      do  10  isp = 1, nsp
C     jsp = 1 in collinear case, isp in the noncollinear; thus
C     (isp,jsp) are indices to the diagonal in the spinor part of gd
      jsp = min(isp,nspc)
      offd = 0
      offpd = 0
      do  10  kpl = 1, glist(-1)
      ipl = glist(kpl-1)
      call gtibpl(ipl,npl,pgplp,ib1,ib2)
      do  20  ib = ib1, ib2
        lm = -1
        do  21  l = 0, nl-1
        do  21  m = -l, l
          lm = lm + 1
          ld = lm
          if (.not. lall) ld = l

C         call zprm('g',2,gd(1,1,isp,jsp,offpd+lm),nz,nz,1)
          call padcof(nz,zp,gd(1,1,isp,jsp,offpd+lm),nz,cof)

          do  25  ie = 1, ndos
            z(1) = dosw(1) + (ie-1)*de
            z(2) = 0
            call pade(1,z,nz,zp,nz,cof,gpad)
            pdos(ie,isp,offd+ld) = pdos(ie,isp,offd+ld) + wt*gpad(2)
   25     continue
C         call prmx('pdos',pdos(1,isp,offd+ld),ndos,ndos,1)

   21   continue
        offd  = offd+ld+1
        offpd = offpd+lm+1
   20 continue
   10 continue

      deallocate(cof)
      end
