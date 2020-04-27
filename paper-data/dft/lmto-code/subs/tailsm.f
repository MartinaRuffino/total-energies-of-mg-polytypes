      subroutine tailsm(lrhot,nr,nrmt,nsp,a,b,rmt,rsm,nxi0,nxi,exi,rofi,
     .  rho,rhot,hfc,hfct)
C- Fit tails of rho to smoothed Hankel functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   lrhot :0 make fit for rho only; 1 make fit for rho and rhot
Ci   nr    :number of radial mesh points
Ci   nrmt  :number of points between 0..rmt
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rmt   :muffin-tin radius, in a.u.
Ci   rsm   :smoothing radius for smoothed hankel fit
Ci   nxi0  :leading dimension of hfc, nfct
Ci   nxi   :number of energies to include in fit
Ci   exi   :energies to include in fit
Ci   rofi  :radial mesh points
Ci   rho   :spherical valence charge density times 4*pi*r*r
Ci   rhot  :total charge density (used if lrhot=1)
Co Outputs
Co   hfc   :coefficients to h.f. fits to rho
Co   hfct  :coefficients to h.f. fits to rhot (if lrhot=1)
Cl Local variables
Cl   rsq   :work array holding rofi**2
Cr Remarks
Cr   A fit is make to tails of the valence charge density for r>rmt
Cr   using smoothed hankel functions of smoothing radius rsm.
Cr   Fit is constrained to agree with integrated charge.
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   19 Apr 02 Move rsq to a local array.  Altered argument list.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lrhot,nsp,nr,nrmt,nxi0,nxi
      double precision a,b,rmt,rsm
      double precision hfc(nxi0,nsp),exi(nxi),rho(nr,nsp),rofi(nr)
      double precision hfct(nxi0,nsp),rhot(nr)
C ... Local parameters
      logical lzero
      integer lx0(20),isp,ie,ir,ipr,stdo,nglob
      double precision sr,x0(0:2),xi(0:2),fpi,qout,qcst(20),qcst0(20),
     .  err,rsq(nr)
      real(8), allocatable :: xil(:),wk(:)
      integer, allocatable :: idx(:)

      fpi = 16*datan(1d0)
      call getpr(ipr)
      stdo = nglob('stdo')

C --- Tabulate smoothed Hankels on a radial mesh ---
      do  ir = 1, nr
        rsq(ir) = rofi(ir)**2
      enddo
      do  ie = 1, nxi
        lx0(ie) = 0
      enddo
      allocate(xil(nr*nxi),idx(nr*2),wk(nr*6))
C     Patch to properly handle case when rsm=0
      lzero = rsq(1) == 0 .and. rsm == 0
      if (lzero) rsq(1) = rsq(2)/100
      call hansr(rsm,0,0,nxi,lx0,exi,rsq,nr,nr,idx,wk,0,xil)
      if (lzero) rsq(1) = 0
      deallocate(idx,wk)

C --- Fit smoothed hankels to rho ---
      do  isp = 1, nsp

        if (isp == 1 .and. ipr >= 30) write (stdo,1) nxi,rmt,rsm
    1   format(/' tailsm: fit tails to',i2,' smoothed hankels, rmt=',
     .    f8.5,', rsm=',f8.5)
        if (isp == 2 .and. ipr >= 20)
     .    write(stdo,'(/'' tailsm: spin 2 ...'')')

C   ... Fitting constraints for smoothed Hankels
        do  ie = 1, nxi
          if (rsm < 1d-9) then
            sr = dsqrt(-exi(ie))*rmt
            qcst(ie) = -dsqrt(fpi)*(sr+1)*dexp(-sr)/exi(ie)
          else
            call hansmr(rmt,0d0,1/rsm,x0,1)
            call hansmr(rmt,exi(ie),1/rsm,xi,1)
            qcst(ie) = dsqrt(fpi)/exi(ie)*(-dexp(rsm**2/4*exi(ie))
     .        - rmt**3*(xi(1)-dexp(rsm**2/4*exi(ie))*x0(1)))
          endif
          qcst0(ie) = qcst(ie)
        enddo

C   ... Fit for this spin
        call hnsmft(rofi,rho(1,isp),nr,qout,a,nrmt,exi,qcst,
     .    xil,hfc(1,isp),nxi,err)

        if (isp == 1 .and. ipr >= 20 .and. ipr < 30)
     .    call awrit3('%N tailsm:  fit tails to %i functions with'
     .    //' rsm=%;7g.  rms error=%;7F',' ',80,stdo,nxi,rsm,err)

C   ... Fit a second time for the full density in rhot
        if (lrhot /= 0) then
          do  ie = 1, nxi
            qcst(ie)=qcst0(ie)
          enddo
          call hnsmft(rofi,rhot(1+(isp-1)*nr),nr,qout,a,nrmt,exi,qcst,
     .      xil,hfct(1,isp),nxi,err)
        endif

CC   ... Evaluate integral of smoothed charge density (printout)
C        qsm = 0
C        qsmr = 0
C        qsmt = 0
C        qsmrt = 0
C        rmt0 = 0
C        do  16  ie = 1, nxi
C          if (rsm < 1d-9) then
C            sr = dsqrt(-exi(ie))*rmt0
C            xx = -dsqrt(fpi)*(sr+1)*dexp(-sr)/exi(ie)
C          else
C            call hansmr(rmt0,0d0,1/rsm,x0,1)
C            call hansmr(rmt0,exi(ie),1/rsm,xi,1)
C            xx = dsqrt(fpi)/exi(ie)*(-dexp(rsm**2/4*exi(ie))
C     .        - rmt0**3*(xi(1)-dexp(rsm**2/4*exi(ie))*x0(1)))
C          endif
C          qsmr = qsmr + hfc(ie,isp)*qcst0(ie)
C          qsm  = qsm  + hfc(ie,isp)*xx
C          qsmrt = qsmrt + hfct(ie,isp)*qcst0(ie)
C          qsmt  = qsmt  + hfct(ie,isp)*xx
Cc|          write (stdo,455) ie,xx,qcst0(ie)
Cc|  455     format('  ie=',i4,'   qall',f14.8,'   qout',f14.8)
C
C   16   continue
C        if (ipr >= 30) print 345, qsm-qsmr,qsm,qsmt-qsmrt,qsmt
C  345   format (' valence:   qsm(r<rmt)=',f10.6,
C     .     '  qsm(all space)=',f10.6/
C     .     ' val+core:  qsm(r<rmt)=',f10.6,'  qsm(all space)=',f10.6)
C
C

      enddo
      deallocate(xil)

      end
