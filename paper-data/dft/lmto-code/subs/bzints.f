      subroutine bzints(n1,n2,n3,ep,wp,nq,nband,nbmx,nsp,emin,emax,
     .  dos,nr,ef,job,ntet,idtet,sumev,sumwp)
C- BZ integration by linear tetrahedron method
C ----------------------------------------------------------------------
Ci Inputs
Ci   nq    :no. of irr. k-points
Ci   ep    :energy bands
Ci   nband :number of bands
Ci   nbmx  :dimensions ep
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   emin  :(job=1) lower bound for energy
Ci   emax  :(job=1) upper bound for energy
Ci   dos   :(job=2) If nr=1, dos(1) used for printout: written as DOS(EF)
Ci   nr    :(job=1) number of points for idos
Ci         :(job=2) set nr=1 to print out DOS(EF); dos not used.
Ci   job   :(job=1) makes idos;
Ci         :(job=2) makes Bloechl weights.
Ci         :sign of job<0:  spin up, spin down bands are coupled.
Ci         :In this case, nsp must be unity!
Ci   ntet  :No. of different tetrahedra
Ci   idtet :idtet(0,i) = number of tetrahedra of the i'th kind
Ci         :idtet(1-4,i) points to the 4 irr. k-points defining tetr.
Cio Inputs/Outputs
Cio  ef    :Fermi energy
Cio        :job=1: output
Cio        :job=2: input
Co Outputs
Co   sumev :sum of eigenvalues (job = 2)
Co   dos   :Integrated density of states (idos) (job = 1)
Co   wp    :Bloechl quadrature weights (job = 2)
Cu Updates
Cu   17 Jan 05 Returns sumwp
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer n1,n2,n3,nq,nband,nbmx,nsp,idtet(0:4,*),nr,job,ntet
      double precision ep(nbmx,nsp,nq),dos(nr,nsp),wp(nband,nsp,nq),
     .  emin,emax,ef,sumev,sumwp
C Local parameters
      integer ib,iq,iq1,iq2,iq3,iq4,isp,itet,jjob
      integer lgunit,ipr,stdo
      double precision ec(4),wc(4,2),ebot,etop,sev1,sev2,sumwm,
     .  volwgt

      call getpr(ipr)
      stdo = lgunit(1)
C     stdl = lgunit(2)
      jjob = iabs(job)
      if (job < 0 .and. nsp == 2 .or.
     .   jjob /= 1 .and. jjob /= 2) call fexit2(-1,111,
     .  ' Exit -1 BZINTS: job=%i and nsp=%i not allowed',job,nsp)
      if (jjob == 1) call dpzero(dos,nsp*nr)
      if (jjob == 2) call dpzero(wp,nband*nsp*nq)
      sev1 = 0d0
      sev2 = 0d0
      volwgt = dble(3-nsp)/(n1*n2*n3*6)
      if (job < 0) volwgt = volwgt/2
      do  isp = 1, nsp
C --- Loop over tetrahedra ---
        do  itet = 1, ntet
         iq1 = idtet(1,itet)
         iq2 = idtet(2,itet)
         iq3 = idtet(3,itet)
         iq4 = idtet(4,itet)
          do  ib = 1, nband
C --- Set up energies at 4 corners of tetrahedron ---
           ec(1) = ep(ib,isp,iq1)
           ec(2) = ep(ib,isp,iq2)
           ec(3) = ep(ib,isp,iq3)
           ec(4) = ep(ib,isp,iq4)
           etop = dmax1(ec(1),ec(2),ec(3),ec(4))
           ebot = dmin1(ec(1),ec(2),ec(3),ec(4))
           if (jjob == 1) then
             if ( ebot < emax )
     .       call slinz(volwgt*idtet(0,itet),ec,emin,emax,dos(1,isp),nr)
           elseif ( ef >= ebot ) then
               call fswgts(volwgt*idtet(0,itet),ec,ef,etop,wc)
               sev1 = sev1 + wc(1,1)*ec(1) + wc(2,1)*ec(2) +
     .                       wc(3,1)*ec(3) + wc(4,1)*ec(4)
               sev2 = sev2 + wc(1,2)*ec(1) + wc(2,2)*ec(2) +
     .                       wc(3,2)*ec(3) + wc(4,2)*ec(4)
               wp(ib,isp,iq1) = wp(ib,isp,iq1) + wc(1,1) + wc(1,2)
               wp(ib,isp,iq2) = wp(ib,isp,iq2) + wc(2,1) + wc(2,2)
               wp(ib,isp,iq3) = wp(ib,isp,iq3) + wc(3,1) + wc(3,2)
               wp(ib,isp,iq4) = wp(ib,isp,iq4) + wc(4,1) + wc(4,2)
             endif
          enddo
        enddo
      enddo
      if (jjob == 2) then
        sumev = sev1 + sev2
        sumwp = 0d0
        do  isp = 1, nsp
          sumwm = 0d0
          do  ib = 1, nband
            do  iq = 1, nq
              sumwm = sumwm + wp(ib,isp,iq)
              sumwp = sumwp + wp(ib,isp,iq)
            enddo
          enddo
        enddo
        if ( ipr >= 10 ) then
C     ... When bands are coupled, moment from bands makes no sense
          if (nsp == 1 .and. nr == 1) then
            write(stdo,1) ef,sumwp,dos(1,1),sev1+sev2,sev2
          elseif (nsp == 1) then
            write(stdo,2) ef,sumwp,sev1+sev2,sev2
C           write(stdl,2) ef,sumwp,sev1+sev2,sev2
          elseif (nsp == 2 .and. job > 0) then
            sumwm = sumwp-2*sumwm
            write(stdo,3) ef,sumwp,sumwm,sev1+sev2,sev2
          endif
          if (dabs(sumwp-nint(sumwp)) > 1d-6 .and. ipr >= 20) then
            call logwarn(10,' warning! non-integral number of electrons ---'//
     .        ' possible band crossing at E_f')
          endif
        endif
C      if (nsp == 1) call awrit4('bzi tetra ef %,6;6d  q %,6;6d  '//
C     .    'sev %,6;6d  Blo %,6;6d',' ',80,stdl,ef,sumwp,sev1+sev2,sev2)
C      if (nsp == 2) call awrit5('bzi tetra ef %,6;6d  q %,6;6d  mom'//
C     .  ' %,6;6d  sev %,6;6d  Blo %,6;6d',' ',80,stdl,ef,sumwp,sumwm,
C     .  sev1+sev2,sev2)
      endif

    1 format(1x,'BZINTS: Fermi energy:',f14.6,';',F11.6,
     .  ' electrons;  D(Ef):',F9.3/
     .       9x,'Sum occ. bands:',f13.7,
     .       '  incl. Bloechl correction:',f12.6)
    2 format(1x,'BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons'/
     .       9x,'Sum occ. bands:',f12.6,
     .       ' incl. Bloechl correction:',f10.6)
    3 format(' BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons',
     .  '  amom=',f8.4/9x,'Sum occ. bands:',f13.7,
     .       ',  incl. Bloechl correction:',f12.6)
    4 format(' (warning): non-integral number of electrons ---',
     .        ' possible band crossing at E_f')
      end

