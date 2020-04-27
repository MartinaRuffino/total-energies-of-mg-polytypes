      subroutine bzjnos(n1,n2,n3,ep,nq,nbmx,nsp,npol,nfilo,nfiup,
     .   nemlo,nemup,emin,emax,esciss,jdos,optmt,dos,nr,ef,ntet,idtet)
C- BZ integration of joint NOS, optionally with matrix element, by linear method
C ----------------------------------------------------------------------
Ci Inputs:
Ci   n1..n3:number of divisions for the k-point mesh
Ci   ep    :energy bands
Ci   nq    :no. of irr. k-points
Ci   nbmx  :dimensions ep
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   npol  :number of polarizations; needed for dimensioning.
Ci         :npol should be 1 for jdos
Ci   nfilo :Loop over occupied bands nfilo, nfiup
Ci   nfiup :Loop over occupied bands nfilo, nfiup
Ci   nemlo :Loop over unoccupied bands nemlo, nemup
Ci   nemup :Loop over unoccupied bands nemlo, nemup
Ci   emin, emax, nr: energy window and nr of points (see remarks)
Ci   esciss:Shift energy of unoccupied states by esciss
Ci          (scissors operator)
Ci   jdos  :T: compute joint DOS, omitting matrix elements optmt
Ci   optmt :matrix elements to include with JDOS, e.g. gradient operator
Ci   dos   :density of states
Ci   nr    :number of radial mesh points
Ci   ef    :Fermi level
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of tetrahedron;
Ci         :idtet(0,i) number of tetrahedra of the i'th kind
Co Outputs:
Co  dos, Integrated Joint Density of States (idos)
Cr Remarks
Cr   Adapted from bzints to make joint density of states.
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
Cu Updates
Cu   31 Dec 10 Modified jb loop to conform with modified optmt
Cu   13 Sep 09 Bug fix, JDOS case with nsp=2
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer n1,n2,n3,nq,nbmx,nsp,idtet(0:4,*),nr,ntet,npol,
     .  nfilo,nfiup,nemlo,nemup
      double precision ep(nbmx,nsp,nq),dos(nr,npol,nsp),
     .  emin,emax,ef,esciss
      double precision optmt(3,nfilo:nfiup,nemlo:nemup,nsp,*)
      logical jdos
C Local parameters
      integer ib,jb,iq1,iq2,iq3,iq4,isp,itet,k,ibf,ibm
      double precision ec(4),ecj(4),ed(4),ebot,etop,
     .  volwgt,eboti,ebotj,etopi,etopj,wt,wt2,wt3,wtm

      call info5(10,0,0,' BZJNOS:  ef=%1;6d  emin=%1;6d;  '//
     .  'emax=%1;6d;  %i bins, %i polarizations',ef,emin,emax,nr-1,npol)
      if (emin < 0) call info0(10,0,0,
     .  ' BZJNOS: (warning) emin<0 for joint DOS')
      if (jdos) call sanrg(.true.,npol,1,1,'bzjnos','npol')
C     if (.not. jdos) call sanrg(.true.,npol,3,3,'bzjnos','npol')
      call dpzero(dos,nr*npol*nsp)
      volwgt = dble(3-nsp)/(n1*n2*n3*6)
      do  isp = 1, nsp
C --- Loop over tetrahedra ---
        do  itet = 1, ntet
        iq1 = idtet(1,itet)
        iq2 = idtet(2,itet)
        iq3 = idtet(3,itet)
        iq4 = idtet(4,itet)
        ibf = 0
          do  ib = nfilo, nfiup
          ibf = ibf+1
C --- Set up energies at 4 corners of tetrahedron ---
          ec(1) = ep(ib,isp,iq1)
          ec(2) = ep(ib,isp,iq2)
          ec(3) = ep(ib,isp,iq3)
          ec(4) = ep(ib,isp,iq4)
          etopi = dmax1(ec(1),ec(2),ec(3),ec(4))
          eboti = dmin1(ec(1),ec(2),ec(3),ec(4))
          if (eboti > ef) cycle
C ...     wt2 cludge for handling near ef
          if (dabs(etopi-eboti) > 1d-8) then
            wt2 = dmin1(1d0,(ef-eboti)/(etopi-eboti))
          else
C             print *,  ' ***WARNING***  etopi=eboti, setting wt2=1'
C             print 100,' isp,itet,ib,etopi,ef=',isp,itet,ib,etopi,ef
C  100        format(a,i3,i7,i5,2f12.6)
            wt2 = 1d0
          endif
          ibm = 0
C         jblo = max0(ib+1,nemlo)
          do  jb = nemlo, nemup
             ibm = ibm+1
             if (jb <= ib .or. jb > nbmx) cycle
C ...        Set up energies at 4 corners of tetrahedron for jb
             ecj(1) = ep(jb,isp,iq1) + esciss
             ecj(2) = ep(jb,isp,iq2) + esciss
             ecj(3) = ep(jb,isp,iq3) + esciss
             ecj(4) = ep(jb,isp,iq4) + esciss
             etopj = dmax1(ecj(1),ecj(2),ecj(3),ecj(4))
             ebotj = dmin1(ecj(1),ecj(2),ecj(3),ecj(4))
             if ( etopj < ef+emin) cycle
C             print 335, ib,ec
C             print 335, jb,ecj
C  335        format(i4,4f12.6)
             ed(1) = ep(jb,isp,iq1)-ep(ib,isp,iq1)
             ed(2) = ep(jb,isp,iq2)-ep(ib,isp,iq2)
             ed(3) = ep(jb,isp,iq3)-ep(ib,isp,iq3)
             ed(4) = ep(jb,isp,iq4)-ep(ib,isp,iq4)
             etop = dmax1(ed(1),ed(2),ed(3),ed(4))
             ebot = dmin1(ed(1),ed(2),ed(3),ed(4))
             if (ebot > emax) cycle
C ...        wt3 cludge for handling near ef+emin
             if (dabs(etopj-ebotj) > 1d-8) then
               wt3 = dmin1(1d0,(etopj-ef-emin)/(etopj-ebotj))
             else
C               print *,  ' ***WARNING***  etopj=ebotj, setting wt3=1'
C               print 100,' isp,itet,jb,etopj,ef+emin=',
C     .           isp,itet,jb,etopj,ef+emin
               wt3 = 1d0
             endif
C             print 336, 'emax,ebot,etop=      ',emax,ebot,etop
C             print 336, 'eboti,etopi, ef=     ',eboti,etopi,ef
C             print 336, 'ebotj,etopj, ef+emin=',ebotj,etopj,ef+emin
C  336        format(a,3f12.6)
C             if ( ebot < emax  .and. (wt2+wt3 /= 2)) then
C               print 336, 'wt,wt2,wt3=',wt,wt2,wt3
C             endif
             if (wt2 > 1 .or. wt3 > 1) call rx('bug in bzjnos')
             do  k = 1, npol
               wt = volwgt*idtet(0,itet)
               if (jdos) then
                 call slinz(wt*wt2*wt3,ed,emin,emax,dos(1,1,isp),nr)
               else
                 wtm = optmt(k,ib,jb,isp,iq1)+optmt(k,ib,jb,isp,iq2)
     .               + optmt(k,ib,jb,isp,iq3)+optmt(k,ib,jb,isp,iq4)
                 wt = wt*wtm / 4d0
                 call slinz(wt*wt2*wt3,ed,emin,emax,dos(1,k,isp),nr)
               endif
             enddo
           enddo
         enddo
       enddo
      enddo

      end
