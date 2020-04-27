      subroutine nesjdos(loptic,nqp,nbmx,nfilm,nempm,nsp,nfilo,nfiup,nemlo,
     .  nemup,wgts,evl,n,w,tol,emin,emax,esciss,jdos,optmt,imref,kT,ndos,dos)
C- Make (possibly nonequilibrium) joint density of states from bands, sampling
C-----------------------------------------------------------------------
Ci  Inputs
Ci   nqp   :number of q-points
Ci   nbmx  :first dimension of evl
Ci   nfilm :second dimension of optmt
Ci   nempm :third dimension of optmt
Ci   nsp   :2 for spin polarised bands, 1 otherwise
Ci   nfilo :Loop over occupied bands nfilo, nfiup
Ci   nfiup :Loop over occupied bands nfilo, nfiup
Ci   nemlo :Loop over unoccupied bands nemlo, nemup
Ci   nemup :Loop over unoccupied bands nemlo, nemup
Ci   wgts  :eigenvalue weights
Ci         :NB: only abs(wtkp) is used; see bzmesh.f
Ci   evl   :eigenvalues
Ci   N, W  :Methfessel-Paxton order and broadening parameter
Ci   tol   :allowed error in DOS due to truncating the gaussian,
Ci         :if negative on entry, range is set to -tol*W
Ci   emin, :emax, ndos; energy range and number of energy mesh points
Ci   esciss:Shift energy of unoccupied states by esciss
Ci         :(scissors operator)
Ci  jdos   :compute joint DOS, omitting matrix elements optmt
Ci  optmt  :matrix elements of gradient operator
Ci  imref  :imref(1) = Fermi energy of holes; imref(2) = Fermi energy of electrons
Ci  ndos   :number of energy mesh points
Co  Outputs
Co   dos  :joint density of states, or joint DOS weighted by optmt,
Co        :is ADDED to dos.  Note dos is NOT initialized in this routine.
Cr  Remarks
Cr    Model nonequilibrium absorption and photoluminescence with two
Cr    quasi-Fermi levels imref(1) (holes) and imref(2) (electrons), with E_F^n > E_F^p
Cr    Absorption: valence states are filled, conduction empty => weight v->c transition probability with:
Cr      f(E-imrefp) * [1-f(E-imrefn)]
Cr    Emission: valence states are empty, conduction filled => weight c->v transition probability with:
Cr      [1-f(E-imrefp)] * f(E-imrefn)
Cr
Cr    Valence and conduction states are occupied with a probability each given by a fermi distribution
Cr    with their respective Fermi energies.  The energy-conserving Delta-function is smoothed by
Cr    a generalized gaussian (M-P integration scheme)
Cr
Cr    Possibly useful formula:  It can readily be shown that
Cr    \[f(E)[1-f(E+\Delta)] = [f(E) - f(E+\Delta)]/[1 - {e^{-\Delta/kT}}]\]
Cr    (Not needed here)
Cu Updates
Cu   17 Dec 14 Adapted from makjdos
Cu   13 Sep 09 Bug fix, JDOS case with nsp=2
Cu   27 Apr 04 Take abs(wgts) in conformity with use of sign as flag
Cu             (see bzmesh 15 Sep 02 modification)
Cu   22 Feb 04 Removed initialization of dos
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer  loptic,nqp,nbmx,nsp,n,ndos,nfilm,nempm,nfilo,nfiup,nemlo,nemup
      double precision w,emin,emax,tol,wt,emesh,imref(2),kT,esciss,wgts(nqp),
     .  evl(nbmx,nsp,nqp),dos(0:ndos-1,3,nsp),optmt(3,nfilm,nempm,nsp,*)
      logical jdos
C ... Local parameters
      integer i,isp,iq,meshpt,mesh1,mesh2,mrange,ib1,ib2,nglob,stdo,npol,k,ibf,ibm
      double precision e,x,range,test,step,d,s,xx,f,ei,ef,faco,facu
!       f(ei,ef)=1/(1+exp((ei-ef)/kT))  ! Fermi function
      f(ei,ef)=exp((ef-ei)/kT)/(1+exp((ef-ei)/kT))  ! Fermi function, avoiding fp overflow
      external delstp

      stdo = nglob('stdo')
      npol = 3
      if (jdos) npol = 1
      step = (emax - emin) / (ndos - 1)
      if (emin < 0) call info0(10,1,0,' NESJDOS: (warning) emin<0 for joint DOS')
      if (n < 0) call rx('NESJDOS: Must use M-P gaussian broadening (BZ_N>0)')

C --- Set tolerance and/or range of gaussians for smoothed delta-function ---
      if (tol > 0d0) then
        do  i = 0, ndos-1
          x = i * step / w
          call delstp(0,x,test,s,xx)
          mrange = i + 1
          if (test < tol) goto 5
        enddo
        call info0(20,1,0,' *** Warning : tol too small in MKJDOS')
    5   continue
        range = 2 * mrange * step
        test = tol
      else
        range = -tol * w
        mrange = range / ( 2 * step )
        call delstp(0,-tol/2,test,s,xx)
      endif

      call info8(31,0,0,' NESJDOS: imref=%2:1;6d  N.W=%1;6d   emin=%d emax=%d  %i bins'//
     .  '%N%10fRange of gaussians=%1;2d W (%i bins)  Est DOS err=%1,2;2g',
     .  imref,n+w,emin,emax,ndos-1,range/w,2*mrange,test)

C --- Loop over spin and k-points ---
      do  isp = 1, nsp
        do  iq = 1, nqp
C     ... Double loop over occupied bands and unoccupied bands
          ibf = 0
          do  ib1 = nfilo, nfiup
            ibf = ibf+1
!           if (evl(ib1,isp,iq) > efermi) cycle ! cutoff handled by Fermi function
            faco = f(evl(ib1,isp,iq),imref(1)) ! For absorption
            if (loptic == 9) faco = 1-faco   ! For emission
            ibm = 0
            do  ib2 = nemlo, nemup
              ibm = ibm+1
              if (ib2 < ib1) cycle
C             if (evl(ib2,isp,iq) < imrefn) cycle  ! cutoff handled by Fermi function
              facu = 1-f(evl(ib2,isp,iq),imref(2)) ! For absorption
              if (loptic == 9) facu = 1-facu     ! For emission
              e = evl(ib2,isp,iq) - evl(ib1,isp,iq)
              if (e > emax+step+range/2) cycle  ! Outside omega window
              meshpt = (e - emin) / step
              mesh1 = max(meshpt-mrange,0)
              mesh2 = min(meshpt+mrange,ndos-1)

C         ... Loop over polarizations (if jdos=T, npol=1)
              do  k = 1, npol
                wt = abs(wgts(iq))/nsp * faco*facu
                if (.not. jdos) wt = wt*optmt(k,ibf,ibm,isp,iq)
                if (wt == 0) cycle
                do  meshpt = mesh1, mesh2
                  emesh = emin + meshpt * step
                  x = (emesh - e) / w
                  call delstp(n,x,d,s,xx) ! Gaussian broadened delta-function
                  if (jdos) then
                    dos(meshpt,isp,1) = dos(meshpt,isp,1) + wt * d / w
                  else
                    dos(meshpt,k,isp) = dos(meshpt,k,isp) + wt * d / w
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      end
