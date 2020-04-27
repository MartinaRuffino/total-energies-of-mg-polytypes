      subroutine mkjdos(nqp,nbmx,nfilm,nempm,nsp,nfilo,nfiup,nemlo,
     .  nemup,wgts,evl,n,w,tol,emin,emax,esciss,jdos,optmt,efermi,
     .  ndos,dos)
C- Make joint density of states from bands, sampling
C-----------------------------------------------------------------------
Ci  Input
Ci    nqp  :number of q-points
Ci    nbmx :first dimension of evl
Ci    nfilm:second dimension of optmt
Ci    nempm:third dimension of optmt
Ci    nsp  :2 for spin polarised bands, 1 otherwise
Ci    nfilo:Loop over occupied bands nfilo, nfiup
Ci    nfiup:Loop over occupied bands nfilo, nfiup
Ci    nemlo:Loop over unoccupied bands nemlo, nemup
Ci    nemup:Loop over unoccupied bands nemlo, nemup
Ci    wgts :eigenvalue weights
Ci         :NB: only abs(wtkp) is used; see bzmesh.f
Ci    evl  :eigenvalues
Ci    N, W : Methfessel-Paxton order and broadening parameters;
Ci    tol  :allowed error in DOS due to truncating the gaussian,
Ci          if negative on entry, range is set to -tol*W
Ci    emin, emax, ndos; energy range and number of energy mesh points
Ci   esciss:Shift energy of unoccupied states by esciss
Ci          (scissors operator)
Ci   jdos  :compute joint DOS, omitting matrix elements optmt
Ci   optmt :matrix elements of gradient operator
Ci   efermi:Fermi energy
Ci   ndos  :number of energy mesh points
Co  Output
Co    dos :joint density of states, or joint DOS weighted by optmt,
Co        :is ADDED to dos.
Co        :It is caller's responsibility to initialize the DOS array

Cr  Remarks
Cr    All transitions between unoccupied and occupied states are summed
Cr    assuming states below Ef are filled and states above Ef are empty.
Cr    The delta-function is broadened with the MP-sampling scheme
Cu Updates
Cu   31 Dec 10 Modified ib2 loop to conform with modified optmt
Cu   13 Sep 09 Bug fix, JDOS case with nsp=2
Cu   27 Apr 04 Take abs(wgts) in conformity with use of sign as flag
Cu             (see bzmesh 15 Sep 02 modification)
Cu   22 Feb 04 Removed initialization of dos
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nqp,nbmx,nsp,n,ndos,nfilm,nempm,nfilo,nfiup,nemlo,nemup
      double precision wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,3,nsp),
     .  optmt(3,nfilm,nempm,nsp,*),w,emin,emax,tol,wt,emesh,efermi,esciss
      logical jdos
C ... Local parameters
      integer i,isp,iq,meshpt,mesh1,mesh2,mrange,iprint,ib1,ib2,
     .  lgunit,stdo,npol,k,ibf,ibm
      double precision e,x,range,test,step,d,s,xx
      external delstp

      stdo = lgunit(1)
      npol = 3
      if (jdos) npol = 1
      step = (emax - emin) / (ndos - 1)
      if (emin < 0) call info0(10,1,0,' MKJDOS: (warning) emin<0 for joint DOS')

C --- Set tolerance and/or range of gaussians ---
      if ( tol > 0d0 ) then
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

      call info8(31,0,0,' MKJDOS:  ef=%1;6d  N.W=%1;6d   emin=%d emax=%d  %i bins'//
     .  '%N%10fRange of gaussians=%1;2d W (%i bins)  Est DOS err=%1,2;2g',
     .  efermi,n+w,emin,emax,ndos-1,range/w,2*mrange,test)

C --- Loop over spin and k-points ---
      do  isp = 1, nsp
        do  iq = 1, nqp
C ...     Double loop over occupied bands and unoccupied bands
          ibf = 0
          do  ib1 = nfilo, nfiup
            ibf = ibf+1
            if (evl(ib1,isp,iq) > efermi) cycle ! Occ state above Ef
            ibm = 0
              do  ib2 = nemlo, nemup
              ibm = ibm+1
              if (ib2 <= ib1) cycle
              if (evl(ib2,isp,iq) < efermi) cycle  ! Unocc state below Ef
              e = evl(ib2,isp,iq) - evl(ib1,isp,iq)+esciss
              if (e > emax+step+range/2) cycle  ! Outside omega window
              meshpt = (e - emin) / step
              mesh1 = max(meshpt-mrange,0)
              mesh2 = min(meshpt+mrange,ndos-1)

C ...         Loop over polarizations (if jdos=T, npol=1)
              do  k = 1, npol
                wt = abs(wgts(iq))/nsp
                if (.not. jdos) then
                  wt = wt*optmt(k,ibf,ibm,isp,iq)
                endif
                do  meshpt = mesh1, mesh2
                  emesh = emin + meshpt * step
                  x = (emesh - e) / w
                  call delstp(n,x,d,s,xx)
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
