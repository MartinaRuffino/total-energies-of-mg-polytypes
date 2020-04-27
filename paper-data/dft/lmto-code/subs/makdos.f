      subroutine makdos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax,
     .                  ndos,dos)
C- Make density of states from bands
C-----------------------------------------------------------------------
Ci  Input
Ci   nqp   :number of q-points
Ci   nband :number of bands
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   wgts  :band weights
Ci   evl   :band eigenvalues
Ci   n,w   :Methfessel-Paxton order and broadening parameters
Ci   tol   :(tol>0) allowed error in DOS due to truncating the gaussian,
Ci         :        to a finite energy range (number of bins)
Ci         :(tol<0) dimensionless energy window specifying truncation
Ci         :        of gaussian.  Energy window for which gaussian is
Ci         :        taken to be nonzero is set to -tol*w
Ci   emin, emax, ndos: energy range and number of energy mesh points
Ci   nbmx  :leading dimension of evl
Co  Ouput
Co    dos: density of states
C-----------------------------------------------------------------------
      implicit none
      integer nqp,nband,nbmx,nsp,n,ndos
      double precision wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp),
     .                 w,emin,emax,tol,wt,emesh
      integer i,isp,iband,iq,meshpt,mesh1,mesh2,mrange,iprint
      double precision e,x,range,test,step,d,s,xx
      external delstp

      call dpzero(dos,nsp*ndos)
      step = (emax - emin) / (ndos - 1)
      if ( tol > 0d0 ) then
        do  i = 0, ndos-1
          x = i * step / w
          call delstp(0,x,test,s,xx)
          if ( test < tol ) then
            mrange = i + 1
            goto 5
          endif
        enddo
        if (iprint() > 30) print *,'makdos (warning) : tol too small'
    5   continue
        range = 2 * mrange * step
        test = tol
      else
        range = -tol * w
        mrange = range / ( 2 * step )
        call delstp(0,-tol/2,test,s,xx)
      endif
      if (iprint() > 30) write (*,100) range/w,2*mrange,test
  100 format(/1x,'MAKDOS :  range of gaussians is ',f5.2,
     .       'W (',i4,' bins).'
     .       /11x,'Error estimate in DOS : ',1pe9.2,' per state.')
      do  iq = 1, nqp
        wt = abs(wgts(iq)) / nsp
        do  iband = 1, nband
          do  isp = 1, nsp
          e = evl(iband,isp,iq)
          meshpt = (e - emin) / step
          mesh1 = meshpt - mrange
          mesh2 = meshpt + mrange
          if (mesh2 >= ndos) mesh2 = ndos-1
          if (mesh1 < 0) mesh1 = 0
            do  meshpt = mesh1, mesh2
            emesh = emin + meshpt * step
            x = (emesh - e) / w
            call delstp(n,x,d,s,xx)
            dos(meshpt,isp) = dos(meshpt,isp) + wt * d / w
            enddo
          enddo
        enddo
      enddo
      end
