      subroutine dostet(nbmx,nsp,nspx,nevmx,nchan,n1,n2,n3,ntet,
     .  idtet,eband,doswt,npts,emin,emax,lidos,wk,zos)
C- Density of states to third order in k spacing by tetrahedron method
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbmx  :leading dimension of eband
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspx  :number of independent spin channels (1 unless nsp=2, independent spins)
Ci   nevmx :number of bands
Ci   nchan :number of DOS channels
Ci   n1..n3:number of k divisions in each direction: only n1*n2*n3 is used
Ci   ntet  :number of tetrahedra and (tetirr.f)
Ci   idtet :idtet(1..4,i) the irreducible k-points defining the ith tetrahedron (tetirr.f)
Ci   eband :energy bands
Ci   doswt :DOS weights
Ci   npts  :number of DOS tabulation points in energy window [emin, emax];
Ci   emin  :energy window [emin, emax]
Ci   emax  :energy window [emin, emax]
Ci   lidos :F zos = dos
Ci         :T zos = energy integral of dos
Ci   wk    :work array
Co Outputs:
Co   zos : DOS (or integrated DOS for lidos=T) for each spin and nchan
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nchan,nsp,nspx,nbmx,npts,ntet,idtet(0:4,*),n1,n2,n3,nevmx
      double precision eband(nbmx,nspx,1),emin,emax,wk(npts),
     .  zos(npts,nsp,nchan),doswt(nchan,nbmx,nsp,*)
      logical lidos
C ... Local parameters
      integer isp,ib,i,itet,ichan,iq1,iq2,iq3,iq4,nspc,jsp,ksp
      double precision bin,eigen(4),v,wt,ebot,dmin1,eg0(4)
C ... External calls
      external dcopy,dpzero,slinz

      if (npts <= 1 .or. npts <= 2 .and. .not. lidos) call rx1(
     .  'dostet: npts(=%i) too small for DOS : require npts>2',npts)
      nspc = nsp / nspx
      call dpzero(zos,npts*nsp*nchan)
      bin = npts - 1
      bin = (emax - emin) / bin
      v = ( 3d0  -  nsp ) / ( n1 * n2 * n3 * 6d0 )

C --- Loop over tetrahedra ---
      do  itet = 1, ntet
        iq1 = idtet(1,itet)
        iq2 = idtet(2,itet)
        iq3 = idtet(3,itet)
        iq4 = idtet(4,itet)

C --- Loop over spins and sum over bands ---
        do  isp = 1, nspx
          do  ib = 1, nevmx
            eigen(1) = eband(ib,isp,iq1)
            eigen(2) = eband(ib,isp,iq2)
            eigen(3) = eband(ib,isp,iq3)
            eigen(4) = eband(ib,isp,iq4)
            call dcopy(4,eigen,1,eg0,1)
            ebot = dmin1(eigen(1),eigen(2),eigen(3),eigen(4))

            if (ebot <= emax) then
            do  jsp = 1, nspc
C       ... ksp is isp for uncoupled spins, and jsp for coupled spins
            ksp = max(jsp,isp)

C       ... Accumulate no. states assuming constant wt from this tet
            do  ichan = 1, nchan
              wt = doswt(ichan,ib,ksp,iq1)
     .           + doswt(ichan,ib,ksp,iq2)
     .           + doswt(ichan,ib,ksp,iq3)
     .           + doswt(ichan,ib,ksp,iq4)
              wt = wt * idtet(0,itet) * v / 4d0
              call slinz(wt,eigen,emin,emax,zos(1,ksp,ichan),npts)
            enddo
            enddo
            endif
          enddo
        enddo
      enddo

      if (lidos) return

C --- DOS from finite difference of NOS ---
      bin = 2d0 * bin
      do  isp = 1, nsp
        do  ichan = 1, nchan
          do  i = 2, npts-1
            wk(i) = (zos(i+1,isp,ichan)-zos(i-1,isp,ichan))/bin
          enddo
        wk(1)    = wk(2)
        wk(npts) = wk(npts-1)
        call dcopy(npts,wk,1,zos(1,isp,ichan),1)
        enddo
      enddo
      end
