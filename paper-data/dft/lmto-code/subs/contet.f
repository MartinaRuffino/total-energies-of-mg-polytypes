      subroutine contet(mode,nbmx,nsp,nspx,nevmx,nchan,n1,n2,n3,ntet,idtet,
     .  alat,plat,qp,ipq0,igstar,iblst,eband,cvec,doswt,npts,emin,emax,lidos,ief,zos)
C- Conductivity or DOS-related quantity by tetrahedron integration
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :selects what function is to be integrated (see mkcond.f)
Ci         :1s digit:
Ci         :0 => dos
Ci         :1 => 1/2 | grad_k E(k) . cvec | (ballistic conductivity)
Ci         :2 => grad_1 E(k) . grad_2 E(k) where _1 and _2 are
Ci         :     direction vectors specified by cvec
Ci         :10s digit:
Ci         :0 normal mode
Ci         :1 Write k-resolved data for tetrahedra which enclose ef
Ci   nbmx  :dimensions eband,doswt
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspx  :number of independent spin channels
Ci          (1 unless nsp=2, and independent spins)
Ci   nevmx :number of bands to sum over
Ci   nchan :number of DOS channels
Ci   n1..n3:no. divisions for the 3 recip. lattice. vectors
Ci   ntet  :no. of different tetrahedra
Ci   idtet :idtet(0,i) =no of tetrahedra of the i'th kind
Ci         :idtet(1-4,i) marks the 4 irreducible qp of i'th tetrahedron
Ci   qp    :list of irreducible k-points: qp(1..3,i) is i'th qp
Ci   ipq0  :ipq as generated by bzmesh for qp of the irr. BZ.
Ci   igstar:table of inverse mapping of ipq for the BZ of reduced
Ci         :symmetry; generated by call to bzmesh with igstar(0)=2
Ci         :and with no symmetry operations.
Ci   iblst :(iblst(1)>0)  => iblst = a prescribed list of energy bands
Ci         :In this mode sum is bands iblst(1)..iblst(nevmx).
Ci         :(iblst(1)<=0) => sum over bands 1..nevmx
Ci   eband :energy bands
Ci   cvec  :direction in which to calc. conductivity
Ci   doswt :number of states, for each of nchan channels, at each energy
Ci         :and irr. qp of whatever to be integrated by the tetrahedron method.
Ci   npts  :number of tabulation points in energy range (emin,emax)
Ci   emin  :lower bound to energy window
Ci   emax  :upper bound to energy window
Ci   ef    :(used only if 10s digit mode = 1)
Ci         :keep k-resolved data for for tetrahedra that straddle ef
Ci   lidos :F zos = conductivity or dos-related quantity
Ci         :T zos = energy integral of this quantity
Ci   wk    :work array of size npts
Co Outputs:
Co   zos   :conductivity or other dos-related quantity on uniform
Co         :mesh in energy window (emin,emax); or energy integral
Co         :of the same (lidos=T); see Remarks
Cr Remarks
Cr   This routine uses the tetrahedron method to integrate quantities
Cr       zos_n(E) = int d^3k delta(E-E(k)) doswt_n (k) f(E,k)
Cr   foreach n = 1..nchan.  doswt_n is the number of states in channel
Cr   n; thus if f(k)=1 contet returns the density of states for each
Cr   channel.  The caller chooses what function f is to be integrated
Cr   by input mode (see mkcond.f)
Cu Updates
Cu   24 Feb 15 Bug fix for noncollinear case
Cu   03 Feb 01 several revisions: different modes for mkcond
Cu             integrate over subset of bands
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nchan,nsp,nspx,nbmx,npts,ntet,idtet(0:4,ntet),n1,n2,
     .  n3,nevmx,igstar(0:*),ipq0(n1,n2,n3),iblst(nevmx),ief
      real(8) eband(nbmx,nspx,*),emin,emax,plat(3,3),alat,
     .  zos(npts,nsp,nchan),doswt(nchan,nbmx,nsp,*),qp(3,*),cvec(3,2)
      logical lidos
C ... Dynamically allocated local arrays
      integer,allocatable :: ncrossi(:,:),ncrosst(:,:),iblstl(:)
      real(8),allocatable :: wk(:),dosk(:,:),dosib(:),qpl(:,:)
C ... Local parameters
      logical lkres
      integer isp,ib,i,itet,ichan,iq(4),jq(4),nspc,jsp,ksp,i1,i2,i3,ibp,
     .  ipass,ncross,ifi,iqmax,mode0
      real(8) bin,eigen(4),v,wt,ebot,etop,dmin1,wc,ef,qpi(3),qlat(3,3),automks,aby2pi
      real(8), parameter :: tol=1d-6
      character*7 strmod(0:3)
      data strmod /'dos ','<v/2.c>','<v^2>','<|v|>'/

      procedure(integer) fopna
C ... External calls
      external dcopy,dpzero,mkcond,slinz

      if (npts <= 1 .or. npts <= 2 .and. .not. lidos) call rx1(
     .  'contet: npts(=%i) too small for DOS : require npts>2',npts)
      nspc = nsp / nspx
      call dpzero(zos,npts*nsp*nchan)
      bin = npts - 1
      bin = (emax - emin) / bin
      v = ( 3d0  -  nsp ) / ( n1 * n2 * n3 * 6d0 )
      lkres = mod(mode/10,10) > 0
      ef = ief*bin + emin
      if (lkres .and. nchan>1) call rx('contet" k-res can only be used with nchan=1')
      mode0 = mod(mode,10); if (mode0 > 3) mode0 = 0

C --- Make passes when 10's digit mode > 0 ---
      do  ipass = 2, 1, -1
      if (ipass == 2) then
        if (.not. lkres) cycle  ! Only one pass in normal mode
        allocate(ncrossi(nevmx,nsp))
        ncross = 0; call iinit(ncrossi,size(ncrossi))
      elseif (lkres) then
        call info5(30,1,0,' CONTET: %i tetra enclose ef=%;6d.  Resolved by band: %n:1i',
     .    ncross,ef,nevmx,ncrossi,4)
        allocate(dosk(nsp,ncross),ncrosst(ncross,2),dosib(nevmx))
        ncross = 0; call iinit(ncrossi,size(ncrossi))
        call dpzero(dosk,size(dosk))
      endif

C --- Loop over irreducible tetrahedra ---
      iqmax = 0
      do  itet = 1, ntet
        iq(1) = idtet(1,itet)
        iq(2) = idtet(2,itet)
        iq(3) = idtet(3,itet)
        iq(4) = idtet(4,itet)
        iqmax = max(iqmax,maxval(iq))
        
C   ... For each of the iq, find mapping to equivalent jq in eband
        do  i = 1, 4
          call cmpi123(0,i1,i2,i3,igstar(iq(i)))
          jq(i) = ipq0(i1,i2,i3)
        enddo

C --- Loop over spins and sum over bands ---
        do  isp = 1, nspx
          do  ib = 1, nevmx
            if (iblst(1) > 0) then
              ibp = iblst(ib)
            else
              ibp = ib
            endif
            eigen(1) = eband(ibp,isp,jq(1))
            eigen(2) = eband(ibp,isp,jq(2))
            eigen(3) = eband(ibp,isp,jq(3))
            eigen(4) = eband(ibp,isp,jq(4))

            ebot = dmin1(eigen(1),eigen(2),eigen(3),eigen(4))
            etop = dmax1(eigen(1),eigen(2),eigen(3),eigen(4))
C       ... Setup for k-resolved info for analysis near Ef.
            if (lkres) then
              if (ef<=etop .and. ef>=ebot) then
                ncross = ncross+1
                ncrossi(ib,isp) = ncrossi(ib,isp) + 1
              endif
            endif

            if (ipass == 1 .and. ebot <= emax) then

              do  jsp = 1, nspc
C       ...   ksp is isp for uncoupled spins, and jsp for coupled spins
              ksp = max(jsp,isp)

C             mkcond reorders eigen, so it should not be called again when jsp=2.
              if (jsp == 1) call mkcond(mode0,plat,eigen,qp,iq(1),iq(2),iq(3),iq(4),cvec,wc)

C         ... Accumulate k-resolved info for analysis near Ef.
              if (lkres .and. ef<=etop .and. ef>=ebot) then
                dosk(ksp,ncross) = wc
                ncrosst(ncross,1) = itet
                ncrosst(ncross,2) = ibp
              endif

C         ... Accumulate no. states assuming constant wt from this tet
              do  ichan = 1, nchan
C           ... This is weight for no. states
                wt = doswt(ichan,ibp,ksp,jq(1))
     .             + doswt(ichan,ibp,ksp,jq(2))
     .             + doswt(ichan,ibp,ksp,jq(3))
     .             + doswt(ichan,ibp,ksp,jq(4))
                wt = wt * idtet(0,itet) * v / 4d0
C         ..  . Add weights from mkcond
                wt = wt*wc
                call slinz(wt,eigen,emin,emax,zos(1,ksp,ichan),npts)
              enddo
            enddo
            endif

          enddo
        enddo
      enddo
      enddo ! Extra pass when 10s digit mode > 0

      if (lidos) return

C --- DOS from finite difference of NOS ---
      allocate(wk(npts))
      bin = 2d0 * bin

C     Shorten q to print midpoints in file dosq
      if (lkres) then
        allocate(qpl(3,iqmax))
        call mkqlat(plat,qlat,wc) ! ... qlat = (plat^-1)^T so that qlat^T . plat = 1
        do  i = 1, iqmax
          call shorbz(qp(1,i),qpl(1,i),qlat,plat)
        enddo
      endif

      do  isp = 1, nsp

        do  ichan = 1, nchan
          do  i = 2, npts-1
            wk(i) = (zos(i+1,isp,ichan)-zos(i-1,isp,ichan))/bin
          enddo
        wk(1)    = wk(2)
        wk(npts) = wk(npts-1)
        call dcopy(npts,wk,1,zos(1,isp,ichan),1)
        enddo

        if (lkres) then
          call dpzero(dosib,size(dosib))
          ifi = fopna('dosq',-1,0)
          call awrit5('# q+band -resolved %?;n==0;DOS;;%-1j%?;n==1;<v>;;%-1j%?;n==2;<v^2>;;'//
     .      '  %i bands  %i qp  ef=%,6;6d  microcell vol %,6;6g',' ',120,ifi,mode0,nevmx,iqmax,ef,v)
          write(ifi,"('#',12x,'-----  q  -----',14x,'tet  ib     DOS')")
          allocate(iblstl(nevmx))
          do  ib = 1, nevmx
            if (iblst(1) > 0) then
              ibp = iblst(ib)
            else
              ibp = ib
            endif
            iblstl(ib) = ibp
            do  i = 1, ncross
              if (ibp /= ncrosst(i,2)) cycle
              wt = dosk(isp,i)
              itet = ncrosst(i,1)
              iq(1) = idtet(1,itet); iq(2) = idtet(2,itet); iq(3) = idtet(3,itet); iq(4) = idtet(4,itet)
              qpi(:) = (qpl(:,iq(1)) + qpl(:,iq(2)) + qpl(:,iq(3)) + qpl(:,iq(4)))/4
              write(ifi,333) qpi, itet, ibp, wt
  333         format(3f12.6,i9,i4,f12.6)
              dosib(ib) = dosib(ib) + wt
            enddo
          enddo
          wt = 0
          do  ib = 1, nevmx
            wt = wt + dosib(ib)/ncross
            if (ncrossi(ib,isp) == 0) cycle
            dosib(ib) = dosib(ib)/ncrossi(ib,isp)
          enddo

          aby2pi = alat/(8*atan(1d0))
          automks = 13.605/6.5822e-16*0.5292d-10 * aby2pi / 1d6  ! convert to 10^6 m/sec

          call info5(10,1,0,' '//trim(strmod(mode0))//' (spin %i) at Ef = %;4d a.u.'//
     .      '%?#n==1|n==3# = %;4d x 10^6 m/s.#%j#'//
     .      '  Resolve by band:',isp,wt,mode0,wt*automks,5)

          call arrprt('  ib   '//strmod(mode0),'%,4i%,4;9D','id',nevmx,0,5,
     .      0,'  | ',iblstl,dosib,qpi,qpi,qpi,qpi,qpi,qpi)
          deallocate(iblstl)
          if (ief>0 .and. ief<npts) then
            call info2(2,0,0,' total DOS(Ef) = %,6;6d',zos(ief,isp,1),2)
          endif
        endif
      enddo
      deallocate(wk)
      if (allocated(ncrossi)) deallocate(ncrossi)

      end
