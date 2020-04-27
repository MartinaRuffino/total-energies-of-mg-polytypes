      subroutine v0intr(nl,nsp,lmax,z,rhozbk,rmax,a,nr,rofi,pnu,qnu,pz,
     .  bscal,v,rhoin,rho,rhoc,g,gp,nmix,niter,qc,lfrz,avw,ekap,mode,vintra)
C- Make derivatives of total energy or C w.r.t. q_j
C ----------------------------------------------------------------
Ci Inputs
Ci   Same as atomsc, with the following additions:
Ci   mode  1  Make dC_i/dq_j
Ci         2  Make d^2 E /dq_i dq_j
Co Outputs
Co   vintra   (mode 1) dC_i/dq_j
Co            (mode 2) d^2 E /dq_i dq_j
Cu Updates
Cu   03 Jun 14 Enable deep local and valence orbitals both to have charge
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nr,nsp,nl,nmix,niter,lmax,mode,lfrz
      double precision rofi(nr,2),v(nr,nsp),rho(nr,nsp),rhoc(nr,nsp),
     .  rhoin(nr,nsp),pnu(nl,nsp),qnu(3,nl,nsp),pz(nl),g(nr),gp(nr,4),
     .  vintra(nl,nl,nsp,nsp),rhozbk,avw,ekap,bscal
C ... Dynamically allocated local arrays
!      include 'structures.h'
      type(dp_wk_vec) vk(2)
      real(8), allocatable :: vl(:)
C ... Local parameters
      integer isp,jsp,il,jl,k,ki,kj,off,kl,ksp,iprint,jmin
      integer, parameter :: n0=10, ncmx=200, nvmx=20
      double precision sumec,sumtc,ekin,utot,vrmax(2),rhrmx,df,sumev,
     .  rhoeps,etot,amgm,qtot,ves,xx(2),thrpv,z,rmax,a,qc,exrmax(2),
     .  thrpvl(n0),etot0,pp(6,n0*2),cnl(2,n0),enl(2,2),fi,fj,qcut,
     .  pl(n0*2),ql(3,n0*2),ec(ncmx),ev(nvmx),t(n0,2),qnur(1),GFr(1)
C     double precision pz(n0,2),qz(3,n0,3)
      integer idmod(n0),idmoz(n0),idu(5)
      parameter (qcut=.05d0)
C ... Offset of index kl,ksp for array of dimension (nl,nsp)
      off(kl,ksp) = kl + nl*(ksp-1)

      call tcn('v0intr')
      if (nl > n0) call rx('v0intr: increase n0')
      call dpzero(vintra,nl*nl*nsp*nsp)
      call iinit(idu,5)
      ves = 0
      xx(1) = 0
      xx(2) = 0
      call dcopy(3*nl*nsp,qnu,1,ql,1)
      idmod(1:nl) = 0
      call iinit(idmoz,n0)
C     npan = 1
      df = .001d0
      allocate(vl(nr*nsp))
      allocate(vk(1)%p(nr*nsp))
      allocate(vk(2)%p(nr*nsp))

      call pshpr(iprint()-40)
      call getqvc(nsp,nl,lmax,z,pnu,qnu,pz,0,0,0,0,xx,qc,qtot,amgm,
     .  xx,xx)
      ec(1) = 0

C --- Derivatives of C_i wrt q_j ---
      if (mode == 1) then
        do  isp = 1, nsp
        do  jl = 1, nl
          do  k = 1, 2
            if (k == 1) fj= df
            if (k == 2) fj=-df
            call dcopy(nl*nsp,pnu,1,pl,1)
            call dcopy(nr*nsp,v,1,vl,1)
            ql(1,off(jl,isp)) = ql(1,off(jl,isp)) + fj

            call atomsc(1,nl,nsp,lmax,z,rhozbk,0,0,0d0,rmax,a,nr,
     .        rofi,ec,ev,pl,ql,qnur,pz,bscal,idmod,vl,rhoin,rho,
     .        rhoc,nmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,
     .        rhrmx,vrmax,qtot,exrmax,'pot',niter,lfrz,1d0)

            call potpar(nl,nsp,lmax,z,rmax,avw,ekap,.false.,.false.,
     .        .false.,a,nr,rofi,vl,pl,idmod,xx,ql,xx,idu,xx,xx,1d0,
     .        thrpv,thrpvl,g,gp,pp,xx,xx,xx,xx,t,GFr)

            ql(1,off(jl,isp)) = ql(1,off(jl,isp)) - fj
            do  il = 1, nl
              cnl(k,il) = pp(2,off(il,isp))
            enddo
          enddo
          do  il = 1, nl
            vintra(jl,il,isp,isp) = (cnl(1,il)-cnl(2,il))/2d0/df
          enddo
        enddo
        enddo

C --- d^2 E /dq_i dq_j ---
      elseif (mode == 2) then

C ... Make etot(dq=0)
      call dcopy(nl*nsp,pnu,1,pl,1)
      call dcopy(nr*nsp,v,1,vl,1)
      call atomsc(1,nl,nsp,lmax,z,rhozbk,0,0,0d0,rmax,a,nr,rofi,ec,
     .  ev,pl,ql,qnur,pz,bscal,idmod,vl,rhoin,rho,rhoc,nmix,qc,
     .  sumec,sumtc,sumev,ekin,utot,rhoeps,etot0,amgm,rhrmx,vrmax,qtot,
     .  exrmax,'pot',niter,lfrz,1d0)

      do  isp = 1, nsp
      do  jsp = isp, nsp
C ... Forget the off-diagonal for now
C     if (jsp /= isp) goto 20

      do  il = 1, nl
        jmin = il
        if (jsp /= isp) jmin = 1
        do  jl = jmin, nl

        vintra(jl,il,jsp,isp) = 0
        vintra(il,jl,isp,jsp) = 0
        if (ql(1,off(il,isp))+ql(1,off(jl,jsp)) < qcut/nsp) cycle

C   ... Accumulate E++, E+-, E-+, E--
        do  ki = 1, 2
          if (ki == 1) fi =  df
          if (ki == 2) fi = -df
          ql(1,off(il,isp)) = ql(1,off(il,isp)) + fi
          do  kj = 1, 2
            if (kj == 1) fj =  df
            if (kj == 2) fj = -df
            call dcopy(nl*nsp,pnu,1,pl,1)
            call dcopy(nr*nsp,v,1,vl,1)
C       ... Make a better guess for initial v on second pass
            if (ki == 2) then
              call dscal(nr*nsp,2d0,vl,1)
              call daxpy(nr*nsp,-1d0,vk(3-kj)%p,1,vl,1)
            endif
            ql(1,off(jl,jsp)) = ql(1,off(jl,jsp)) + fj
C       ... Self-consistent total energy for these moments
            call atomsc(1,nl,nsp,lmax,z,rhozbk,0,0,0d0,rmax,a,nr,
     .        rofi,ec,ev,pl,ql,qnur,pz,bscal,idmod,vl,rhoin,rho,
     .        rhoc,nmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,
     .        rhrmx,vrmax,qtot,exrmax,'pot',niter,1,1d0)
C            call awrit4('%x  ki,kj %i %i il jl %i %i',
C     .        outs,120,0,ki,kj,il,jl)
C            call awrit3('%a ql %,4;4d %,4;4d %,4;4d',
C     .        outs,120,0,ql(1,1),ql(1,2),ql(1,3))
C            if (nsp == 2)
C     .        call awrit3('%a %,4;4d %,4;4d %,4;4d',
C     .        outs,120,0,ql(1,4),ql(1,5),ql(1,6))
C            call awrit1('%a e %,10;10d',outs,120,-6,etot)
C       ... Setup for a better guess for initial v for second pass
            if (ki == 1) call dcopy(nr*nsp,vl,1,vk(kj)%p,1)
            ql(1,off(jl,jsp)) = ql(1,off(jl,jsp)) - fj
            enl(ki,kj) = etot
          enddo
          ql(1,off(il,isp)) = ql(1,off(il,isp)) - fi
        enddo

C   ... E'' by finite difference E++ - E+- - E-+ + E--
        vintra(jl,il,jsp,isp) = (enl(1,1) - enl(1,2) -
     .                           enl(2,1) + enl(2,2))/(2*df)**2
        vintra(il,jl,isp,jsp) = vintra(jl,il,jsp,isp)

C        call pshpr(10)
C        call info5(0,0,1,' il %i jl %i isp %i vintra %,6;6d',
C     .    il,jl,isp,vintra(jl,il,jsp,isp),0)
C        call poppr

        enddo
      enddo
      enddo
      enddo

      else
        call fexit(-1,111,' Exit -1 v0intr: bad mode',mode)
      endif

      deallocate(vl,vk(1)%p,vk(2)%p)
      call poppr

C      ltmp = aiovla(alabl,vintra,nl,nl-1,nsp,-6)
C      stop

      call tcx('v0intr')
      end
