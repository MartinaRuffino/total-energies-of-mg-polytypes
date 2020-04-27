      subroutine addbkgsm(smrho,k1,k2,k3,nsp,qbg,vol,fac)
C- add uniform background to smooth charge density
Ci Inputs: smrho: smooth density
Ci         qbg: background charge
Ci         k1,k2,k3, dimensions of smrho mesh
Ci         nsp: number of spins
Ci         vol: volume per cell
Ci         fac: fac*qbg/vol//nsp is added to
C-------------------------------------------
      implicit none
      integer j1,j2,j3,k1,k2,k3,nsp,is
      double precision qbg,rhobg,vol,fac
      double complex smrho(k1,k2,k3,2)
      rhobg=qbg/vol
      do is=1,nsp
       do j1=1,k1
        do j2=1,k2
         do j3=1,k3
            smrho(j1,j2,j3,is)=smrho(j1,j2,j3,is)+rhobg*fac/nsp
         enddo
        enddo
       enddo
      enddo
      end
      subroutine adbkql(nbas,nsp,qbg,vol,fac,s_spec,s_site)
C- Add uniform bkg charge density to local smooth rho
C----------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr rmt lmxl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2
Cio    Passed to:  *
Ci nbas: number of atoms in basis
Ci qbg: background charge
Ci nsp: spins
Ci vol: vol of cell
Ci fac: fac * backg density is added
Cu Updates
Cu   11 Dec 12 Migrated to f90 structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 Zero-radius sites skipped over
C----------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nrmx,nlmx,nlml,lmxl,nbas,nsp
      parameter (nrmx=5001,nlmx=64)
      double precision qbg,fac
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,nr,is
      double precision rhobkg,vol,a,rmt,rofi(nrmx)

      rhobkg = fac*qbg/vol
      do  ib = 1, nbas
       is = s_site(ib)%spec
       a = s_spec(is)%a
       nr = s_spec(is)%nr
       rmt = s_spec(is)%rmt
       lmxl = s_spec(is)%lmxl
       if (lmxl == -1) goto 10
       nlml=(lmxl+1)**2
C      nrml=nr*nlml
       call rxx(nr > nrmx,  'addbkgloc: increase nrmx')
       call rxx(nlml > nlmx,'addbkgloc: increase nlmx')
       call radmsh(rmt,a,nr,rofi)
       call addbkgl(s_site(ib)%rho1,s_site(ib)%rho2,rhobkg,
     .   nr,nsp,rofi,nlml)
   10  continue
      enddo
      end
      subroutine addbkgl(rho1,rho2,rhobkg,nr,nsp,rofi,nlml)
C adds uniform background to local smooth density at this site for
C l=0 component (ilm=1)
C for each spin
      implicit none
      integer nsp,is,nr,nlml,i
      double precision rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rofi(nr)
      double precision rhobkg,pi,srfpi
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
C     y0 = 1d0/srfpi
      do is = 1, nsp
        do i = 1, nr
          rho1(i,1,is) = rho1(i,1,is)+srfpi*rofi(i)**2*rhobkg/nsp
          rho2(i,1,is) = rho2(i,1,is)+srfpi*rofi(i)**2*rhobkg/nsp
        enddo
      enddo
      end


