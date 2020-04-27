      subroutine hs1cx(nds,alat,plat,pos,
     .  exi,nbas,nl,ips,rsm,rsma,kmx,hcr,
     .  iax,ntab,s,cy,cg,indxcg,jcg,hsvx,hslx,hvx,hlx)
C- Value and Laplacian of screened solid smoothed Hankels at
C  the surface of augmentation spheres using one-center expansion,
C  all clusters
C ----------------------------------------------------------------
Ci Inputs
Ci   nds   :leading dimensions of s
Ci   rtab  :site positions corresponding to entries in iax table (ppair9)
Ci   exi   :energy of Hankel functions
Ci   nbas  :size of basis
Ci   nl    :max angular momentum + 1 (for dimensioning)
Ci   rsm   :Hankel smoothing radii
Ci   rsma  :smoothing radii for polynomial expansion
Ci   hcr   :augmentation sphere radii
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci         :iax(9,i) = (lmx+1)**2 for basis atom iax(2,i) (mkiaxd)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   s     :real-space structure constant matrix:
Ci         :H^a_RL = H_R'L' S_R'L',RL
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordan coefficients
Co Outputs
Co   hsvx,hslx :L-resolved value/Laplacian of screened smoothed Hankels
Co             :at the surface of augmentation spheres.
Co             :first index corresponds to the L-channel of the augmentation sphere,
Co             :second index enumerates screened functions for given cluster
Co             :third index runs over the atoms of the cluster as in iax table
Co   hvx,hlx   :same as above but for unsmoothed Hankels (for debugging)
Cl Local variables
Ci   iat   :Site R for which screened Hankel functions are calculated
Cl   ixi   :Site R' around which the value and Laplacian are obtained
Cr Remarks
Cr   Evaluates screened envelope function K_mRL at a point x-R
Cr   The function is defined as
Cr     Hs^a_mRL(x) = sum_m'R'L' Hs_m'R'L' * S_m'R'L',mRL      (*)
Cr   where Hs_R'(x) is the smoothed solid Hankel function Hs(x-R')
Cr
Cr   So far only the "single kappa" case is implemented, ie (*) assumes m = 1
Cr
Cr   Note: 2nd generation strux (called "Salpha") are -B^alpha here.
Cr
Cr   The spherical harmonics for solid Hankels are defined wrt the
Cr   head of the respective cluster as this is the convention used
Cr   in constructing strux matrix S.
Cr
Cu Updates
Cu   18 Dec 08 (S Lozovoi) Adapted from strg1c.f
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nbas,nds,nl
      integer, parameter :: niax=10, nkap=1
      integer iax(niax,*),ntab(nbas),ips(nbas),kmx(*)
      double precision alat,plat(3,3),pos(3,*),exi,
     .  rsm(0:nl-1,*),hcr(0:nl-1,*),rsma(*),
     .  s(nds,nds,nkap,nkap,*),
     .  cy(*)
      double precision hsvx(nds,nds,*),hslx(nds,nds,*)
C Clebsch Gordan arrays
      integer indxcg(*),jcg(*)
      double precision cg(*)
C Local variables
      integer nclus,iat,offR,offRp1,offRj,iatRp,jat,ll,il,lj,
     .  ixi,ixj,isg,isj
      integer nlmh,nlmg,nlmj,kmaxg,ndim
      integer, parameter :: n0 = 10, nclusm = 200
      double precision rsmj(0:n0-1),hcrg(0:n0-1),ej(0:n0-1)
      double precision rsmag
      double precision pclus(3,nclusm)
c     double precision sx
C Automatic arrays
      double precision hvl(nds,nds,0:1),gvl(nds,nds,0:1)
c...deb
c... for unsmoothed Hankels
      double precision rsmj0(0:n0-1)
      double precision hvx(nds,nds,*),hlx(nds,nds,*)
      integer ideb
c     integer ih,ig,ij,ii,ixjj
c...deb

      call tcn('hs1cx')
C --- Checks ---
      if (nds > n0**2)
     .  call rxi('hs1cx: nds is bigger than n0**2. nds = ',nds)
      if (ntab(1) /= 0)
     .  call rxi('hs1cx: ntab(1) should be zero, not ',ntab(1))

      ndim = nds*nds*ntab(nbas+1)
      call dpzero(hsvx,ndim)
      call dpzero(hslx,ndim)
c...deb
c... for unsmoothed Hankels
      rsmj0(:) = 1d-3
c     rsmj0(:) = 1d-2
      call dpzero(hvx,ndim)
      call dpzero(hlx,ndim)
c...deb

C --- begin big loop over clusters
      do  iat = 1, nbas

C ... Offset to iax table for cluster connected to R
        offR  = ntab(iat)
        nclus = ntab(iat+1) - offR
        nlmh = iax(9,offR+1)

C --- More checks ---
        if (nclus <= 0)
     .    call rxi('hs1cx: empty cluster encountered for iat = ',iat)
        if (nclus > nclusm)
     .    call rxi('hs1cx: size of the cluster exceeds the maximum.'//
     .    ' nclus = ',nclus)

c ... find coordinates of all atoms in the cluster
        do  ixi = 1, nclus
          call acoord(iat,ixi,alat,plat,pos,iax,ntab,pclus(1,ixi))
        enddo

C --- begin loop over expansion sites R' = ixi
        do  ixi = 1, nclus
C         offRp1 is (offset to iax table for site R') + 1
          offRp1 = offR + ixi
C         Basis index, Lmax, kmax, and aug. radii for site R'
          iatRp = iax(2,offRp1)
          nlmg  = iax(9,offRp1)
          isg = ips(iatRp)
          kmaxg = kmx(isg)
          rsmag = rsma(isg)
          do  il = 0, ll(nlmg)
            hcrg(il) = hcr(il,isg)
          enddo

C --- H^a(R') = sum_R'' H(R' - R'') * S_R''L'', R'L' ---
          do  ixj = 1, nclus
            offRj = offR+ixj
            jat = iax(2,offRj)
            isj = ips(jat)
            nlmj  = iax(9,offRj)
            lj = ll(nlmj)
            do  il = 0, lj
              rsmj(il) = rsm(il,isj)
            enddo
c So far Hankel energies do not depend on either l or species
            call dvset(ej(0),1,lj+1,exi)

cC       Use x - R'' =  (x-R') + (R'-R) - (R''-R)
c          xmrpp = xmrp - rtab(:,offR+ixj) + rtab(:,offRp1)

c ... 1c-decomposition of smoothed Hankels around ixi
c         call hs1c(pclus(1,ixj),rsmj,nlmj,exi,
c    .      pclus(1,ixi),rsmag,nlmg,kmaxg,hcrg,
c    .      cy,cg,indxcg,jcg,nds,hsv,hsl)
c in fact we don't need Gaussians, but leave as it is for now
            call rx('rework call to gh1c')
C            call gh1c(11,pclus(1,ixj),rsmj,nlmj,ej,xx,
C     .        pclus(1,ixi),rsmag,nlmg,kmaxg,hcrg,
C     .        cy,cg,indxcg,jcg,nds,gvl,hvl)

c           call tcn('direct summation')
c            do  ij = 1, nlmj
c              do  ih = 1, nlmh
cc only single kappa case
c                sx = s(ij,ih,1,1,offR+ixj)
cc...deb
cC       ixjj points to (R'',R) pair in strux
cc               ixjj = iax(8,offR+ixj)
cc               sx = s(ih,ij,1,1,ixjj)
cc...deb
c                do  ig = 1, nlmg
c                  hsvx(ig,ih,offRp1) = hsvx(ig,ih,offRp1) +
c     .              hvl(ij,ig,0)*sx
c                  hslx(ig,ih,offRp1) = hslx(ig,ih,offRp1) +
c     .              hvl(ij,ig,1)*sx
c                enddo
c              enddo
c            enddo
c           call tcx('direct summation')
c Same as above but using dgemm
            call tcn('using dgemm for hs')
            call dgemm('N','N',nlmg,nlmh,nlmj,1d0,hvl(1,1,0),nds,
     .        s(1,1,1,1,offR+ixj),nds,1d0,hsvx(1,1,offRp1),nds)
            call dgemm('N','N',nlmg,nlmh,nlmj,1d0,hvl(1,1,1),nds,
     .        s(1,1,1,1,offR+ixj),nds,1d0,hslx(1,1,offRp1),nds)
            call tcx('using dgemm for hs')
c...deb
c Also make unsmoothed Hankels (= Hsm with small rsm)
            ideb = 1
            if (ideb == 0) then
c             call gh1c(1,pclus(1,ixj),rsmj0,nlmj,ej,
            call rx('rework call to gh1c')
C              call gh1c(11,pclus(1,ixj),rsmj0,nlmj,ej,xx,
C     .          pclus(1,ixi),rsmag,nlmg,kmaxg,hcrg,
C     .          cy,cg,indxcg,jcg,nds,gvl,hvl)
            else
              call h1c(pclus(1,ixj),nlmj,ej,pclus(1,ixi),nlmg,hcrg,
     .          cy,cg,indxcg,jcg,nds,hvl)
c                 print '(a,3f12.4,i5)',' hs1cx: return from h1c'
            endif

            call tcn('using dgemm for h')
            call dgemm('N','N',nlmg,nlmh,nlmj,1d0,hvl(1,1,0),nds,
     .        s(1,1,1,1,offR+ixj),nds,1d0,hvx(1,1,offRp1),nds)
            call dgemm('N','N',nlmg,nlmh,nlmj,1d0,hvl(1,1,1),nds,
     .        s(1,1,1,1,offR+ixj),nds,1d0,hlx(1,1,offRp1),nds)
            call tcx('using dgemm for h')
          enddo
c...deb
c ...   end loop over expansion sites
        enddo
c ... end loop over clusters
      enddo

      call tcx('hs1cx')
      end
