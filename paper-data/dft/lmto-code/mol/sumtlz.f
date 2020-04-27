      subroutine sumtlz(rhoi,poti,pot0,nxi,lxi,exi,n0,rint,
     .  nbas,nbasp,lmaxp,alat,pos,ips,ioff,ioffv0,cg,jcg,indxcg,cy,
     .  ib,r1,nr,nr1,rofi,nlml,nsp,ixp,dxp,vval,rl,vl,xi,sxi0)
C- Sum tails of interstitial estatic pot and density at site ib
C ----------------------------------------------------------------------
Ci Inputs
Ci   rhoi  :coefficients to density in atom-centered density basis
Ci   poti  :coefficients to estat potential in atom-centered density basis
Ci   pot0  :coefficients to estat potential, homogeneous contribution
Ci   nxi   :number of functions per l for c.d. basis, by species
Ci   lxi   :l cutoffs for each function in c.d. basis, by species
Ci   exi   :Hankel energies for each function in c.d. basis, by species
Ci   n0    :dimensioning parameter
Ci   rint  :range of interstial charge density
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis.
Ci         :Sites nbas+1..nbasp are treat as classical point multipoles.
Ci         :sumtlz adds to the electrostatic potential the contributions
Ci         :from these multipoles, and makes vval for them
Ci   lmaxp :(nbasp>nbas) l-cutoff for point multipoles
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   pos   :basis vectors
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   ioff  :table of offsets to rhoi and poti for each site
Ci   ioffv0:table of offsets to vval and pot0 for each site
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   cy    :Normalization constants for spherical harmonics
Ci   r1    :(augmentation radius for site ib, ib<=nbas
Ci   nr    :number of radial mesh points at site ib, ib<=nbas
Ci   nr1   :leading dimension of vl(ib)
Ci   rofi  :radial mesh points at ib
Ci   nlml  :L cutoff for vl
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ixp   :0 for molecules branch
Ci         :Otherwise, parameters related to Ewald summations
Ci   dxp   :alat for molecules branch
Ci         :Otherwise, parameters related to Ewald summations
Ci   xi    :work array of dimension nr*(lmaxl+1)
Co Outputs
Co   rl    :envelope density on radial mesh
Co   vl    :envelope electrostatic potential on radial mesh
Co   vval  :(ib<=nbas) estat potential at r1
Co         :(ib>nbas) estat potential at site ib / r**l for r->0
Co   sxi0  :integral of fit function in the interstitial
Cu Updates
Cu   27 Jun 01
Cu   30 Dec 94 Made spin pol (MvS).
Cu   22 Jun 01 Adds zero energy bessel tails from point multipoles (nbasp)
C ----------------------------------------------------------------------
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(n0,1),nxi(1),ips(1),lxi(n0,1),ioff(1),ioffv0(1),
     .  poti(1),pos(3,1),rint(1),cg(1),jcg(1),indxcg(1),cy(1),dr(3),
     .  pot0(1),vl(nr1,1),rl(nr,1),rofi(nr),rhoi(1),vval(nlml),
     .  sxi0(1),xi(nr,0:1),ojkl(0:20),ixp(1),dxp(1)
      integer w(1)
      common /w/ w
      call getpr(ipr)
      lmaxl = ll(nlml)
      nri = ioff(nbas+1)
      if (ib <= nbas) then
        call dpzero(vl,    nlml*nr1)
        call dpzero(rl,    nlml*nr*nsp)
        call defrr(oh,     nr)
        call defrr(oy,     nr)
      else
        call defrr(oh,     1)
      endif

C ... Setup for strux
      call strxsu(nlml,nxi,lxi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,          nlml*nlmbx)

C --- Start loop over other atoms ---
      do  20  jb = 1, nbasp
        call dpdist(pos(1,ib),pos(1,jb),3,dd)
        if(alat*dd < 1d-4 .and. ixp(1) == 0) goto 20
        dr(1) = alat*(pos(1,jb)-pos(1,ib))
        dr(2) = alat*(pos(2,jb)-pos(2,ib))
        dr(3) = alat*(pos(3,jb)-pos(3,ib))
        lmax  =  lmaxp
        if (jb <= nbas) then
          lmax = -1
          js = ips(jb)
          do 2 je = 1,nxi(js)
            lmax = max0(lmax,lxi(je,js))
    2     continue
C     ... Loop over energies, add bessel tails
          call defrr(ohl,nlmp*nxi(js))
          if(alat*dd < rint(js)+r1) then
            j1 = ioff(jb)+1
            call rstr0(nxi(js),lxi(1,js),exi(1,js),nlmp,1,
     .        dr(1),dr(2),dr(3),lmaxl,1,w(ohl),w(ohl))
            do 11 je = 1,nxi(js)
              e = exi(je,js)
              lb = lxi(je,js)
              nlmb = (lb+1)**2
              if(nlmb > nlmbx) call rx('sumtlz: increase nlmbx')
C|    call mstrux(e,dr,w(os),nlml,nlmb,nlml,cg,indxcg,jcg,cy)
              ojj = ojkl(lb)
C     call nstrux(e,dr,nlml,nlmb,nlmp,npow,w(oikl),w(ojj),w(oip),
C    .   w(ocf),cy,w(os))
C              ih = nlmp*(je-1)
C              call hstrux(e,nlml,nlmb,nlmp,npow,ih,w(oikl),w(ojj),
C     .                    w(oip),w(ocf),w(ohl),w(os))
              call hstrux(e,nlml,nlmb,nlmp,npow,je,je,w(oikl),w(ojj),
     .                     w(oip),w(ocf),w(ohl),w(os))
              if (ib <= nbas) then
              call ropbes(rofi,e,lmaxl,w(oy),w(oh),xi,nr,1)
              call sumtl1(nlml,nlmb,w(os),nr,nr1,poti(j1),xi,w(oh),rofi,
     .                    vl)
              call sumtl1(nlml,nlmb,w(os),nr,nr,rhoi(j1),xi,w(oh),rofi,
     .                    rl)
              if (nsp == 2) call sumtl1(nlml,nlmb,w(os),nr,nr,
     .                         rhoi(j1+nri),xi,w(oh),rofi,rl(1,1+nlml))
              call sumtl2(nlml,nlmb,w(os),e,r1,sxi0(j1))
              else
                call sumtl3(nlml,nlmb,w(os),poti(j1),vval)
              endif
              j1 = j1+nlmb
   11       continue
          endif
          call rlse(ohl)
        endif
C   --- Add bessel tails with energy zero ---
        j0 = ioffv0(jb)+1
        nlmb = (lmax+1)**2
C|    call mstrux(0d0,dr,w(os),nlml,nlmb,nlml,cg,indxcg,jcg,cy)
        ojj = ojkl(lmax)
C      call nstrux(0d0,dr,nlml,nlmb,nlmp,0   ,w(oikl),w(ojj),w(oip),
C     .   w(ocf),cy,w(os))
        call defrr(ohl,nlmp)
        call rstr0(1,lmax,0d0,nlmp,1,dr(1),dr(2),dr(3),
     .    lmaxl,1,w(ohl),w(ohl))
C        call hstrux(0d0,nlml,nlmb,nlmp,0,0,w(oikl),w(ojj),w(oip),
C     .    w(ocf),w(ohl),w(os))
        call hstrux(0d0,nlml,nlmb,nlmp,0,1,1,w(oikl),w(ojj),w(oip),
     .    w(ocf),w(ohl),w(os))
        call rlse(ohl)
        if (ib <= nbas) then
          call ropbes(rofi,0d0,lmaxl,w(oy),w(oh),xi,nr,1)
          call sumtl1(nlml,nlmb,w(os),nr,nr1,pot0(j0),xi,w(oh),rofi,vl)
        else
          call sumtl3(nlml,nlmb,w(os),pot0(j0),vval)
        endif
   20 continue

      if (ib <= nbas) then
        do  ilm = 1, nlml
          vval(ilm) = vval(ilm)+vl(nr1,ilm)
        enddo
      endif
      call rlse(oh)
      end
      subroutine sumtl1(nlml,nlmb,s,nr,nr1,poti,xi,h,rofi,vl)
C- One-center expansion at ib of interstitial pot from Hankels at jb
      implicit real*8 (a-h,p-z), integer (o)
      dimension vl(nr1,1),rofi(1),xi(nr,0:1),s(nlml,nlmb),
     .   poti(1),h(1)
      lmaxl = ll(nlml)
      do  1  i = 1, nr1
    1 h(i) = 1d0
      ilml = 0
      do  16  l = 0, lmaxl
        do  14  m = -l, l
          ilml = ilml+1
          sum = 0d0
          do  15  ilmb = 1, nlmb
   15     sum = sum + s(ilml,ilmb)*poti(ilmb)
          do  3  i = 1, nr1
    3     vl(i,ilml) = vl(i,ilml) + sum*h(i)*xi(i,l)
   14   continue
        do  2  i = 1, nr1
    2   h(i) = h(i)*rofi(i)
   16 continue
      end

      subroutine sumtl2(nlml,nlmb,s,e,r,sxi0)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s(nlml,nlmb),sxi0(1),psi(0:3),phi(0:3)
      srfpi = dsqrt(16d0*datan(1d0))
      call bessl(e*r*r,1,phi,psi)
      sumj = srfpi*phi(1)*r**3
      do 1 ilm = 1,nlmb
  1   sxi0(ilm) = sxi0(ilm) - s(1,ilm)*sumj
      end
      subroutine sumtl3(nlml,nlmb,s,poti,vval)
C- One-center expansion at ib, r->0, of istl pot from Hankels at jb
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlml  :L cutoff for vval
Ci   nlmb  :L cutoff for poti
Ci   s     :real-space structure constant matrix
Ci   poti  :coefficients to estat potential in atom-centered density basis
Co Outputs
Co   vval  :L component of potential / r^l for r->0 added to site ib
Cr Remarks
Cr   sumtl3 is the analog of sumtl1 for spheres of zero radius, but with
Cr   nonzero multipole moments.
Cr   It returns true vval/r**l, which is finite as r->0
Cu Updates
Cu   27 Jun 01
C ----------------------------------------------------------------------
      implicit none
      integer nlml,nlmb
      double precision s(nlml,nlmb),poti(nlmb),vval(nlml)
      integer lmaxl,l,m,ilml,ilmb,ll
      double precision sum,df
      lmaxl = ll(nlml)
      df = 1d0
      ilml = 0
      do  16  l = 0, lmaxl
        df = df*(2*l+1)
        do  14  m = -l, l
          ilml = ilml+1
          sum = 0d0
          do  15  ilmb = 1, nlmb
   15     sum = sum + s(ilml,ilmb)*poti(ilmb)/df
          vval(ilml) = vval(ilml) + sum
   14   continue
   16 continue
      end
