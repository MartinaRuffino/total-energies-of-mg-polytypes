      subroutine subtls(nxi,lxi,exi,n0,rint,nbas,nbasp,lmaxp,alat,pos,
     .  ips,ioff,cg,jcg,indxcg,cy,r1,nr1,nr2,rofi,rwgt,nlml,nlmx,
     .  ixp,dxp,vl,vxc1,rl1,dvxc2,zeta,zetxc,ib,rhoi,rhois,poti,pot0,
     .  ioffv0,vval,qmom,upm,for1,for2,for4)
C- Subtracts integrals of xi-tails times pot from zeta,zetxc
C ----------------------------------------------------------------------
Ci Inputs
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
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   cy    :Normalization constants for spherical harmonics
Ci   r1    :augmentation radius at site ib
Ci   nr1   :number of radial mesh points at site ib
Ci   nr2   :number of radial mesh points in rim region at site ib
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nlml  :(local density l-cutoff + 1)**2
Ci   nlmx  :L cutoff istl density, site ib (dimensions qmom, pot0)
Ci   lmaxp  :dimensions local arrays.  Must be >= (max(lxi(ib)+1)**2
Ci   ixp   :0 for molecules branch
Ci         :Otherwise, parameters related to Ewald summations
Ci   dxp   :alat for molecules branch
Ci         :Otherwise, parameters related to Ewald summations
Ci   vl    :one-center expansion of smooth es potential at site ib
Ci   vxc1  :one-center expansion of smooth xc potential at site ib
Ci   rl1   :one-center expansion of smooth density at site ib
Ci   dvxc2 :difference in xc potential, true and smoothed istl density
Ci   ib    :site for which xi-tails are subtracted
Ci   rhoi  :coefficients to density in atom-centered density basis
Ci   rhois :rhoi, spin-up and spin down components combined
Ci   poti  :coefficients to estat potential in atom-centered density basis
Ci   pot0  :coefficients to estat potential, homogeneous contribution
Ci   ioffv0:table of offsets to vval and pot0 for each site
Ci   vval  :electrostatic potential at MT boundaries
Ci   qmom  :multipole moments of on-site densities
Co Outputs
Co   zeta  :integrals of density functions with electrostatic potential
Co         :Contribution to zeta from augmentation at ib is subtracted
Co   zetxc :integrals of density functions with xc potential.
Co         :Contribution to zetaxc from augmentation at ib is subtracted.
Co   upm   :interaction energy between total electrostatic potential
Co         :and point multipole at site ib(>nbas) added to upm.
Co   for1  :f1 + f7; see local variables for definitions
Co   for2  :(f2+f3+f4-f5)/2; see local variables for definitions
Co   for4  :f6 ; see local variables for definitions
Cl Local variables
Cl   Contributions to force:
Cl   f1    :d integral( rho * eps )
Cl   f2    :integral ( ves * drho)
Cl   f3    :integral( rho * dves )
Cl   f4    :force from smooth pot * change in qval inside sphere
Cl   f5    :force from qmom * change in es pot inside sphere
Cl   f6    :derivative of interstitial charge
Cl         :f6 gets scaled with vint later.
Cl   f7    :force from interstitial xc ring integrals
Cr Remarks
Cr
Cu Updates
Cu   4 Jan 95 spin polarized (MvS)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n0,nbas,nbasp,lmaxp,nlml,nlmx,nr1,nr2,ib,nxi(1),ips(1),
     .  lxi(n0,1),ioff(1),jcg(1),indxcg(1),ioffv0(1),ixp(1)
      double precision alat,r1,exi(n0,1),pos(3,1),rint(1),cg(1),cy(1),
     .  vl(nr1,nlml),rl1(nr1,nlml),vxc1(nr1,nlml,*),dvxc2(nr2,
     .  nlml,*),zeta(1),zetxc(1),rhoi(1),rhois(1),poti(1),pot0(1),
     .  vval(nlml),qmom(nlml),rwgt(1),rofi(1),for1(3,nbasp),
     .  for2(3,nbasp),for4(3,nbasp),dxp(1),upm
C ... Local parameters
      integer stdo,nlmxx,i0,i1,ih,ipr,is,j0,j1,jb,je,js,kr,lb,lgunit,ll,
     .  lmax,lmaxib,lmaxl,lmxh,lsp,m,nlmb,nlmbx,nlmh,nlmp,npow,nri,nrwk,
     .  nsp,ojkl(0:20)
      integer ocf,ogd,ogs,oh,ohd,ohl,oikl,oip,ojj,os,osd,oxi,oy
      parameter(nlmxx=49)
      double precision dd,e,fes,x0,sum(nlmxx),xval(nlmxx),xmom(nlmxx),
     .  dr(3),f1(3),f2(3),f3(3),f4(3),f5(3),f6(3),f7(3)
C ... Heap
      integer w(1)
      common /w/ w

      stdo = lgunit(1)
      call getpr(ipr)
      lmaxl = ll(nlml)
      lmaxib = ll(nlmx)
      nsp = lsp()+1
      nri = ioff(nbas+1)
      i0 = ioffv0(ib)+1
      if (ib <= nbas) then
        nrwk = max0(nr1,nr2)
        call defrr(oxi,nrwk*(lmaxl+1))
        call defrr(oh, nrwk)
        call defrr(oy, nrwk)
      else
        call defrr(oxi,1)
      endif

C ... Setup for strux
      call strxsu(nlmx,nxi,lxi,n0,0,nbas,ips,cg,jcg,indxcg,
     .  nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      if (nlmp < (lmaxp+1)**2) call rx('subtls: lmxp lt lmaxp')
      call defrr(os,  nlmx*nlmbx)
      call defrr(ogs, nlmx*nlmbx*3)
      call defrr(osd, nlmx*nlmbx)
      call defrr(ogd, nlmx*nlmbx*3)

C ... One-center electrostatic terms
      if (ib <= nbas) then
        i1 = ioff(ib)+1
        is = ips(ib)
        call stles1(r1,nxi(is),lxi(1,is),exi(1,is),poti(i1),pot0(i0),
     .    zeta(i1))
      endif

C --- Start loop over other spheres ---
      do  jb = 1, nbasp
        lmax = lmaxp
        call dpdist(pos(1,ib),pos(1,jb),3,dd)
        if (alat*dd < 1d-4 .and. ixp(1) == 0) goto 20
        dr(1) = alat*(pos(1,jb)-pos(1,ib))
        dr(2) = alat*(pos(2,jb)-pos(2,ib))
        dr(3) = alat*(pos(3,jb)-pos(3,ib))
        do  m = 1, 3
          f7(m) = 0d0
          f1(m) = 0d0
          f2(m) = 0d0
          f3(m) = 0d0
          f4(m) = 0d0
          f5(m) = 0d0
          f6(m) = 0d0
        enddo
        if (jb > nbas) then
          lmxh = lmaxib+lmax+1
          nlmh = (lmxh+1)**2
          goto 50
        endif

        js = ips(jb)
        lmax = -1
        do  je = 1, nxi(js)
          lmax = max0(lmax,lxi(je,js))
        enddo

C --- Loop over non-zero xi energies ---
        if (alat*dd <= rint(js)+r1) then
          j1 = ioff(jb)+1
          kr = nr1+1
          lmxh = lmaxib+lmax+1
          nlmh = (lmxh+1)**2
          call defrr(ohl,nlmh*nxi(js))
          call defrr(ohd,nlmh*nxi(js))
          call rstr0(nxi(js),lxi(1,js),exi(1,js),nlmh,1,
     .      dr(1),dr(2),dr(3),lmaxib+1,0,w(ohl),w(ohd))
C         lmax=-1
          do  je = 1, nxi(js)
            e = exi(je,js)
            lb = lxi(je,js)
C           lmax=max0(lmax,lb)
            nlmb = (lb+1)**2
            if (nlmb > nlmxx) call rx('subtlz: increase nlmbx')
            ojj = ojkl(lb)
C            call nstrpg(e,dr,nlmx,nlmb,nlmp,npow,w(oikl),w(ojj),w(oip),
C     .        w(ocf),cy,w(os),w(osd),w(ogs),w(ogd))
            ih = nlmh*(je-1)
            call hstrpg(e,nlmx,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
C       ... 2-center electrostatic terms
            if (ib <= nbas) then
              call stles2(r1,nxi(is),lxi(1,is),exi(1,is),nlmx,nlmb,
     .          w(os),w(osd),e,poti(i1),poti(j1),pot0(i0),zeta(i1),
     .          zeta(j1),rhois(i1),rhois(j1),w(ogs),w(ogd),f2)
C         ... zetxc and force from rho shifting against eps
              call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,vxc1,
     .          sum)
              call stlxc1(nlmx,nlml,nlmb,sum,w(os),-1d0,zetxc(j1))
              call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1),f1)
              if (nsp == 2) then
                call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,
     .            vxc1(1,1,2),sum)
                call stlxc1(nlmx,nlml,nlmb,sum,w(os),-1d0,zetxc(j1+nri))
                call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1+nri),f1)
              endif
              call stlset(e,lmaxl,nr2,w(oxi),w(oh),w(oy),rofi(kr),
     .          rwgt(kr),dvxc2,sum)
              call stlxc1(nlmx,nlml,nlmb,sum,w(os),1d0,zetxc(j1))
              call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1),f7)
              if (nsp == 2) then
                call stlset(e,lmaxl,nr2,w(oxi),w(oh),w(oy),rofi(kr),
     .            rwgt(kr),dvxc2(1,1,2),sum)
                call stlxc1(nlmx,nlml,nlmb,sum,w(os),1d0,zetxc(j1+nri))
                call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1+nri),f7)
              endif
C         ... Zeta and force from rho shifting against estatic pot vl
              call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,vl,
     .          sum)
              call stlxc1(nlmx,nlml,nlmb,sum,w(os),-1d0,zeta(j1))
              call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhois(j1),f2)
C         ... Force from ves shifting against rho
              call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,rl1,
     .          sum)
              call stlfor(nlmx,nlml,nlmb,sum,w(ogs),poti(j1),f3)
C         ... Force from vval and qmom, also deriv of sphere charge
              call stlqi1(e,r1,lmaxl,vval,qmom,xval,xmom,x0)
              call stlfor(nlmx,nlml,nlmb,xval,w(ogs),rhois(j1),f4)
              call stlfor(nlmx,nlml,nlmb,xmom,w(ogs),poti(j1),f5)
              call stlqi2(nlmx,nlml,nlmb,x0,w(ogs),rhois(j1),f6)

C       --- Special case : classical multipoles ---
            else

C         ... 2-center electrostatic contribution to force
C              call stles3(nlmx,lmaxib,nlmb,e,pot0(i0),rhois(j1),w(ogs),
C     .          f2)
C         ... Force and energy from qmom
              call stlqi3(lmaxp,qmom,xmom)
              call stlfor(nlmx,nlml,nlmb,xmom,w(ogs),poti(j1),f5)
              call stlepm(nlmx,nlml,nlmb,xmom,w(os),poti(j1),upm)
            endif
            j1 = j1+nlmb
          enddo
          call rlse(ohl)
        endif

C   --- Functions with energy zero ---
   50   continue
        j0 = ioffv0(jb)+1
        nlmb = (lmax+1)**2
        lmxh = lmaxib+lmax+1
        nlmh = (lmxh+1)**2
        ojj = ojkl(lmax)
C       call nstrug(0d0,dr,nlmx,nlmb,nlmp,0   ,w(oikl),w(ojj),w(oip),
C    .    w(ocf),cy,w(os),w(ogs))
        call defrr(ohl,nlmh)
        call rstr0(1,lmax,0d0,nlmh,1,dr(1),dr(2),dr(3),
     .    lmaxib+1,1,w(ohl),w)
        call hstrug(0d0,nlmx,nlmb,nlmp,0,0,w(oikl),w(ojj),w(oip),
     .    w(ocf),w(ohl),w(os),w(ogs))
        call rlse(ohl)
C   ... 2-center electrostatic terms for spheres
        if (ib <= nbas) then
          call stles0(r1,nxi(is),lxi(1,is),exi(1,is),nlmx,nlmb,w(os),
     .      pot0(j0),zeta(i1),rhois(i1),w(ogs),f2)
C     ... Force from ves shifting against rho
          call stlset(0d0,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,rl1,
     .      sum)
          call stlfor(nlmx,nlml,nlmb,sum,w(ogs),pot0(j0),f3)
C     ... Force from vval and qmom
          call stlqi1(0d0,r1,lmaxl,vval,qmom,xval,xmom,x0)
          call stlfor(nlmx,nlml,nlmb,xmom,w(ogs),pot0(j0),f5)
C   ... 2-center electrostatic terms for classical multipoles
        else
          call stlqi3(lmaxp,qmom,xmom)
          call stlfor(nlmx,nlml,nlmb,xmom,w(ogs),pot0(j0),f5)
          call stlepm(nlmx,nlml,nlmb,xmom,w(os),pot0(j0),upm)
        endif

C   --- Add together force contributions ---
        if (ipr >= 60) write(stdo,809) ib,jb,f1,f2,f3,f4,f5,f6
  809   format(' ib,jb=',2i3,'   f1=',3f12.6/13x,'   f2=',3f12.6/13x,
     .    '   f3=',3f12.6/13x,'   f4=',3f12.6/13x,'   f5=',3f12.6
     .    /13x,'   f6=',3f12.6)
        if (ipr >= 60) write(stdo,899) f7
  899   format(13x,'   f7=',3f12.6)
        do  m = 1, 3
          for1(m,ib) = for1(m,ib) - f1(m) + f7(m)
          for1(m,jb) = for1(m,jb) + f1(m) - f7(m)
          fes = 0.5d0*(f2(m)+f3(m)+ f4(m) - f5(m))
          for2(m,ib) = for2(m,ib) - fes
          for2(m,jb) = for2(m,jb) + fes
          for4(m,ib) = for4(m,ib) - f6(m)
          for4(m,jb) = for4(m,jb) + f6(m)
        enddo
   20   continue
      enddo

      call rlse(oxi)
      end
      subroutine stlset(e,lmaxl,nr1,xi,h,y,rofi,rwgt,vl,sum)
C- Set up integrals bessels*potential for one energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :Energy of interstitial function, bessel as one-center exp.
Ci   lmaxl :l cutoff for local density
Ci   nr1   :number of radial mesh points at site ib
Ci   xi    :work array to hold one-center expansion on radial mesh
Ci   h     :work array of length nr1
Ci   y     :work array of length nr1
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   vl    :potential expanded to lmaxl, on radial mesh
Co Outputs
Co   sum   :integral potential * smoothed function
Cu Updates
Cu   22 Jun 01
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmaxl,nr1
      double precision e,vl(nr1,1),rofi(1),xi(nr1,0:1),sum(1),
     .  y(1),h(1),rwgt(1)
C ... Local parameters
      integer i,ilm,l,m
      double precision sam

      call ropbes(rofi,e,lmaxl,y,h,xi,nr1,1)
      do  i = 1, nr1
        h(i) = rofi(i)*rwgt(i)
      enddo
      ilm = 0
      do  l = 0, lmaxl
        do  i = 1, nr1
          h(i) = h(i)*rofi(i)
        enddo
        do  m = -l,l
          ilm = ilm+1
          sam = 0d0
          do  i = 1, nr1
            sam = sam + h(i)*xi(i,l)*vl(i,ilm)
          enddo
          sum(ilm) = sam
        enddo
      enddo
      end
      subroutine stlxc1(nlmx,nlml,nlmb,sum,s,fac,zetj)
C- Add to zetj
      implicit none
C ... Passed parameters
      integer nlmx,nlml,nlmb
      double precision fac,sum(1),s(nlmx,nlmb),zetj(1)
C ... Local parameters
      integer ilm,jlm
      do  15  jlm = 1, nlmb
      do  15  ilm = 1, nlml
   15 zetj(jlm) = zetj(jlm) + fac*s(ilm,jlm)*sum(ilm)
      end
      subroutine stlfor(nlmx,nlml,nlmb,sum,gs,rhoi,f)
C- Add contribution to force from a shift in (rho,pot)_j x (pot,rho)_i
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlmx  :leading dimension of gs
Ci   nlml  :(density l-cutoff+1)**2
Ci   nlmb  :l cutoff of (pot or rho) at site i
Ci   sum   :either potential or density at site i
Ci   gs    :gradient of strux connecting site j to site i
Ci   rhoi  :coefficients to density at site j
Co Outputs
Co   f     :contribution to force added
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmx,nlml,nlmb
      double precision sum(1),rhoi(1),f(3),gs(nlmx,nlmb,3)
C ... Local parameters
      integer ilm,jlm
      do  jlm = 1, nlmb
      do  ilm = 1, nlml
        f(1) = f(1) + gs(ilm,jlm,1)*rhoi(jlm)*sum(ilm)
        f(2) = f(2) + gs(ilm,jlm,2)*rhoi(jlm)*sum(ilm)
        f(3) = f(3) + gs(ilm,jlm,3)*rhoi(jlm)*sum(ilm)
      enddo
      enddo
      end
      subroutine stlepm(nlmx,nlml,nlmb,qmom,s,potj,u)
C- Add contribution to energy from pot_j x q_i
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlmx  :leading dimension of s
Ci   nlml  :(density l-cutoff+1)**2
Ci   nlmb  :l cutoff of (pot or rho) at site i
Ci   qmom  :moment at site i
Ci   s     :strux connecting site j to site i
Ci   potj  :coefficients to density at site j
Co Outputs
Co   u     :contribution to energy added
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmx,nlml,nlmb
      double precision qmom(1),potj(1),u,s(nlmx,nlmb)
C ... Local parameters
      integer ilm,jlm
      do  jlm = 1, nlmb
      do  ilm = 1, nlml
        u = u + s(ilm,jlm)*potj(jlm)*qmom(ilm)
      enddo
      enddo
      end

      subroutine stlqi2(nlmx,nlml,nlmb,x0,gs,rhoi,f)
C- Only for derivative of instl charge
      implicit none
C ... Passed parameters
      integer nlmb,nlml,nlmx
      double precision x0,rhoi(1),f(3),gs(nlmx,nlmb,3)
C ... Local parameters
      integer jlm

      do  15  jlm = 1, nlmb
      f(1) = f(1) - gs(1,jlm,1)*rhoi(jlm)*x0
      f(2) = f(2) - gs(1,jlm,2)*rhoi(jlm)*x0
      f(3) = f(3) - gs(1,jlm,3)*rhoi(jlm)*x0
   15 continue
      end
      subroutine stlqi1(e,r1,lmaxl,vval,qmom,xval,xmom,x0)
C- Setup for force terms from qmom and vval
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :energy of hankel at remote site, expanded at ib as bessel
Ci   r1    :augmentation radius
Ci   lmaxl :l cutoff for local density
Ci   vval  :value of electrostatic potential at site ib
Ci   qmom  :multipole moments of true density at site ib
Co Outputs
Co   xval  :vval * (mpol moment of bessel) / r**l
Co   xmom  :qmom * (val of bessel) / r**l
Co   x0    :integral of bessel, l=0, inside site ib
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmaxl
      double precision e,r1,x0,xval(1),xmom(1),vval(1),qmom(1)
C ... Local parameters
      integer ilm,l,m
      double precision bmom,bval,phi(0:20),psi(0:20)

      call bessl(e*r1*r1,lmaxl+1,phi,psi)
      ilm = 0
      do  l = 0, lmaxl
C       bval = val of bessel at r1; bmom = multipole moment of bessel
        bval = phi(l)*r1**l
        bmom = phi(l+1)*r1**(2*l+3)
        if (l == 0) x0 = bmom*3.544907702d0
        do  m = -l, l
          ilm = ilm+1
          xval(ilm) = vval(ilm)*(bmom/r1**l)
          xmom(ilm) = qmom(ilm)*(bval/r1**l)
        enddo
      enddo
      end

      subroutine stlqi3(lmaxp,qmom,xmom)
C- Setup for force terms from qmom, for classical multipoles
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmaxp :l cutoff for classical multipole
Ci   qmom  :multipole moments of true density at site ib
Co Outputs
Co   xmom  :qmom * (val of bessel) / r**l
Cr Remarks
Cr   stlqi3 is analog of stlqi1 for point multipoles, r->0.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmaxp
      double precision qmom(1),xmom(1)
C ... Local parameters
      integer ilm,l,m
      double precision df

      ilm = 0
      df = 1
      do  l = 0, lmaxp
        df = df*(2*l+1)
        do  m = -l, l
          ilm = ilm+1
          xmom(ilm) = qmom(ilm)/df
        enddo
      enddo
      end

      subroutine stles2(r1,nxi,lxi,exi,nlmx,nlmb,s,sd,ej,poti,potj,
     .  pot0i,zeti,zetj,rhoi,rhoj,gs,gd,f2)
C- 2-center terms in sphere contr to zeta, one pair of sites, one energy ej
C ----------------------------------------------------------------------
Ci Inputs
Ci   r1    :augmentation radius at site ib
Ci   nxi   :number of functions per l for c.d. basis at ib
Ci   lxi   :l cutoffs for each function in c.d. basis at ib
Ci   exi   :Hankel energies for each function in c.d. basis at ib
Ci   nlmx  :leading dimension of s,sd,gs,gd
Ci   nlmb  :(lmax+1)**2 for site j
Ci   s     :structure constants connecting sites i (energies exi) and j(ej)
Ci   sd    :energy derivative of s
Ci   ej    :energy at site j
Ci   poti  :estat coeffs for site i
Ci   potj  :estat coeffs for site j
Ci   pot0i :coefficients to estat potential, homogeneous contribution
Ci   rhoi  :density coeffs for site i
Ci   rhoj  :density coeffs for site j
Ci   gs    :gradient of s
Ci   gd    :gradient of sd
Co Outputs
Co   zeti  :integrals of density with estat potential, site i:
Co         :contribution to zeti from (j,ej) contr. to potj subtracted
Co   zetj  :integrals of density with estat potential, site j
Co         :contribution to zetj from (j,ej) contr. to potj subtracted
Co   f2    :integral ( ves * drho) added into f2
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lxi(1),nlmb,nlmx,nxi
      double precision  ej,r1
      double precision poti(1),potj(1),s(nlmx,nlmb),zeti(1),zetj(1),
     .  exi(1),sd(nlmx,nlmb),pot0i(1),rhoi(1),rhoj(1),gs(nlmx,nlmb,3),
     .  gd(nlmx,nlmb,3),f2(3)
C ... Local parameters
      integer i,ie,ilm,jlm,l,lmax,lx,m
      double precision ex,fac,sx,wk,xxx,
     .  fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20)

      i = 0
      lmax = -1
      do  10  ie = 1, nxi
      lx = lxi(ie)
      lmax = max0(lmax,lx)
      ex = exi(ie)
      call wronkj(ex,ej,r1,lx,fkk,fkj,fjk,fjj)
C ... For ex /= ej
      if (dabs(ej-ex) > 1d-8) then
        ilm = 0
        do  l = 0, lx
          fac = fkj(l) + 1d0/(ej-ex)
          do  m = -l, l
            i = i+1
            ilm = ilm+1
            do  jlm = 1, nlmb
              xxx = -fac*(rhoj(jlm)*poti(i) + rhoi(i)*potj(jlm))
              f2(1) = f2(1) + xxx*gs(ilm,jlm,1)
              f2(2) = f2(2) + xxx*gs(ilm,jlm,2)
              f2(3) = f2(3) + xxx*gs(ilm,jlm,3)
            enddo
            do  jlm = 1, nlmb
              zetj(jlm) = zetj(jlm) + s(ilm,jlm)*fac*poti(i)
              zeti(i  ) = zeti(i  ) + s(ilm,jlm)*fac*potj(jlm)
            enddo
          enddo
        enddo
C ... For ej == ex
      else
        ilm = 0
        do  l = 0, lx
        fac = fkj(l)
        do  m = -l, l
          i = i+1
          ilm = ilm+1
          do  jlm = 1, nlmb
            xxx = -(rhoj(jlm)*poti(i) + rhoi(i)*potj(jlm))
            f2(1) = f2(1)+xxx*(fac*gs(ilm,jlm,1) + 0.5d0*gd(ilm,jlm,1))
            f2(2) = f2(2)+xxx*(fac*gs(ilm,jlm,2) + 0.5d0*gd(ilm,jlm,2))
            f2(3) = f2(3)+xxx*(fac*gs(ilm,jlm,3) + 0.5d0*gd(ilm,jlm,3))
          enddo
C   ... This code causes problems when ib=jb.
C          do  12  jlm = 1, nlmb
C          sx=s(ilm,jlm)*fac+sd(ilm,jlm)*0.5d0
C          zetj(jlm)=zetj(jlm)+sx*poti(i)
C  12      zeti(i  )=zeti(i  )+sx*potj(jlm)
C     ... Use this instead
          wk = 0d0
          do  jlm = 1, nlmb
            sx = s(ilm,jlm)*fac + sd(ilm,jlm)*0.5d0
            zetj(jlm) = zetj(jlm) + sx*poti(i)
            wk = wk + sx*potj(jlm)
          enddo
          zeti(i) = zeti(i) + wk
        enddo
        enddo
      endif
   10 continue

C ... Terms from v0(ilm) rhoj(jlm)
      call wronkj(0d0,ej,r1,lmax,fkk,fkj,fjk,fjj)
      ilm = 0
      do  l = 0, lmax
        fac = fkj(l) + 1d0/ej
        do  m = -l, l
          ilm = ilm+1
          do  jlm = 1, nlmb
            xxx = -fac*rhoj(jlm)*pot0i(ilm)
            f2(1) = f2(1) + xxx*gs(ilm,jlm,1)
            f2(2) = f2(2) + xxx*gs(ilm,jlm,2)
            f2(3) = f2(3) + xxx*gs(ilm,jlm,3)
          enddo
          do    jlm = 1, nlmb
            zetj(jlm) = zetj(jlm) + s(ilm,jlm)*fac*pot0i(ilm)
          enddo
        enddo
      enddo
      end

      subroutine stles3(nlmx,lmax,nlmb,ej,pot0i,rhoj,gs,f2)
C- 2-center terms in sphere contr to force, classical multipoles
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlmx  :leading dimension of gs
Ci   lmax  :density l cutoff for augmentation site,
Ci         :corresponding to pot0
Ci   ej    :energy at site j
Ci   pot0i :coefficients to estat potential, homogeneous contribution
Ci   rhoj  :density coeffs for site j
Ci   gs    :gradient of strux sites i, e=0 to site j
Co Outputs
Co   f2    :integral ( ves * drho) added into f2
Cr Remarks
Cr   Contributions to force from v0(ilm) rhoj(jlm), where v0 is centered
Cr   at a classical multipole site and rhoj centered at a sphere.
Cr   Adapted from stles2.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmb,lmax,nlmx
      double precision  ej
      double precision pot0i(1),rhoj(1),gs(nlmx,nlmb,3),f2(3)
C ... Local parameters
      integer ilm,jlm,l,m
      double precision fac,xxx,r1,
     .  fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20)

      r1 = .7d0
      call wronkj(0d0,ej,r1,lmax,fkk,fkj,fjk,fjj)
      ilm = 0
      do  l = 0, lmax
        fac = fkj(l) + 1d0/ej
        do  m = -l, l
          ilm = ilm+1
          do  jlm = 1, nlmb
            xxx = -fac*rhoj(jlm)*pot0i(ilm)
            f2(1) = f2(1) + xxx*gs(ilm,jlm,1)
            f2(2) = f2(2) + xxx*gs(ilm,jlm,2)
            f2(3) = f2(3) + xxx*gs(ilm,jlm,3)
          enddo
        enddo
      enddo
      end

      subroutine stles1(r1,nxi,lxi,exi,poti,pot0,zeti)
C- Add 1-center terms to zeta
C ----------------------------------------------------------------------
Ci Inputs
Ci   r1    :augmentation radius
Ci   nxi   :number of functions per l for c.d. basis, by species
Ci   lxi   :l cutoffs for each function in c.d. basis, by species
Ci   exi   :Hankel energies for each function in c.d. basis, by species
Ci   poti  :coefficients to estat potential in atom-centered density basis
Ci   pot0  :coefficients to estat potential, homogeneous contribution
Co Outputs
Co   zeti  :The on-site contribution to es pot * zeta is subtracted
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nxi,lxi(1)
      double precision r1,poti(1),zeti(1),exi(1),pot0(1)
C ... Local parameters
      integer i0,ie,ilm,j0,je,l,li,lj,ll,lmax,lx
      double precision fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20)

      i0 = 0
      do  ie = 1, nxi
        li = lxi(ie)
        j0 = 0
        lmax = -1
        do  je = 1, nxi
          lj = lxi(je)
          lmax = max0(lmax,lj)
          lx = min0(li,lj)
          call wronkj(exi(ie),exi(je),r1,lx,fkk,fkj,fjk,fjj)
          do  ilm = 1, (lx+1)**2
            l = ll(ilm)
            zeti(ilm+i0) = zeti(ilm+i0) + poti(ilm+j0)*fkk(l)
          enddo
          j0 = j0+(lj+1)**2
        enddo
        lx = min0(li,lmax)
        call wronkj(exi(ie),0d0,r1,lx,fkk,fkj,fjk,fjj)
        do  ilm = 1,(lx+1)**2
          l = ll(ilm)
          zeti(ilm+i0) = zeti(ilm+i0) + pot0(ilm)*fkk(l)
        enddo
        i0 = i0+(li+1)**2
      enddo
      end
      subroutine stles0(r1,nxi,lxi,exi,nlmx,nlmb,s,pot0j,zeti,rhoi,gs,
     .  f2)
C- 2-center terms in sphere contr to zeta, one pair of sites, pot0j
C ----------------------------------------------------------------------
Ci Inputs
Ci   r1    :augmentation radius
Ci   nxi   :number of functions per l for c.d. basis, by species
Ci   lxi   :l cutoffs for each function in c.d. basis, by species
Ci   exi   :Hankel energies for each function in c.d. basis, by species
Ci   nlmx  :leading dimension of s,sd,gs,gd
Ci   nlmb  :(lmax+1)**2 for site j
Ci   s     :strux connecting sites i (energies exi) and j(e=0)
Ci   pot0j  :estat coeffs for site j, energy zero
Ci   rhoi  :density coeffs for site i
Ci   gs    :gradient of s
Co Outputs
Co   zeti  :integrals of density with estat potential, site i:
Co         :contribution to zeti from (j,e=0) contr. to potj subtracted
Co   f2    :integral ( ves * drho) added into f2
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmb,nlmx,nxi,lxi(1)
      double precision f2(3),r1,pot0j(1),s(nlmx,nlmb),zeti(1),
     .  gs(nlmx,nlmb,3),rhoi(1),exi(1)
C ... Local parameters
      integer i,ie,ilm,jlm,l,lx,m
      double precision ex,fac,xxx
      double precision fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20)

      i = 0
      do  ie = 1, nxi
        lx = lxi(ie)
        ex = exi(ie)
        call wronkj(ex,0d0,r1,lx,fkk,fkj,fjk,fjj)
        ilm = 0
        do  l = 0, lx
          fac = fkj(l) - 1d0/ex
          do  m = -l, l
            i = i+1
            ilm = ilm+1
            do    jlm = 1, nlmb
              xxx = -fac*rhoi(i)*pot0j(jlm)
              f2(1) = f2(1) + xxx*gs(ilm,jlm,1)
              f2(2) = f2(2) + xxx*gs(ilm,jlm,2)
              f2(3) = f2(3) + xxx*gs(ilm,jlm,3)
            enddo
            do  jlm = 1, nlmb
              zeti(i) = zeti(i) + s(ilm,jlm)*fac*pot0j(jlm)
            enddo
          enddo
        enddo
      enddo
      end
