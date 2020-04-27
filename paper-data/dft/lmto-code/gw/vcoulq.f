      subroutine vcoulq(s_lat,lgrad,q,ntpb,npwmb,ntmb,nbas,nat,
     .  lcutmx,npb,mxnpbl,qlat,strx,rojp,rojb,sgbb,sgpb,fouvb,E,
     .  vcoul)
C- Coulomb matrix (Yukawa potential) for one q.
C ----------------------------------------------------------------------
Ci Inputs
Ci   q     :calculate vcoulq(q)
Ci   ntpb  :Total total number of local product basis functions within MTs
Ci         :Formerly called nbloch in old GW code
Ci   npwmb :number of IPW's for this q
Ci   ntmb  :largest dimension of mixed product basis functions = ntpb + npmbx
Ci         :Formerly called nblochpmx in old GW code
Ci         :Maybe not needed?
Ci   nbas  :size of basis
Ci   nat   :Number of atoms with augmentation spheres, needed for dimensioning
Ci   lcutmx:maximum l-cutoff for product basis functions
Ci         :Dimensions npb,rprodx,rojb,sgbb
Ci         :Called lxx in old GW code.
Ci   npb   :indexing for rprodx.
Ci         :Functions npb(l,iat):npb(l+1,iat) are the family of radial B functions
Ci         :stored in rprodx for quantum number l and site iat.
Ci   mxnpbl:maximum number of radial B functions; dimensions rojb,sgbb
Ci         :Called nxx in old GW code.
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   vol   :cell volume
Ci   strx  :Structure constants (strxq)
Ci   rojb  :Onsite integral rojb(BRlμ) = 1/(2l+1)!! int dr r*Jl(E,r) r*BRlμ(r)
Ci         :Needed to evaluate integrals <B^k_RLμ(r)|V_l(r,r')|B^k_R'L'μ'(r')>
Ci         :when (R,T) and (R',T') are different since v expanded in J_L(E,r)
Ci         :Here V_l(r,r') = (r<)J_l(E,r<) (r>)H_l(E,r>)
Ci   sgbb  :Onsite integrals <B|v(onsite)|B> corresponding to rojb when (R,T) = (R',T')
Ci         :sgbb = sigma integral, given by Eq 43 in ecalj manual (June 19, 2015)
Ci         :Follows from the 1-center expansion of the Yukawa potential
Ci         :sgbb(B_Rlμ, B_Rlν)
Ci         :  = 4*pi/(2l+1) * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r B_Rlμ(r)) (r' B_Rlν(r'))
Ci         :  = 4*pi/(2l+1) int dr dr' (r B_Rlμ(r)) V_l(r,r') (r' B_Rlν(r'))
Ci         :where V_l(r,r') = (r<)J_l(E,r<) (r>)H_l(E,r>)
Ci   rojp  :Onsite integral <j_L(E,r) | P(q+G)_L > where
Ci         :  |P(q+G)_L> = projection of exp(i (q+G) r) to channnel L
Ci         :  |j_L(E=0,r)> = r^l/(2l+1)!! Y_L [check]
Co         : Needed to evaluate <P_G|v|B> since v expanded in J_L(E,r)
Ci   sgpb  :Onsite integrals <J|v(onsite)|B> corresponding to rojp when (R,T) = (R',T')
Ci         :Here J=J(G^2) --- 1-center expansion of plane waves
Ci         :  sgpb(J(G^2), B_Rlν) =
Ci         :  4*pi/(2l+1) (pjyl)+ * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r J_L(r,G^2)) (r' B_Rlν(r'))
Ci   fouvb :Fourier integral of PW and product basis inside augmentation sphere
Ci         :  4*pi/(G^2-E) * integral(J*B) ... note E is negative
Ci   E     :Convergence parameter (bare coulomb -> Yukawa with decay E)
Co Outputs
Co   vcoul :coulomb interaction, in Hartrees
Co         :To convert to Ry units, e^2=2 => vcoul (Ry) = 2*vcoul (Hartree)
Cl Local variables
Cl   lcuta :this may eventually be site based.  Right now it is lcutmx for all sites
Cl   rkpr  :r*Jl(E,r),  Jbar definition (see mkrjk)
Cl   rkmr  :r*Hl(E,r),  Hbar definition (see mkrjk)
Cl   rofi  :radial mesh points
Cr Remarks
Cr
Cu Updates
Cu   05 May 15
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_lat)::   s_lat
C ... Passed parameters
      integer :: lgrad,ntpb,ntmb,npwmb,nbas,nat,lcutmx,mxnpbl,npb(0:lcutmx*(nat+1))
      real(8) :: alat,qlat(3,3),vol,q(3),E
      complex(8) :: strx((lcutmx+1)**2,nbas,(lcutmx+1)**2,nbas)
C     rho-type onsite integral
      real(8)    :: rojb(mxnpbl,0:lcutmx,0:lgrad,nat)
      complex(8) :: rojp(npwmb,(lcutmx+1)**2,nat)
C     sigma-type onsite integral
      real(8)    :: sgbb(mxnpbl,mxnpbl,0:lcutmx,nat)
      complex(8) :: sgpb(npwmb,mxnpbl,(lcutmx+1)**2,nat)
C     Fourier
      complex(8) :: fouvb(npwmb,mxnpbl,(lcutmx+1)**2,nat)
C     Matrix elements
      complex(8) :: vcoul(ntmb,ntmb)

C ... Dynamically allocated arrays
      real(8), allocatable :: qpg2(:)

C ... Local parameters
      logical, parameter :: debug=.false.
      integer iat,ig,l,m,n,i1,i2,nlmbas,ibl1,ibl2,m1,m2,n1,n2,ilm,jobBes,lcuta
      real(8):: pi,fpi,tpiba,qpg(3)
      complex(8), parameter :: img=(0d0,1d0)
      integer :: ipb,ilmb,nf1,nf2
      complex(8), allocatable :: rojblm(:,:), tmp(:,:)
      complex(8), parameter :: zone = 1.0_8, znull = 0.0_8
!MvS OMP related variables
!     integer omppid,idalloc,nthreads,tid

      procedure(integer) :: idalloc,allocvb
      procedure(real(8)) :: ddot

C#ifdefC ptest
C      integer nblochnpwmb
C      complex(8),allocatable :: matp(:),matp2(:)
C#endif

C     double precision delwc,swalltime

      call tcn('vcoulq')

      pi     = 4d0*datan(1d0)
      fpi    = 4*pi
      alat   = s_lat%alat
      vol    = s_lat%vol
      tpiba  = 2*pi/alat
      jobBes = 0  ! mod(GWversion()/1000,10)
      lcuta  = lcutmx   ! If lcuta is site-dependent, need to rework rojp below

      nlmbas = nat * (lcutmx+1)**2

C ... Assemble rojblm : robp with radial functions sequential and resolved by m
C     Is this the best way to represent?  This array is really sparse!
      allocate(rojblm(ntpb,nlmbas))
      rojblm = znull
      ipb = 0; ilmb = 0
      do  iat = 1, nat
        do  l  = 0, lcuta
          nf1 = npb(l+1+(lcutmx+1)*(iat-1)) - npb(l+(lcutmx+1)*(iat-1))
          do  n  = 1, nf1
            ilm = l*l
            do  m  = 1, 2*l+1
              ipb  = ipb + 1
              ilm  = ilm+1
              rojblm(ipb,ilmb+ilm) = rojb(n,l,0,iat)
            enddo
          enddo
        enddo
        ilmb = ilmb + (lcuta+1)**2
      enddo
      if (ipb /= ntpb) call rx('vcoulq: product basis mismatch')

C --- <B|v|B> ---
      call tcn('BvB')
C ... Offsite term
      allocate(tmp(ntpb,nlmbas))
      call zgemm('n','n',ntpb,nlmbas,nlmbas,zone,rojblm,ntpb,strx,nlmbas,znull,tmp,ntpb)
      call zgemm('n','t',ntpb,ntpb,nlmbas,zone,tmp,ntpb,rojblm,ntpb,znull,vcoul,ntmb)
      deallocate(tmp)

C ... Onsite term (diagonal in m), in full product basis order
      ibl1 = 1
      i1 = 1
      do  iat = 1, nat
        do  l = 0, lcuta
          nf1 = npb(l+1+(lcutmx+1)*(iat-1)) - npb(l+(lcutmx+1)*(iat-1))
          do  n1 = 1, nf1
            do  m1 = 1, 2*l+1
              ibl2 = i1
              nf2 = npb(l+1+(lcutmx+1)*(iat-1)) - npb(l+(lcutmx+1)*(iat-1))
              do  n2 = 1, nf2
                do  m2 = 1, 2*l+1
                  if (m1 == m2) vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l,iat)
                  ibl2 = ibl2 + 1
                end do
              end do
              ibl1 = ibl1 + 1
            end do
          end do
          i1 = ibl2
        end do
      end do

C#ifdefC DEBUG
C      call info2(0,0,0,' sumcheck vcoulq(bb)  q=%s,%3d %s %2;18,6D',q,sum(vcoul(1:ntpb,1:ntpb)))
C#ifdefC DEBUG2
C      call zprm0('(1p9e22.12)'); call zprm('vcoul(bb)',2,vcoul,ntmb,ntpb,ntpb)
C#endifC
C#endif

      call tcx('BvB')

      if (npwmb == 0) goto 99
C
C --- <P_G|v|B> ---
      call tcn('PvB')

      allocate(tmp(nlmbas,npwmb))
      call zgemm('t','c',nlmbas,npwmb,nlmbas,zone,strx,nlmbas,rojp,npwmb,znull,tmp,nlmbas)
      call zgemm('t','t',npwmb,ntpb,nlmbas,-zone,tmp,nlmbas,rojblm,ntpb,znull,vcoul(ntpb+1,1),ntmb)

      ibl1 = 1
      do  iat = 1, nat
        do  l = 0, lcuta
          nf1 = npb(l+1+(lcutmx+1)*(iat-1)) - npb(l+(lcutmx+1)*(iat-1))
          do  n  = 1, nf1
            do  m = 1, 2*l+1
              ilm = l*l+m
              call zaxpy(npwmb,zone,fouvb(1,n,ilm,iat),1,vcoul(ntpb+1,ibl1),1)
              call zaxpy(npwmb,-zone,sgpb(1,n,ilm,iat),1,vcoul(ntpb+1,ibl1),1)
              ibl1 = ibl1 + 1
            end do
          end do
        end do
      end do

C#ifdefC DEBUG
C      call info2(0,0,0,' sumcheck vcoulq(pb)  q=%s,%3d %s %2;18,6D',q,sum(vcoul(ntpb+1:ntpb+npwmb,1:ntpb)))
C#ifdefC DEBUG2
C      call zprm0('(1p9e22.12)'); call zprm('vcoul(bb+pb)',2,vcoul(ntpb+1,1),ntmb,npwmb,ntpb)
C#endifC
C#endif

      call tcx('PvB')

C --- <P_G|v|P_G> ---
      call tcn('PvP')

      allocate(qpg2(npwmb))
      do  ig = 1, npwmb
        qpg(1:3) = tpiba * (q(1:3) + matmul(qlat,s_lat%igvcc(1:3,ig)))
        qpg2(ig) = ddot(3,qpg,1,qpg,1)
      enddo

C ... Make (strx)^T * rojp^+ rojp^T.
C     Note: requires rojp to have contiguous elements => lcuta must be lcutmx for all sites
      call zgemm('t','t',npwmb,npwmb,nlmbas,zone,tmp,nlmbas,rojp,npwmb,zone,vcoul(ntpb+1,ntpb+1),ntmb)
      deallocate(tmp)

C ... Diagonal term in interstitial vcoulq
      do  ig = 1, npwmb
        vcoul(ntpb+ig,ntpb+ig) = vcoul(ntpb+ig,ntpb+ig) + fpi*vol/(qpg2(ig) - E) !E must be negative
      enddo

      call tcx('PvP')

C ... Re-entry point if npwmb=0
   99 continue
      deallocate(rojblm)

C --- Right-upper triangle of vcoul ---
      do  i1 = 1, ntpb+npwmb
        do  i2 = 1, i1-1
          vcoul(i2,i1) = dconjg(vcoul(i1,i2))
        enddo
        vcoul(i1,i1) = dble(vcoul(i1,i1))
      enddo

      call info0(60,0,0,'%9fDiagonal Vcoul')

      do  i1 = 1, ntpb+npwmb
        if (i1==11) call info0(60,0,0,' ... ')
        if(i1>10 .and. i1<ntpb+npwmb .and. mod(i1,max(10,(ntpb+npwmb)/20))/=0) cycle
        call info2(60,0,0,' %,5i %2:2;9F',i1,vcoul(i1,i1))
      enddo

C#ifdefC DEBUG
C      call info2(0,0,0,' sumcheck vcoulq/1000 q=%s,%3d %s %2;18,6D',q,sum(vcoul(1:ntpb+npwmb,1:ntpb+npwmb))/1000)
C#ifdefC DEBUG2
C      call zprm0('(1p9e22.12)'); call zprm('vcoul',2,vcoul,ntmb,ntpb+npwmb,ntpb+npwmb)
C#endifC
C#endif

      call tcx('vcoulq')
      end
      subroutine onsitecoul(s_lat,s_site,s_spec,mode,lgrad,nbas,lcutmx,mxnpbl,npb,E,nrmx,rprodx,
     .  q,ntpb,npwmb,ntmb,rojb,sgbb,sgbc,rojp,sgpb,fouvb,vcoul)
C- Make onsite Coulomb integrals in all augmentation spheres
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa a nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt
Cio    Passed to:  *
Ci Inputs
Ci   mode  :controls what is made
Ci         :1s digit
Ci         :0 make none of rojb, sgbb, sgbc
Ci         :1 make rojb
Ci         :2 make sgbb
Ci         :4 make sgbc (requires lgrad also be set)
Ci         :Any combination is allowed
Ci         :10s digit
Ci         :1 make rojp,sgpb,fouvb, on-site of PW part of vcoul
Ci   lgrad :0 make coulomb integrals
Ci         :1 Also make integrals needed for gradients of coulomb interaction (not finished)
Ci   nbas  :size of basis
Ci   lcutmx:maximum l-cutoff for product basis functions, typically 2*(nl-1)
Ci         :Dimensions npb,rprodx,rojb,sgbb
Ci         :Called lxx in old GW code.
Ci   lcuta :site-dependent lcut for integrals rojb,sgbb,...
Ci         :Called lx in old GW code.  For now, lcuta = lcutmx
Ci   mxnpbl:maximum number of radial B functions; dimensions rojb,sgbb
Ci         :Called nxx in old GW code.
Ci   npb   :Table of offsets to radial functions in rprodx
Ci         :Functions npb(l,iat):npb(l+1,iat) are the family of radial B functions
Ci         :stored in rprodx for quantum number l and site iat.
Ci         :Formerly called nx in old GW code
Ci   E     :Use Yukawa potential for coulomb interaction, 1/sqrt(-E) is screening length
Ci   nrmx  :Maxiumum number of radial mesh points, and leading dimension of rprodx
Ci         :Formerly called nrx in old GW code
Ci   rprodx:Product basis functions * r: r*BRlμ(r)
Ci         :Normalized such that \int_0^rmax (r*BRlμ(r))^2 dr = 1
Ci   q     :Make integrals for wave number q
Ci   ntpb  :Total total number of local product basis functions within MTs
Ci         :Formerly called nbloch in old GW code
Ci   npwmb :number of G vectors for Coulomb interaction and
Ci         :number of interstitial product functions for this q
Ci         :Called ngc in old GW code.
Ci   ntmb  :largest dimension of mixed product basis functions = ntpb + npmbx
Ci         :Used here solely as leading dimension of vcoul
Ci         :Formerly called nblochpmx in old GW code
Co Outputs
Co   rojb  :Onsite integral rojb(BRlμ) = 1/(2l+1)!! int dr r*Jl(E,r) r*BRlμ(r)
Co         :Needed to evaluate integrals <B^k_RLμ(r)|V_l(r,r')|B^k_R'L'μ'(r')>
Co         :when (R,T) and (R',T') are different since v expanded in J_L(E,r)
Co         :Here V_l(r,r') = (r<)J_l(E,r<) (r>)H_l(E,r>)
Co         :Not made, but can do so (see below), if lgrad=1 (not sure if needed)
Co         :rojb(:,:,1,:) = 1/(2l+1)!! int dr r*J'_l(E,r) r*BRlμ(r)
Co   sgbb  :Onsite integrals <B|v(onsite)|B> corresponding to rojb when (R,T) = (R',T')
Co         :sgbb = sigma integral, given by Eq 43 in ecalj manual (June 19, 2015)
Co         :Follows from the 1-center expansion of the Yukawa potential
Co         :sgbb(B_Rlμ, B_Rlν)
Co         :  = 4*pi/(2l+1) * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r B_Rlμ(r)) (r' B_Rlν(r'))
Co         :  = 4*pi/(2l+1) int dr dr' (r B_Rlμ(r)) V_l(r,r') (r' B_Rlν(r'))
Co         :where V_l(r,r') = (r<)J_l(E,r<) (r>)H_l(E,r>)
Co   rojp  :Onsite integral <j_L(E,r) | P(q+G)_L > where
Co         :  | P(q+G)_L> = projection of exp(i (q+G) r) to channnel L
Co         :  | j_L(E=0,r)> = r^l/(2l+1)!! Y_L [check]
Co         : rojp(G,L) = pjyl(L,G) * W{J_L(E), J_L(G^2)}
Co         : pjyl is defined in Eq 49 in ecalj manual (June 19, 2015)
Co         : Needed to evaluate <P_G|v|B> since v expanded in J_L(E,r)
Co   sgpb  :Onsite integrals <J|v(onsite)|B> corresponding to rojp when (R,T) = (R',T')
Co         :Here J=J(G^2) --- 1-center expansion of plane waves
Co         :  sgpb(J(G^2), B_Rlν) =
Co         :  4*pi/(2l+1) (pjyl)+ * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r J_L(r,G^2)) (r' B_Rlν(r'))
Co   fouvb :Fourier integral of PW and product basis inside augmentation sphere
Co         :  4*pi/(G^2-E) * integral(J*B) ... note E is negative
Co Outputs
Co   vcoul :On-site + Fourier contribution to PW-PW part of vcoul is computed (see Remarks)
Cl Local variables
Cl   rkpr  :r*Jl(E,r),  Jbar definition (see Remarks)
Cl         :if lgrad is 1, r*J'l(E,r) is also computed
Cl   rkmr  :r*Hl(E,r),  Hbar definition (see Remarks)
Cl         :if lgrad is 1, r*H'l(E,r) is also computed
Cl   jG    :r*Jl(G^2,r), MSM Standard definition (besslr)  [called ajr in old mkjp_4]
Cl   vB    :Coulomb integral : int v(r,r') (r' * B(r'))
Cl   vjG   :Coulomb integral : int v(r,r') jG(r')
Cr Remarks
Cr
Cr *Expansion theorem for bare coulomb interaction:
Cr   exp(-lamda |r-r'|)/|r-r'| =
Cr   sum_lm  4*pi/(2l+1) (r<) J_l(E,r<) (r>) H_l(E,r>) (-1)^m Y_(l,m)(r) Y_(l,-m)(r')
Cr
Cr *Here the Jl and Hl are "bar Bessel and Hankel functions"
Cr  with the property Jl -> r^l  and Hl -> r^-l-1 for E=0,
Cr  are proportional to the radial part of customary spherical Hankels and Bessels.
Cr  See mkrjk for definition.
Cr
Cr *Expansion theorem is readily derived from the generating function for Legendre polynomials
Cr   1/sqrt(1+h^2-2*h*cos(gamma)) = sum_l h^l P_l(cos(gamma))
Cr *and the addition theorem, using the angle between r and r' for gamma
Cr   P_l(cos(gamma)) = 4*pi/(2l+1) sum_m (-1)^m Y_(l,m)(r) Y_(l,-m)(r')
Cr
Cr  See Section 18.4 of Takao's ecalj manual (June 2015)
Cr
Cr *Remarks on sgpp,fouvp, and vcoul.
Cr  For consistency, sgpp and fouvp should be returned as arrays together with
Cr  sgpb and fouvb, and vcoulq should not be generated here.
Cr
Cr  However, since sgpp and fouvp use up lots of memory, which moreover scales as N^3,
Cr  they are made locally and their contribution directly added to vcoulq.
Cr
Cr  If it is desired that the on-site contribution vcoulq be made elsewhere,
Cr  one alternative is to return jG and vjG from this routine and make it on the fly.
Cr  jG and vjG would have to be decomposed by site.
Cu Updates
Cu   20 Dec 17 Adapted from Kotani's mkjb_4, mkjp_4 in old GW code.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)

C ... Passed parameters
      integer :: mode,lgrad,nbas,nrmx,lcutmx,ntpb,npwmb,ntmb,mxnpbl,npb(0:lcutmx*(nbas+1))
!     integer lcuta(nbas)
      real(8) :: E,rprodx(nrmx,*)
      real(8) :: rojb(mxnpbl,0:lcutmx,0:lgrad,nbas) ! rho-type onsite integral
      real(8) :: sgbb(mxnpbl,mxnpbl,0:lcutmx,nbas) ! sigma-type onsite integral
      real(8) :: sgbc(mxnpbl,0:1,nbas) ! sigma-type onsite integral for nuc+core density and B
C     Specific to PW-related integrals
      real(8) :: q(3)
      complex(8) :: rojp(npwmb,(lcutmx+1)**2,nbas) ! rho-type onsite integral
      complex(8) :: sgpb(npwmb,mxnpbl,(lcutmx+1)**2,nbas)
      complex(8) :: fouvb(npwmb,mxnpbl,(lcutmx+1)**2,nbas)
      complex(8) :: vcoul(ntmb,ntmb)
!     complex(8) :: sgpp(npwmb,npwmb,nat)

C ... Dynamically allocated local arrays
      real(8), allocatable :: rkpr(:,:,:),rkmr(:,:,:)
      real(8), allocatable :: JBint(:),HBint(:),vB(:,:,:)
      real(8), allocatable :: qpg(:,:),yl(:,:),qpg2(:)
      complex(8),allocatable :: pjyl(:,:),phase(:)
C     Help arrays to make fast Bessel functions
      real(8), allocatable :: xi(:,:),y(:),h(:),jG(:,:,:),vjG(:,:,:),sig(:,:),wjj(:,:,:),rhoc(:)

C ... Local parameters
      integer, parameter :: jobBes=0
      real(8), parameter :: rgnuc = .01d0  ! For now.
      integer lcuta(nbas)  ! for now just copy from lcutmx
      integer :: l,n1,n2,ib,iat,is,offr,nf,intopt,ir,nr,mode0,mode1,ig,ig2,ilm,m,nsp
      double precision fac,fpi,pi,a,errmx,tpiba,rmax,r2s,wv12,wv21,qtrue
      double precision rofi(nrmx),rwgt(nrmx),wjjE(0:lcutmx),g0l(0:1)
      complex(8), parameter :: img=(0d0,1d0)
      complex(8) :: p12
      procedure(integer) :: nglob
      procedure(real(8)) :: ddot,dot3,dlength
C     Debugging
      double precision hs(0:1)

      call tcn('onsite-vc')

C --- Initialize ---
      pi  = 4d0*datan(1d0)
      fpi = 4*pi
      tpiba = 2*pi/s_lat%alat
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      call sanrg(.true.,lgrad,0,1,'onsitecoul','lgrad')

      if (mod(mode0,2) > 0) call dpzero(rojb,size(rojb))
      if (mod(mode0/2,2) > 0) call dpzero(sgbb,size(sgbb))
      if (mod(mode0/4,2) > 0 .and. lgrad > 0) call dpzero(sgbc,size(sgbc))

      allocate(rkpr(nrmx,0:lcutmx,0:lgrad),rkmr(nrmx,0:lcutmx,0:lgrad)) ! Barred H,J and radial derivatives
!     allocate(JBint(nrmx),Hbint(nrmx),vB(nrmx,mxnpbl,0:lgrad))
      allocate(JBint(nrmx),Hbint(nrmx),vB(nrmx,mxnpbl,0:2)) ! 2 => extra space for debugging
      if (mode1 > 0) then
        allocate(qpg(3,npwmb),yl(npwmb,(lcutmx+1)**2),qpg2(npwmb),pjyl((lcutmx+1)**2,npwmb))
        forall (ig = 1:npwmb) qpg(1:3,ig) = tpiba * (q(1:3) + matmul(s_lat%qlat,s_lat%igvcc(1:3,ig)))
        call ropyln(npwmb,qpg(1,:),qpg(2,:),qpg(3,:),lcutmx,npwmb,yl,qpg2)
        allocate(jG(nrmx,npwmb,0:lcutmx),vjG(nrmx,npwmb,0:lcutmx))
      endif

C --- For each augmentation sphere, do ---
      call info0(30,1,0,'')
      iat = 0
      do  ib = 1, nbas

        is = s_site(ib)%spec
        if (s_spec(is)%lmxa < 0) cycle
        iat = iat+1
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmax = s_spec(is)%rmt
        if (nr > nrmx) call rx('increase nrmx in tprodbas')
        call radmsh(s_spec(is)%rmt,a,nr,rofi)
        intopt = 10*nglob('lrquad')
        call radwgt(intopt,s_spec(is)%rmt,a,nr,rwgt)
        lcuta(iat) = lcutmx
        call mkrjk(E,nr,rofi,nrmx,lcuta(iat),lgrad,rkpr,rkmr)

C       Debugging: confirm normalization for l=0 orbitals
C       offr = npb(0+(lcutmx+1)*(iat-1))
C       nf   = npb(0+1+(lcutmx+1)*(iat-1))-offr
C       do  n1 = 1, nf
C         qtrue = dot3(nr,rprodx(1,offr+n1),rprodx(1,offr+n1),rwgt)
C         print *, n1,qtrue
C       enddo
C       call prrmsh('r*prodx',rofi,rprodx,nrmx,nr,npb(1))

C   ... Make sgbc
        if (mod(mode0/4,2) > 0 .and. lgrad > 0) then
          nsp = nglob('nsp')
C         call prrmsh('rhoc',rofi,s_site(ib)%rhoc,nr,nr,nsp)
          allocate(rhoc(nr))
          call dcopy(nr,s_site(ib)%rhoc,1,rhoc,1)
          if (nsp == 2) call daxpy(nr,1d0,s_site(ib)%rhoc(1,2),1,rhoc,1)
          qtrue = ddot(nr,rhoc,1,rwgt,1)
          fac = s_spec(is)%qc/qtrue
          call dscal(nr,fac,rhoc,1)
          call info5(30,0,0,' onsitecoul : iat = %i  z = %d  qcore = %;6,6d scaled to %;6,6d',iat,s_spec(is)%z,qtrue,qtrue*fac,5)
C         Render rhoc into r * (true core density)
          forall (ir = 2:nr) rhoc(ir) = rhoc(ir)/fpi/rofi(ir)
C         Add smoothed nuclear contribution
!         print *, '!!'; rhoc = 0
          do  ir = 1, nr
            call radgkl(rofi(ir),rgnuc,0,1,0,g0l)
            rhoc(ir) = rhoc(ir) - s_spec(is)%z * g0l(0) * rofi(ir)
          enddo
C         qtrue = fpi * dot3(nr,rofi,rhoc,rwgt); print *, qtrue; stop

C         Integrals 4*pi * int dr dr' (r B_R1μ(r)) V_l(r,r') (r' rhoc(r'))
          offr = npb(1+(lcutmx+1)*(iat-1))
          nf = npb(1+1+(lcutmx+1)*(iat-1))-offr
          if (nf == 0) goto 199
          l = 0; n1 = 1  ! number of core densities and l
C         vB(r) = 4*pi/(2l+1) * int dr' V_l(r,r') (r' rhoc(r'))
C         vB is proportional to r * [true potential * rhoc]
          call rintgp(11,rofi,a,rkpr(1,l,0),rhoc,-1d0,nrmx,nr,1,9,1,errmx,JBint)
          call rintgp(11,rofi,a,rkmr(1,l,0),rhoc,-1d0,nrmx,nr,1,9,1,errmx,HBint)
C         r*Hl(E,r) int_0^r dr' r'Jl(r') r'rhoc(r')  +  r*Jl(E,r) int_r^rmax dr' r'Hl(r') r'rhoc(r')
          forall (ir = 1:nr) vB(ir,n1,0) = rkmr(ir,l,0)*(JBint(1)-JBint(ir)) + rkpr(ir,l,0)*HBint(ir)
C         Radial derivative has 4 terms:
C         rkmr' * (int JB) + rkpr'*HBint + rkmr * (JBint)' + rkpr*(HBint)'
C         However (JBint)' = r J(r) r B(r)   and  (HBint)' = - r H(r) r B(r)
C         Thus 3rd and 4th terms cancel (not surprising: upper bound of one joins lower bound of the other)
C         Radial derivatives of rkpr and rkmr have been computed in call to mkrjk
          forall (ir = 1:nr) vB(ir,n1,1) = rkmr(ir,l,1)*(JBint(1)-JBint(ir)) + rkpr(ir,l,1)*HBint(ir)

C#ifdefC DEBUG3
CC         For this test to work, remove core density from rhoc
CC         Debugging: convert vb(:,1,0) to int v(r,r') rhoc(r')
CC         Compare against analytic answer [applies to nuclear charge only], which store in vb(:,1,1)
C          do  ir = 2, nr
C            vb(ir,1,0) = fpi/(2*l+1)*vb(ir,1,0)/rofi(ir)  ! convert to v(r)
C            call hansmr(rofi(ir),E,1/rgnuc,hs,1)  ! Sm-hankel / r^l
C            forall (l = 0:1) hs(l) = hs(l)*rofi(ir)**l
C            vb(ir,1,2) = -s_spec(is)%z * hs(0)  ! vb(:,1,2) = work array
C          enddo
C          print *, 'l=',l;call prrmsh('analytic, numerical intgrl vB (cols 2,4)',rofi,vb(:,1,0:2),nrmx,nr,3)
C
CC         Debugging: convert vb(:,1,1) to grad_r int v(r,r') r'B(r')
CC         Compare against results calculated two different ways:
CC         Overwrite column 0 with numerically differentiated vB
CC         Overwrite column 2 with analytic derivative of gaussian
C          vB(1,1,0) = (vB(2,1,0)*rofi(3)-vB(3,1,0)*rofi(2))/(rofi(3)-rofi(2))
C          call poldvm(rofi,vB,nr,8,.false.,1d-8,ir,vb(1,1,2))
C          vb(:,1,0) = vb(:,1,2)
C          do  ir = 2, nr
C            vb(ir,1,1) = fpi/(2*l+1)*vb(ir,1,1)/rofi(ir)  ! convert to grad int v(r,r') B(r')
C            call hansmr(rofi(ir),E,1/rgnuc,hs,1)  ! Sm-hankel / r^l
C            forall (l = 0:1) hs(l) = hs(l)*rofi(ir)**l
C            vb(ir,1,2) = s_spec(is)%z * hs(1) ! hs' = hs(l)*l/r - hs(l+1)
C          enddo
C          print *, 'l=',l;call prrmsh('grad intgrl vB (see src)',rofi,vb(:,1,0:2),nrmx,nr,3)
C#endif

C         Do all (nf) integrals with single dgemm call
          do  m = 0, 1
C           Fold in radial integration weights
            forall (ir = 1:nr) vB(ir,1,m) = vB(ir,1,m)*rwgt(ir)
            call dgemm('T','N',nf,1,nr,fpi,rprodx(1,offr+1),nrmx,vB(1,1,m),nrmx,0d0,sgbc(1,m,iat),mxnpbl)
          enddo

  199     continue
          deallocate(rhoc)
        endif

C   ... rojb
        if (mod(mode0,2) > 0) then
          fac = 1d0
          do  l = 0, lcuta(iat)
            fac = fac/(2*l+1)
            offr = npb(l+(lcutmx+1)*(iat-1))
            nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
            do  n1 = 1, nf
              do  m = 0, lgrad
                rojb(n1,l,m,iat) = fac*dot3(nr,rkpr(1,l,m),rprodx(1,offr+n1),rwgt)
              enddo
            enddo
          enddo
        endif
C       call yprmi('rojb for site ib=%i iat=%i',ib,iat,1+0,rojb(1,0,0,iat),0,mxnpbl,mxnpbl,lcuta(iat)+1)

C   ... sgbb = 4pi/(2l+1) int_0^r dr (r B_l(r)) * [...]  where
C       [...] = r*Hl(E,r) int_0^r dr' r'J(r') r'Bl(r') + r*Jl(E,r) int_r^rmax dr' r'H(r') r'Bl(r')
        if (mod(mode0/2,2) > 0) then
        do  l = 0, lcuta(iat)
          offr = npb(l+(lcutmx+1)*(iat-1))
          nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
          if (nf == 0) cycle
C         vB(r) = 4*pi/(2l+1) * int dr' V_l(r,r') (r' B_Rlν(r'))
C         vB is proportional to r * [true potential * B]
          do  n1 = 1, nf
            call rintgp(11,rofi,a,rkpr(1,l,0),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,JBint)
            call rintgp(11,rofi,a,rkmr(1,l,0),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,HBint)
C           r*Hl(E,r) int_0^r dr' r'Jl(r') r'Bl(r')  +  r*Jl(E,r) int_r^rmax dr' r'Hl(r') r'Bl(r')
            forall (ir = 1:nr) vB(ir,n1,0) = rkmr(ir,l,0)*(JBint(1)-JBint(ir)) + rkpr(ir,l,0)*HBint(ir)
          enddo

C          Comment this block and group integrals into a single dgemm call, below
C          do  n1 = 1, nf
C          do  n2 = 1, nf
C            sgbb(n1,n2,l,iat) = fpi/(2*l+1)*dot3(nr,vB(1,n1,0),rprodx(1,offr+n2),rwgt)
C          enddo
C          enddo

C         Fold in radial integration weights
          m = 0
          forall (ir = 1:nr, n1=1:nf) vB(ir,n1,m) = vB(ir,n1,m)*rwgt(ir)
C         Do all (nf,nf) integrals with single dgemm call
          call dgemm('T','N',nf,nf,nr,fpi/(2*l+1),vB(1,1,m),nrmx,rprodx(1,offr+1),nrmx,0d0,sgbb(1,1,l,iat),mxnpbl)

C         Symmetrize
          do  n1 = 1, nf
          do  n2 = n1, nf
            fac = sgbb(n1,n2,l,iat)/2 + sgbb(n2,n1,l,iat)/2
            errmx = max(errmx,abs(sgbb(n1,n2,l,iat)/2-sgbb(n2,n1,l,iat)/2))
            sgbb(n1,n2,l,iat) = fac; sgbb(n2,n1,l,iat) = fac
          enddo
          enddo
          call info2(60,0,0,' symmetrize sgbb for l=%i: errmx=%,3;3g',l,errmx)
C         call yprmi('sgbb for site ib=%i l=%i',ib,l,1+0,sgbb(1,1,l,iat),0,mxnpbl,nf,nf)
        enddo
        endif

C   --- Make ropjp, sgpb ---
        if (mode1 > 0) then

          allocate(phase(npwmb))
          do  ig = 1, npwmb
            phase(ig) = exp(img*s_lat%alat*ddot(3,qpg(1,ig),1,s_lat%pos(1,ib),1))
          enddo

C     ... Make scaling factors pjyl; see Eq 49 in ecalj manual (June 19, 2015)
C         <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
C         q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
          do  ig = 1, npwmb
!           phase(ig) = exp(img*s_lat%alat*ddot(3,qpg(1,ig),1,s_lat%pos(1,ib),1))
            ilm = 0
            do  l = 0, lcuta(iat)
              do  m = -l, l
                ilm = ilm+1
                pjyl(ilm,ig) = fpi*img**l*phase(ig)*yl(ig,ilm)
              enddo
            enddo
          enddo ! G vectors

C#ifdefC DEBUG
C          call info2(1,0,0,' sumcheck pjyl site %i %2;18,6D',ib,sum(pjyl))
C#ifdefC DEBUG2
C          call yprm0('(1p9e22.12)')
C          call yprmi('pjyl for site ib=%i q=%s,%3d',ib,q,3,pjyl,1,(lcutmx+1)**2,(lcutmx+1)**2,npwmb)
C#endifC
C#endif

C     ... Make rojp(G,ilm) = pjyl * W{J(E), J(G^2)}
          do  ig = 1, npwmb
            call wronjje(jobBes,qpg2(ig),E,rmax,1,1,lcutmx,wjjE)
            ilm = 0
            do  l = 0, lcuta(iat)
              do  m = -l, l
                ilm = ilm+1
                rojp(ig,ilm,iat) = (-wjjE(l))*pjyl(ilm,ig)
              enddo
            enddo
          enddo ! G vectors
C#ifdefC DEBUG
C          call info2(1,0,0,' sumcheck rojp site %i %2;18,6D',ib,sum(rojp(:,:,iat)))
C#endif

C     ... Setup for sgbp and sgpp
C         Note : the coulomb integral with J can be performed analytically.
C         But there is an error anyway because of sgpb is an integral of J with B
C         so we do it numerically.
C         On the other hand it would be possible for sigpp, i.e. <J | v | J>,
C         to be evaluated analytically.

C         Tabulate Bessel jG = j_l(sqrt(E)r) * r / (sqrt(E))**l
C         Note: slightly different convention than Jbar: jG = Jbar/(2l+1)!! [check]
          allocate(xi(nr,0:lcuta(iat)),y(nr),h(nr))
          do  ig = 1, npwmb
            call ropbes(rofi,qpg2(ig),lcuta(iat),y,h,xi,nr,1) ! MSM Standard def J/r^l
            do  ir = 1, nr
              r2s = rofi(ir)
              do  l = 0, lcuta(iat)
C               call besslr(qpg2(ig)*rofi(ir)**2,10,0,lcuta(iat),phi,psi)
C               jG(ir,ig,l,iat) = phi(l) * rofi(ir)**(l+1)
                jG(ir,ig,l) = xi(ir,l) * r2s    ! r*J(G^2)
                r2s = r2s*rofi(ir)
              enddo
            enddo
          enddo
C#ifdefC DEBUG
C          call info2(1,0,0,' sumcheck jg site %i %;18,6D',ib,sum(jG(nr,:,:)))
C#endif
          deallocate(xi,y,h)

C     ... sgpb = Int dr r*J_l  Int dr' v(r,r') r'*B_l(r')

C         This version is faster, but less accurate and versatile
C         ilm = 0
C         do  l = 0, lcuta(iat)
C           offr = npb(l+(lcutmx+1)*(iat-1))
C           nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
C           allocate(sig(npwmb,nf))
CC          Integral v(r,r') B(r') for each of nf functions
C           do  n1 = 1, nf
C             call rintgp(11,rofi,a,rkpr(1,l),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,JBint)
C             call rintgp(11,rofi,a,rkmr(1,l),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,HBint)
C             forall (ir = 1:nr) vB(ir,n1) = rkmr(ir,l)*(JBint(1)-JBint(ir)) + rkpr(ir,l)*HBint(ir)
C             do  ig = 1, npwmb
C               sig(ig,n1) = dot3(nr,vB(1,n1),jG(1,ig,l),rwgt)
C             enddo
C           enddo
C           do  m = -l, l
C             ilm = ilm+1
C             do  n1 = 1, nf
C             do  ig = 1, npwmb
C               sgpb(ig,n1,ilm,iat) = dconjg(pjyl(ilm,ig)) * sig(ig,n1)/(2*l+1)*fpi
C             enddo
C             enddo
C           enddo
C           deallocate(sig)
C         enddo
C         call info2(1,0,0,' sumcheck sgpb site %i %2;18,6D',ib,sum(sgpb(:,:,:,iat)))

C         Alternative integration of sgpb : Int dr r*B(r) Int dr' v(r,r')  r'*J(r')
C         More integrals than previous, maybe 2x slower than fast form above,
C         but perhaps more accurate since J is smooth.
C         Also integrals can be used for sgpp; they are accumulated for both sgpb and sgpp

          do  l = 0, lcuta(iat)
            do  ig = 1, npwmb
C             These integrals could be done analytically
              call rintgp(11,rofi,a,rkpr(1,l,0),jG(1,ig,l),-1d0,nrmx,nr,1,9,1,errmx,JBint)
              call rintgp(11,rofi,a,rkmr(1,l,0),jG(1,ig,l),-1d0,nrmx,nr,1,9,1,errmx,HBint)
C             vjG(r) = Int dr' v(r,r') r'*J(r')
              forall (ir = 1:nr) vjG(ir,ig,l) = rkmr(ir,l,0)*(JBint(1)-JBint(ir)) + rkpr(ir,l,0)*HBint(ir)
C             Fold in radial integration weights
              forall (ir = 1:nr) vjG(ir,ig,l) = vjG(ir,ig,l)*rwgt(ir)
            enddo
          enddo
          do  l = 0, lcuta(iat)
            offr = npb(l+(lcutmx+1)*(iat-1))
            nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
            allocate(sig(npwmb,nf))
C           Replace the following with faster dgemm call
C            do  n1 = 1, nf
C              do  ig = 1, npwmb
C                sig(ig,n1) = ddot(nr,vjG(1,ig,l),1,rprodx(1,offr+n1),1)/(2*l+1)*fpi
C              enddo
C            enddo
            call dgemm('T','N',npwmb,nf,nr,1d0,vjG(1,1,l),nrmx,rprodx(1,offr+1),nrmx,0d0,sig,npwmb)
            call dscal(npwmb*nf,fpi/(2*l+1),sig,1)

            do  ilm = l*l+1, (l+1)**2
              do  n1 = 1, nf
              do  ig = 1, npwmb
                sgpb(ig,n1,ilm,iat) = dconjg(pjyl(ilm,ig)) * sig(ig,n1)
              enddo
              enddo
            enddo
            deallocate(sig)
          enddo
C#ifdefC DEBUG
C          call info2(1,0,0,' sumcheck sgpb site %i %2;18,6D',ib,sum(sgpb(:,:,:,iat)))
C#endif

C     ... Make fouvb = 4*pi/(G^2-E) * integral(J*B) ... note E is negative here
          do  ig = 1, npwmb
            ilm = 0
            do  l = 0, lcuta(iat)
              offr = npb(l+(lcutmx+1)*(iat-1))
              nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
              allocate(sig(nf,1))
              do  n1 = 1, nf
                sig(n1,1) = dot3(nr,jG(1,ig,l),rprodx(1,offr+n1),rwgt)
              enddo
              do  m = -l, l
                ilm = ilm+1
                do  n1 = 1, nf
                  fouvb(ig,n1,ilm,iat) = fpi/(qpg2(ig)-E)*dconjg(pjyl(ilm,ig))*sig(n1,1)
                enddo
              enddo
              deallocate(sig)
            enddo
          enddo
C#ifdefC DEBUG
C          call info2(1,0,0,' sumcheck fouvb/10 site %i %2;18,6D',ib,sum(fouvb(:,:,:,iat))/10)
C#endif

C     --- Add sgpp + fouvp into vcoul (to save memory) ---
          allocate(sig(npwmb,npwmb),wjj(npwmb,npwmb,0:lcuta(iat)))
C     ... Wronskian wjj = W{qpg2(ig),qpg2(ig2)} at rmax for l=0:lcuta
          call wronjje(jobBes,qpg2,qpg2,rmax,npwmb,npwmb,lcuta(iat),wjj)
          do  l = 0, lcuta(iat)
            call dgemm('T','N',npwmb,npwmb,nr,1d0,vjG(1,1,l),nrmx,jG(1,1,l),nrmx,0d0,sig,npwmb)
            call dscal(npwmb*npwmb,fpi/(2*l+1),sig,1)

            do  ig = 1, npwmb

C             if (ib == 2 .and. ig == 2) call snot

              do  ig2 = 1, ig
                p12 = 0
                do  ilm = l*l+1, (l+1)**2
C                 p12 = p12 + dconjg(pjyl(ilm,ig)*phase(ig)) * pjyl(ilm,ig2)*phase(ig2)
                  p12 = p12 + dconjg(pjyl(ilm,ig)) * pjyl(ilm,ig2)
                enddo

                wv12 = fpi/(qpg2(ig)-E)*wjj(ig,ig2,l)
                wv21 = fpi/(qpg2(ig2)-E)*wjj(ig2,ig,l) ! Could exploit wjj(ig,ig2,l) = wjj(ig2,ig,l)

C               Add onsite contribution (prop sig) and Fourier contribution (prop wv12 and wv21) to vcoul
                vcoul(ntpb+ig,ntpb+ig2) = vcoul(ntpb+ig,ntpb+ig2) + (sig(ig,ig2) + wv12 + wv21) * p12
              enddo ! ig2
            enddo ! ig
          enddo ! l
          deallocate(sig,wjj)

          deallocate(phase)

        endif  ! mode1 > 0
      enddo  ! Loop over sites
      deallocate(rkpr,rkmr)
      deallocate(JBint,Hbint,vB)

      if (mode1 > 0) then

C#ifdefC DEBUG
C        call info2(1,0,0,' sumcheck onsite vcoul(PW-PW)/1000 %2;18,6D',sum(vcoul(ntpb:ntmb,ntpb:ntmb))/1000,2)
C#ifdefC DEBUG2
C        call zprm0('(1p9e22.12)')
C        call zprm('vcoul (pp)',2,vcoul(ntpb+1,ntpb+1),ntmb,npwmb,npwmb)
C#endifC
C#endif

        deallocate(qpg,yl,qpg2,pjyl)
        deallocate(jG,vjG)

      endif

      call tcx('onsite-vc')

      end

      subroutine mkrjk(E,nr,rofi,nrmx,lcut,lgrad,rkpr,rkmr)
C- Make scaled spherical Bessel and Hankel functions on a radial mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   E     :Energy of Bessel, Hankel (must be <= 0)
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   nrmx  :leading dimension of rofi,rkpr,rkmr
Ci   lcut  :Return rkpr, rkmr for l=0:lcut
Ci   lgrad :if 1, return grdients of rkpr, rkmr for l=0:lcut (see Remarks)
Co Outputs
Ci   rkpr  :r*Jbarl(E,r) (Spherical Bessel functions; see Remarks)
Ci   rkmr  :r*Hbarl(E,r) (Spherical Hankel functions; see Remarks)
Ci   ... if lgrad=1, r*Jbarl'(E,r) and r*Hbarl'(E,r)  are returned
Cr Remarks
Cr   rkpr, rkmr are real functions.
Cr   They are related to the true j_l and h_l as:
Cr     rkpr = (2l+1)!! * j_l(i sqrt(abs(E)) r) * r / (i sqrt(abs(E)))**l
Cr     rkmr = (2l-1)!! * h_l(i sqrt(abs(E)) r) * r * i*(i sqrt(abs(E)))**(l+1)
Cr   They have this property:
Cr     rkpr -> r**l*r      as E->0 (OKA definition * (2(2l+1)); see besslr)
Cr     rkmr -> r**(-l-1)*r as E->0 (OKA definition; see besslr)
Cr   See OKA=2 description in radkj.f
Cu Updates
Cu   20 Dec 17 Adapted from the old GW genhj
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer:: nr,nrmx,lcut,lgrad
      real(8):: E,rofi(nrmx),rkpr(nrmx,0:lcut,0:lgrad),rkmr(nrmx,0:lcut,0:lgrad)
C ... Local parameters
      integer ir,l
      double precision ak(0:lcut),aj(0:lcut),dk(0:lcut),dj(0:lcut)

      rkpr(1,:,:) = 0; rkmr(1,:,:) = 0
      do  ir  = 2, nr
        call radkj(E,rofi(ir),lcut,ak,aj,dk,dj,20)
        do  l = 0, lcut
          rkpr(ir,l,0) = rofi(ir)*aj(l)
          rkmr(ir,l,0) = rofi(ir)*ak(l)
          if (lgrad > 0) then
            rkpr(ir,l,1) = rofi(ir)*dj(l)
            rkmr(ir,l,1) = rofi(ir)*dk(l)
          endif
        enddo
      enddo

C      OLD
C      real(8):: dfac,rl,psi(0:lcut),phi(0:lcut)
C      do  ir  = 2, nr
C        call besslr(E*rofi(ir)**2,0,0,lcut,phi,psi)
C        call radkj(E,rofi(ir),lcut,ak,aj,dk,dj,20)
C        rl = 1
C        dfac = 1
C        do  l = 0, lcut
C          rkpr(ir,l,0) = phi(l)*(rl*rofi(ir))*(dfac*(2*l+1))
C          rkmr(ir,l,0) = psi(l)/rl/dfac
C          dfac = dfac*(2*l+1)
C          rl = rl*rofi(ir)
C          print *, l, rofi(ir)*aj(l)/rkpr(ir,l,0)-1, rofi(ir)*ak(l)/rkmr(ir,l,0)-1
C          rkpr(ir,l,0) = rofi(ir)*aj(l)
C          rkmr(ir,l,0) = rofi(ir)*ak(l)
C          if (lgrad > 0) then
C            rkpr(ir,l,1) = rofi(ir)*dj(l)
C            rkmr(ir,l,1) = rofi(ir)*aj(l)
C          endif
C        enddo
C      enddo
      end

C      subroutine dbgonsitecoul(s_lat,s_site,s_spec,mode,lgrad,nbas0,lcutmx0,mxnpbl0,npb,E,nrmx,rprodx,
C     .  q,ntpb,npwmb,ntmb,rojb,sgbb,sgbc,rojp,sgpb,fouvb,vcoul)
CC- Make onsite Coulomb integrals in all augmentation spheres
CC ----------------------------------------------------------------------
CCio Structures
CCio  s_site :struct for site-specific data; see structures.h
CCi     Elts read:  spec
CCo     Stored:     *
CCo     Allocated:  *
CCio    Elts passed:*
CCio    Passed to:  *
CCio  s_spec :struct for species-specific data; see structures.h
CCi     Elts read:  lmxa a nr
CCo     Stored:     *
CCo     Allocated:  *
CCio    Elts passed:rmt
CCio    Passed to:  *
CCi Inputs
CCi   mode  :controls what is made
CCi         :1s digit
CCi         :0 make none of rojb, sgbb, sgbc
CCi         :1 make rojb
CCi         :2 make sgbb
CCi         :4 make sgbc (requires lgrad also be set)
CCi         :Any combination is allowed
CCi         :10s digit
CCi         :1 make rojp,sgpb,fouvb, on-site of PW part of vcoul
CCi   lgrad :0 make coulomb integrals
CCi         :1 Also make integrals needed for gradients of coulomb interaction (not finished)
CCi   nbas  :size of basis
CCi   lcutmx:maximum l-cutoff for product basis functions, typically 2*(nl-1)
CCi         :Dimensions npb,rprodx,rojb,sgbb
CCi         :Called lxx in old GW code.
CCi   lcuta :site-dependent lcut for integrals rojb,sgbb,...
CCi         :Called lx in old GW code.  For now, lcuta = lcutmx
CCi   mxnpbl:maximum number of radial B functions; dimensions rojb,sgbb
CCi         :Called nxx in old GW code.
CCi   npb   :Table of offsets to radial functions in rprodx
CCi         :Functions npb(l,iat):npb(l+1,iat) are the family of radial B functions
CCi         :stored in rprodx for quantum number l and site iat.
CCi         :Formerly called nx in old GW code
CCi   E     :Use Yukawa potential for coulomb interaction, 1/sqrt(-E) is screening length
CCi   nrmx  :Maxiumum number of radial mesh points, and leading dimension of rprodx
CCi         :Formerly called nrx in old GW code
CCi   rprodx:Product basis functions * r: r*BRlμ(r)
CCi         :Normalized such that \int_0^rmax (r*BRlμ(r))^2 dr = 1
CCi   q     :Make integrals for wave number q
CCi   ntpb  :Total total number of local product basis functions within MTs
CCi         :Formerly called nbloch in old GW code
CCi   npwmb :number of G vectors for Coulomb interaction and
CCi         :number of interstitial product functions for this q
CCi         :Called ngc in old GW code.
CCi   ntmb  :largest dimension of mixed product basis functions = ntpb + npmbx
CCi         :Used here solely as leading dimension of vcoul
CCi         :Formerly called nblochpmx in old GW code
CCo Outputs
CCo   rojb  :Onsite integral rojb(BRlμ) = 1/(2l+1)!! int dr r*Jl(E,r) r*BRlμ(r)
CCo         :Needed to evaluate integrals <B^k_RLμ(r)|V_l(r,r')|B^k_R'L'μ'(r')>
CCo         :when (R,T) and (R',T') are different since v expanded in J_L(E,r)
CCo         :Here V_l(r,r') = (r<)J_l(E,r<) (r>)H_l(E,r>)
CCo         :Not made, but can do so (see below), if lgrad=1 (not sure if needed)
CCo         :rojb(:,:,1,:) = 1/(2l+1)!! int dr r*J'_l(E,r) r*BRlμ(r)
CCo   sgbb  :Onsite integrals <B|v(onsite)|B> corresponding to rojb when (R,T) = (R',T')
CCo         :sgbb = sigma integral, given by Eq 43 in ecalj manual (June 19, 2015)
CCo         :Follows from the 1-center expansion of the Yukawa potential
CCo         :sgbb(B_Rlμ, B_Rlν)
CCo         :  = 4*pi/(2l+1) * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r B_Rlμ(r)) (r' B_Rlν(r'))
CCo         :  = 4*pi/(2l+1) int dr dr' (r B_Rlμ(r)) V_l(r,r') (r' B_Rlν(r'))
CCo         :where V_l(r,r') = (r<)J_l(E,r<) (r>)H_l(E,r>)
CCo   rojp  :Onsite integral <j_L(E,r) | P(q+G)_L > where
CCo         : | P(q+G)_L> = projection of exp(i (q+G) r) to channnel L
CCo         :  | j_L(E=0,r)> = r^l/(2l+1)!! Y_L [check]
CCo         : rojp(G,L) = pjyl(L,G) * W{J_L(E), J_L(G^2)}
CCo         : pjyl is defined in Eq 49 in ecalj manual (June 19, 2015)
CCo         : Needed to evaluate <P_G|v|B> since v expanded in J_L(E,r)
CCo   sgpb  :Onsite integrals <J|v(onsite)|B> corresponding to rojp when (R,T) = (R',T')
CCo         :Here J=J(G^2) --- 1-center expansion of plane waves
CCo         :  sgpb(J(G^2), B_Rlν) =
CCo         :  4*pi/(2l+1) (pjyl)+ * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r J_L(r,G^2)) (r' B_Rlν(r'))
CCo   fouvb :Fourier integral of PW and product basis inside augmentation sphere
CCo         :  4*pi/(G^2-E) * integral(J*B) ... note E is negative
CCo Outputs
CCo   vcoul :On-site + Fourier contribution to PW-PW part of vcoul is computed (see Remarks)
CCl Local variables
CCl   rkpr  :r*Jl(E,r),  Jbar definition (see Remarks)
CCl         :if lgrad is 1, r*J'l(E,r) is also computed
CCl   rkmr  :r*Hl(E,r),  Hbar definition (see Remarks)
CCl         :if lgrad is 1, r*H'l(E,r) is also computed
CCl   jG    :r*Jl(G^2,r), MSM Standard definition (besslr)  [called ajr in old mkjp_4]
CCl   vB    :Coulomb integral : int v(r,r') (r' * B(r'))
CCl   vjG   :Coulomb integral : int v(r,r') jG(r')
CCr Remarks
CCr
CCr *Expansion theorem for bare coulomb interaction:
CCr   exp(-lamda |r-r'|)/|r-r'| =
CCr   sum_lm  4*pi/(2l+1) (r<) J_l(E,r<) (r>) H_l(E,r>) (-1)^m Y_(l,m)(r) Y_(l,-m)(r')
CCr
CCr *Here the Jl and Hl are "bar Bessel and Hankel functions"
CCr  with the property Jl -> r^l  and Hl -> r^-l-1 for E=0,
CCr  are proportional to the radial part of customary spherical Hankels and Bessels.
CCr  See mkrjk for definition.
CCr
CCr *Expansion theorem is readily derived from the generating function for Legendre polynomials
CCr   1/sqrt(1+h^2-2*h*cos(gamma)) = sum_l h^l P_l(cos(gamma))
CCr *and the addition theorem, using the angle between r and r' for gamma
CCr   P_l(cos(gamma)) = 4*pi/(2l+1) sum_m (-1)^m Y_(l,m)(r) Y_(l,-m)(r')
CCr
CCr  See Section 18.4 of Takao's ecalj manual (June 2015)
CCr
CCr *Remarks on sgpp,fouvp, and vcoul.
CCr  For consistency, sgpp and fouvp should be returned as arrays together with
CCr  sgpb and fouvb, and vcoulq should not be generated here.
CCr
CCr  However, since sgpp and fouvp use up lots of memory, which moreover scales as N^3,
CCr  they are made locally and their contribution directly added to vcoulq.
CCr
CCr  If it is desired that the on-site contribution vcoulq be made elsewhere,
CCr  one alternative is to return jG and vjG from this routine and make it on the fly.
CCr  jG and vjG would have to be decomposed by site.
CCu Updates
CCu   20 Dec 17 Adapted from Kotani's mkjb_4, mkjp_4 in old GW code.
CC ----------------------------------------------------------------------
C      use structures
C      implicit none
CC ... For structures
C!     include 'structures.h'
C      type(str_lat)::   s_lat
C      type(str_site)::  s_site(*)
C      type(str_spec)::  s_spec(*)
C
CC ... Turn passed args into parameters for debugging
C      integer, parameter :: nbas=1, mxnpbl=1, lcutmx=2
C      double precision rsm,g0l(0:lcutmx),hs(0:lcutmx),dhs(0:lcutmx),ddhs(0:lcutmx)
C      double precision xx(10)
C
CC ... Passed parameters
C      integer :: mode,lgrad,nbas0,nrmx,lcutmx0,ntpb,npwmb,ntmb,mxnpbl0,npb(0:lcutmx*(nbas+1))
C!     integer lcuta(nbas)
C      real(8) :: E,rprodx(nrmx,*)
C      real(8) :: rojb(mxnpbl,0:lcutmx,0:lgrad,nbas) ! rho-type onsite integral
C      real(8) :: sgbb(mxnpbl,mxnpbl,0:lcutmx,nbas) ! sigma-type onsite integral
C      real(8) :: sgbc(mxnpbl,0:1,nbas) ! sigma-type onsite integral for nuc+core density and B
CC     Specific to PW-related integrals
C      real(8) :: q(3)
C      complex(8) :: rojp(npwmb,(lcutmx+1)**2,nbas) ! rho-type onsite integral
C      complex(8) :: sgpb(npwmb,mxnpbl,(lcutmx+1)**2,nbas)
C      complex(8) :: fouvb(npwmb,mxnpbl,(lcutmx+1)**2,nbas)
C      complex(8) :: vcoul(ntmb,ntmb)
C!     complex(8) :: sgpp(npwmb,npwmb,nat)
C
CC ... Dynamically allocated local arrays
C      real(8), allocatable :: rkpr(:,:,:),rkmr(:,:,:)
C      real(8), allocatable :: JBint(:),HBint(:),vB(:,:,:)
C      real(8), allocatable :: qpg(:,:),yl(:,:),qpg2(:)
C      complex(8),allocatable :: pjyl(:,:),phase(:)
CC     Help arrays to make fast Bessel functions
C      real(8), allocatable :: xi(:,:),y(:),h(:),jG(:,:,:),vjG(:,:,:),sig(:,:),wjj(:,:,:),rhoc(:)
C
CC ... Local parameters
C      integer, parameter :: jobBes=0
C      integer lcuta(nbas)  ! for now just copy from lcutmx
C      integer :: l,n1,n2,ib,iat,is,offr,nf,intopt,ir,nr,mode0,mode1,ig,ig2,ilm,m,nsp
C      double precision fac,fpi,pi,a,errmx,tpiba,rmax,r2s,wv12,wv21,qtrue,rgnuc
C      double precision rofi(nrmx),rwgt(nrmx),wjjE(0:lcutmx)
C      complex(8), parameter :: img=(0d0,1d0)
C      complex(8) :: p12
C      procedure(integer) :: nglob
C      procedure(real(8)) :: ddot,dot3,dlength
C
C      call tcn('onsite-vc')
C
CC --- Special initialization for debugging ---
CC     One function per l
C      forall (ib = 1:lcutmx+1) npb(ib) = ib
C      npb(lcutmx+2) = npb(lcutmx+1)
CC     Radial integration points and weights
C      ib = 1
C      is = s_site(ib)%spec
C      a = s_spec(is)%a
C      nr = s_spec(is)%nr
C      rmax = s_spec(is)%rmt
C      call radmsh(s_spec(is)%rmt,a,nr,rofi)
C      intopt = 10*nglob('lrquad')
C      call radwgt(intopt,s_spec(is)%rmt,a,nr,rwgt)
CC     Use g_0l for radial function.  Note rprodx = r * true fn
C      rsm = 2*rmax/3 ;  print *, '!!';     rsm = 1*rmax/4
C      pi  = 4d0*datan(1d0)
C      fpi = 4*pi
C      do  ir = 1, nr
C        call radgkl(rofi(ir),rsm,0,lcutmx,0,g0l)
C        forall (l = 0:lcutmx) rprodx(ir,l+1) = g0l(l) * rofi(ir)**(l+1)
C      enddo
CC     Debugging : confirm that nabla hs = -4*pi*g0l
C      do  ir = nr, nr
C        call hanszd(2,rofi(ir),0d0,-rsm,lcutmx,hs,dhs,ddhs,xx,xx,xx)
C        call info5(1,0,0,' r = %;8,8F   ddhs     =%n:2;12,12F',rofi(ir),1+lcutmx,ddhs,4,5)
C        call info5(1,0,0,' r = %;8,8F  -4*pi*g0l =%n:2;12,12F',rofi(ir),1+lcutmx,-fpi/rofi(ir)*rprodx(nr,1:1+lcutmx),4,5)
CC       Check hansmr
C        call hansmr(rofi(ir),E*0,1/rsm,ddhs,1)
C        forall (l = 0:1) ddhs(l) = ddhs(l)*rofi(ir)**l
C      enddo
CC    call prrmsh('r*g',rofi,rprodx,nrmx,nr,1+lcutmx)
C
CC --- Initialize ---
C      pi  = 4d0*datan(1d0)
C      fpi = 4*pi
C      tpiba = 2*pi/s_lat%alat
C      mode0 = mod(mode,10)
C      mode1 = mod(mode/10,10)
C      call sanrg(.true.,lgrad,0,1,'onsitecoul','lgrad')
C
C      if (mod(mode0,2) > 0) call dpzero(rojb,size(rojb))
C      if (mod(mode0/2,2) > 0) call dpzero(sgbb,size(sgbb))
C
CC     allocate(rkpr(nrmx,0:lcutmx,0:lgrad),rkmr(nrmx,0:lcutmx,0:lgrad)) ! Barred H,J and radial derivatives
C      allocate(rkpr(nrmx,0:lcutmx+1,0:lgrad),rkmr(nrmx,0:lcutmx+1,0:lgrad))  ! Extra l probably needed
C      allocate(JBint(nrmx),Hbint(nrmx),vB(nrmx,mxnpbl,0:2))  ! Space for work array
C      if (mode1 > 0) then
C        allocate(qpg(3,npwmb),yl(npwmb,(lcutmx+1)**2),qpg2(npwmb),pjyl((lcutmx+1)**2,npwmb))
C        forall (ig = 1:npwmb) qpg(1:3,ig) = tpiba * (q(1:3) + matmul(s_lat%qlat,s_lat%igvcc(1:3,ig)))
C        call ropyln(npwmb,qpg(1,:),qpg(2,:),qpg(3,:),lcutmx,npwmb,yl,qpg2)
C        allocate(jG(nrmx,npwmb,0:lcutmx),vjG(nrmx,npwmb,0:lcutmx))
C      endif
C
CC --- For each augmentation sphere, do ---
C      iat = 0
C      do  ib = 1, nbas
C
C        is = s_site(ib)%spec
C        if (s_spec(is)%lmxa < 0) cycle
C        iat = iat+1
C        a = s_spec(is)%a
C        nr = s_spec(is)%nr
C        rmax = s_spec(is)%rmt
C        if (nr > nrmx) call rx('increase nrmx in tprodbas')
C        call radmsh(s_spec(is)%rmt,a,nr,rofi)
C        intopt = 10*nglob('lrquad')
C        call radwgt(intopt,s_spec(is)%rmt,a,nr,rwgt)
C        lcuta(iat) = lcutmx
C        call mkrjk(E,nr,rofi,nrmx,lcuta(iat)+1,lgrad,rkpr,rkmr)
C
CC       Debugging: confirm normalization for l=0 orbitals
CC       offr = npb(0+(lcutmx+1)*(iat-1))
CC       nf   = npb(0+1+(lcutmx+1)*(iat-1))-offr
CC       do  n1 = 1, nf
CC         qtrue = dot3(nr,rprodx(1,offr+n1),rprodx(1,offr+n1),rwgt)
CC         print *, n1,qtrue
CC       enddo
CC       call prrmsh('r*prodx',rofi,rprodx,nrmx,nr,npb(1))
C
CC   ... Make sgbc
C        if (mod(mode0/4,2) > 0 .and. lgrad > 0) then
C          nsp = nglob('nsp')
CC         call prrmsh('rhoc',rofi,s_site(ib)%rhoc,nr,nr,nsp)
C          allocate(rhoc(nr))
C          call dcopy(nr,s_site(ib)%rhoc,1,rhoc,1)
C          if (nsp == 2) call daxpy(nr,1d0,s_site(ib)%rhoc(1,2),1,rhoc,1)
C          qtrue = ddot(nr,rhoc,1,rwgt,1)
C          fac = s_spec(is)%qc/qtrue
C          call dscal(nr,fac,rhoc,1)
CC         Render rhoc into r * (true core density)
C          forall (ir = 2:nr) rhoc(ir) = rhoc(ir)/fpi/rofi(ir)
CC         Add smoothed nuclear contribution
C          print *, '!!'; rhoc = 0
C          rgnuc = .02d0
C          do  ir = 1, nr
C            call radgkl(rofi(ir),rgnuc,0,1,0,g0l)
C            rhoc(ir) = rhoc(ir) - s_spec(is)%z * g0l(0) * rofi(ir)
C          enddo
C          qtrue = fpi * dot3(nr,rofi,rhoc,rwgt)
C          print *, qtrue
C
CC         Integrals 4*pi * int dr dr' (r B_R1μ(r)) V_l(r,r') (r' rhoc(r'))
C          offr = npb(1+(lcutmx+1)*(iat-1))
C          nf = npb(1+1+(lcutmx+1)*(iat-1))-offr
C          if (nf == 0) goto 199
C          l = 0; n1 = 1
CC         vB(r) = 4*pi/(2l+1) * int dr' V_l(r,r') (r' rhoc(r'))
CC         vB is proportional to r * [true potential * rhoc]
C          call rintgp(11,rofi,a,rkpr(1,l,0),rhoc,-1d0,nrmx,nr,1,9,1,errmx,JBint)
C          call rintgp(11,rofi,a,rkmr(1,l,0),rhoc,-1d0,nrmx,nr,1,9,1,errmx,HBint)
CC         r*Hl(E,r) int_0^r dr' r'Jl(r') r'rhoc(r')  +  r*Jl(E,r) int_r^rmax dr' r'Hl(r') r'rhoc(r')
C          forall (ir = 1:nr) vB(ir,n1,0) = rkmr(ir,l,0)*(JBint(1)-JBint(ir)) + rkpr(ir,l,0)*HBint(ir)
CC         Radial derivative has 4 terms:
CC         rkmr' * (int JB) + rkpr'*HBint + rkmr * (JBint)' + rkpr*(HBint)'
CC         Note (JBint)' = r J(r) r B(r)   and  (HBint)' = - r H(r) r B(r)
CC         Thus 3rd and 4th terms cancel (not surprising: upper bound of one joins lower bound of the other)
CC         Radial derivative of of rkpr and rkmr are pre-computed in mkrjk.
C          forall (ir = 1:nr) vB(ir,n1,1) = rkmr(ir,l,1)*(JBint(1)-JBint(ir)) + rkpr(ir,l,1)*HBint(ir)
C
CC         Debugging: convert vb(:,1,0) to int v(r,r') rhoc(r')
CC         Compare against analytic answer [applies to nuclear charge only], which store in vb(:,1,1)
C          do  ir = 2, nr
C            vb(ir,1,0) = fpi/(2*l+1)*vb(ir,1,0)/rofi(ir)  ! convert to v(r)
C            call hansmr(rofi(ir),E,1/rgnuc,hs,1)  ! Sm-hankel / r^l
C            forall (l = 0:1) hs(l) = hs(l)*rofi(ir)**l
C            vb(ir,1,2) = -s_spec(is)%z * hs(0)  ! vb(:,1,2) = work array
C          enddo
C          print *, 'l=',l;call prrmsh('analytic, numerical intgrl vB (cols 2,4)',rofi,vb,nrmx,nr,3)
C
CC         Debug: convert vb(:,1,1) to grad_r int v(r,r') r'B(r')
CC         Compare against results calculated two different ways:
CC         Overwrite column 0 with numerically differentiated vB
CC         Overwrite column 2 with analytic derivative of gaussian
C          vB(1,1,0) = (vB(2,1,0)*rofi(3)-vB(3,1,0)*rofi(2))/(rofi(3)-rofi(2))
C          call poldvm(rofi,vB,nr,8,.false.,1d-8,ir,vb(1,1,2))
C          vb(:,1,0) = vb(:,1,2)
C          do  ir = 2, nr
C            vb(ir,1,1) = fpi/(2*l+1)*vb(ir,1,1)/rofi(ir)  ! convert to grad int v(r,r') B(r')
C            call hansmr(rofi(ir),E,1/rgnuc,hs,1)  ! Sm-hankel / r^l
C            forall (l = 0:1) hs(l) = hs(l)*rofi(ir)**l
C            vb(ir,1,2) = s_spec(is)%z * hs(1) ! hs' = hs(l)*l/r - hs(l+1)
C          enddo
C          print *, 'l=',l;call prrmsh('grad intgrl vB (see src)',rofi,vb,nrmx,nr,3)
C
CC         Fold in radial integration weights
C          forall (ir = 1:nr, m=0:1, n1=1:nf) vB(ir,n1,m) = vB(ir,n1,m)*rwgt(ir)
CC         Do all (nf) integrals with single dgemm call
C          do  m = 0, 1
C          call dgemm('T','N',1,nf,nr,fpi,vB(1,1,m),nrmx,rprodx(1,offr+1),nrmx,0d0,sgbc(1,m,iat),mxnpbl)
C          enddo
C          call rx0('done')
C
C  199     continue
C        endif
C
CC   ... rojb
C        if (mod(mode0,2) > 0) then
C          fac = 1d0
C          do  l = 0, lcuta(iat)
C            fac = fac/(2*l+1)
C            offr = npb(l+(lcutmx+1)*(iat-1))
C            nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
C            do  n1 = 1, nf
C              do  m = 0, lgrad
C                rojb(n1,l,m,iat) = fac*dot3(nr,rkpr(1,l,m),rprodx(1,offr+n1),rwgt)
C              enddo
C            enddo
C          enddo
C        endif
CC       call yprmi('rojb for site ib=%i iat=%i',ib,iat,1+0,rojb(1,0,0,iat),0,mxnpbl,mxnpbl,lcuta(iat)+1)
C
CC   ... sgbb = 4pi/(2l+1) int_0^r dr (r B_l(r)) * [...]  where
CC       [...] = r*Hl(E,r) int_0^r dr' r'J(r') r'Bl(r') + r*Jl(E,r) int_r^rmax dr' r'H(r') r'Bl(r')
C        if (mod(mode0/2,2) > 0) then
C        do  l = 0, lcuta(iat)
C          offr = npb(l+(lcutmx+1)*(iat-1))
C          nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
C          if (nf == 0) cycle
CC         vB(r) = 4*pi/(2l+1) * int dr' V_l(r,r') (r' B_Rlν(r'))
CC         vB is proportional to r * [true potential * B]
C          do  n1 = 1, nf
C            call rintgp(11,rofi,a,rkpr(1,l,0),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,JBint)
C            call rintgp(11,rofi,a,rkmr(1,l,0),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,HBint)
CC           r*Hl(E,r) int_0^r dr' r'Jl(r') r'Bl(r')  +  r*Jl(E,r) int_r^rmax dr' r'Hl(r') r'Bl(r')
C            forall (ir = 1:nr) vB(ir,n1,0) = rkmr(ir,l,0)*(JBint(1)-JBint(ir)) + rkpr(ir,l,0)*HBint(ir)
CC           Radial derivative has 4 terms:
CC           rkmr' * (int JB) + rkpr'*HBint + rkmr * (JBint)' + rkpr*(HBint)'
CC           Note (JBint)' = r J(r) r B(r)   and  (HBint)' = - r H(r) r B(r)
CC           Thus 3rd and 4th terms cancel (not surprising: upper bound of one joins lower bound of the other)
CC           For radial derivative of h, see JMP 39, 3393, Eq. 4.7   h' = l/r h_l - h_l+1
CC           For radial derivative of j, see radkj    j' = l/r*J_l - e*J_l+1
CC           Derivatives require r * d/dr (rkmr/r) and r * d/dr (rkpr/r) since rkmr,rkpr are h,j scaled by r
C            if (lgrad == 1) then
C              vB(1,n1,1) = 0
C              do  ir = 2, nr
C                vB(ir,n1,1) = rkmr(ir,l,1)*(JBint(1)-JBint(ir)) + rkpr(ir,l,1)*HBint(ir)
C              enddo
C            endif
C          enddo
CC         Debug: convert vb(:,1,0) to int v(r,r') r'B(r')
CC         Compare against analytic answer, which store in vb(:,1,1)
CC         Note constant shift if gaussian is not entirely contained in sphere.
C          do  ir = 2, nr
C            vb(ir,1,0) = fpi/(2*l+1)*vb(ir,1,0)/rofi(ir)  ! convert to v(r)
C            call hanszd(2,rofi(ir),0d0,-rsm,lcutmx,hs,dhs,ddhs,xx,xx,xx)
C            vb(ir,1,2) = hs(l) ! vb(:,1,2) = work array
C          enddo
C          print *, 'l=',l;call prrmsh('analytic, numerical intgrl vB (cols 2,4)',rofi,vb,nrmx,nr,3)
C
CC         Debug: convert vb(:,1,1) to grad_r int v(r,r') r'B(r')
CC         Compare against results calculated two different ways:
CC         Overwrite column 0 with numerically differentiated vB
CC         Overwrite column 2 with analytic derivative of gaussian
C          vB(1,1,0) = (vB(2,1,0)*rofi(3)-vB(3,1,0)*rofi(2))/(rofi(3)-rofi(2))
C          call poldvm(rofi,vB,nr,8,.false.,1d-8,ir,vb(1,1,2))
C          vb(:,1,0) = vb(:,1,2)
C          do  ir = 2, nr
C            vb(ir,1,1) = fpi/(2*l+1)*vb(ir,1,1)/rofi(ir)  ! convert to grad int v(r,r') B(r')
C            call hanszd(2,rofi(ir),0d0,-rsm,lcutmx,hs,dhs,ddhs,xx,xx,xx)
C            vb(ir,1,2) = dhs(l) ! Analytic derivative for potential from Gaussian
C          enddo
C          print *, 'l=',l;call prrmsh('grad intgrl vB (see src)',rofi,vb,nrmx,nr,3)
C
C          cycle
C
CC         Fold in radial integration weights
C          m = 0
C          forall (ir = 1:nr, n1=1:nf) vB(ir,n1,m) = vB(ir,n1,m)*rwgt(ir)
CC         Do all (nf,nf) integrals with single dgemm call
C          call dgemm('T','N',nf,nf,nr,fpi/(2*l+1),vB(1,1,m),nrmx,rprodx(1,offr+1),nrmx,0d0,sgbb(1,1,l,iat),mxnpbl)
C
CC         Symmetrize
C          do  n1 = 1, nf
C          do  n2 = n1, nf
C            fac = sgbb(n1,n2,l,iat)/2 + sgbb(n2,n1,l,iat)/2
C            errmx = max(errmx,abs(sgbb(n1,n2,l,iat)/2-sgbb(n2,n1,l,iat)/2))
C            sgbb(n1,n2,l,iat) = fac; sgbb(n2,n1,l,iat) = fac
C          enddo
C          enddo
C          call info2(60,0,0,' symmetrize sgbb for l=%i: errmx=%,3;3g',l,errmx)
CC         call yprmi('sgbb for site ib=%i l=%i',ib,l,1+0,sgbb(1,1,l,iat),0,mxnpbl,nf,nf)
C        enddo
C
C        call rx0('done')
C
C        endif
C
CC   --- Make ropjp, sgpb ---
C        if (mode1 > 0) then
C
C          allocate(phase(npwmb))
C          do  ig = 1, npwmb
C            phase(ig) = exp(img*s_lat%alat*ddot(3,qpg(1,ig),1,s_lat%pos(1,ib),1))
C          enddo
C
CC     ... Make scaling factors pjyl; see Eq 49 in ecalj manual (June 19, 2015)
CC         <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
CC         q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
C          do  ig = 1, npwmb
C!           phase(ig) = exp(img*s_lat%alat*ddot(3,qpg(1,ig),1,s_lat%pos(1,ib),1))
C            ilm = 0
C            do  l = 0, lcuta(iat)
C              do  m = -l, l
C                ilm = ilm+1
C                pjyl(ilm,ig) = fpi*img**l*phase(ig)*yl(ig,ilm)
C              enddo
C            enddo
C          enddo ! G vectors
C
CC#ifdefC DEBUG
CC          call info2(1,0,0,' sumcheck pjyl site %i %2;18,6D',ib,sum(pjyl))
CC#ifdefC DEBUG2
CC          call yprm0('(1p9e22.12)')
CC          call yprmi('pjyl for site ib=%i q=%s,%3d',ib,q,3,pjyl,1,(lcutmx+1)**2,(lcutmx+1)**2,npwmb)
CC#endifC
CC#endif
C
CC     ... Make rojp(G,ilm) = pjyl * W{J(E), J(G^2)}
C          do  ig = 1, npwmb
C            call wronjje(jobBes,qpg2(ig),E,rmax,1,1,lcutmx,wjjE)
C            ilm = 0
C            do  l = 0, lcuta(iat)
C              do  m = -l, l
C                ilm = ilm+1
C                rojp(ig,ilm,iat) = (-wjjE(l))*pjyl(ilm,ig)
C              enddo
C            enddo
C          enddo ! G vectors
CC#ifdefC DEBUG
CC          call info2(1,0,0,' sumcheck rojp site %i %2;18,6D',ib,sum(rojp(:,:,iat)))
CC#endif
C
CC     ... Setup for sgbp and sgpp
CC         Note : the coulomb integral with J can be performed analytically.
CC         But there is an error anyway because of sgpb is an integral of J with B
CC         so we do it numerically.
CC         On the other hand it would be possible for sigpp, i.e. <J | v | J>,
CC         to be evaluated analytically.
C
CC         Tabulate Bessel jG = j_l(sqrt(E)r) * r / (sqrt(E))**l
CC         Note: slightly different convention than Jbar: jG = Jbar/(2l+1)!! [check]
C          allocate(xi(nr,0:lcuta(iat)),y(nr),h(nr))
C          do  ig = 1, npwmb
C            call ropbes(rofi,qpg2(ig),lcuta(iat),y,h,xi,nr,1) ! MSM Standard def J/r^l
C            do  ir = 1, nr
C              r2s = rofi(ir)
C              do  l = 0, lcuta(iat)
CC               call besslr(qpg2(ig)*rofi(ir)**2,10,0,lcuta(iat),phi,psi)
CC               jG(ir,ig,l,iat) = phi(l) * rofi(ir)**(l+1)
C                jG(ir,ig,l) = xi(ir,l) * r2s    ! r*J(G^2)
C                r2s = r2s*rofi(ir)
C              enddo
C            enddo
C          enddo
CC#ifdefC DEBUG
CC          call info2(1,0,0,' sumcheck jg site %i %;18,6D',ib,sum(jG(nr,:,:)))
CC#endif
C          deallocate(xi,y,h)
C
CC     ... sgpb = Int dr r*J_l  Int dr' v(r,r') r'*B_l(r')
C
CC         This version is faster, but less accurate and versatile
CC         ilm = 0
CC         do  l = 0, lcuta(iat)
CC           offr = npb(l+(lcutmx+1)*(iat-1))
CC           nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
CC           allocate(sig(npwmb,nf))
CCC          Integral v(r,r') B(r') for each of nf functions
CC           do  n1 = 1, nf
CC             call rintgp(11,rofi,a,rkpr(1,l),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,JBint)
CC             call rintgp(11,rofi,a,rkmr(1,l),rprodx(1,offr+n1),-1d0,nrmx,nr,1,9,1,errmx,HBint)
CC             forall (ir = 1:nr) vB(ir,n1) = rkmr(ir,l)*(JBint(1)-JBint(ir)) + rkpr(ir,l)*HBint(ir)
CC             do  ig = 1, npwmb
CC               sig(ig,n1) = dot3(nr,vB(1,n1),jG(1,ig,l),rwgt)
CC             enddo
CC           enddo
CC           do  m = -l, l
CC             ilm = ilm+1
CC             do  n1 = 1, nf
CC             do  ig = 1, npwmb
CC               sgpb(ig,n1,ilm,iat) = dconjg(pjyl(ilm,ig)) * sig(ig,n1)/(2*l+1)*fpi
CC             enddo
CC             enddo
CC           enddo
CC           deallocate(sig)
CC         enddo
CC         call info2(1,0,0,' sumcheck sgpb site %i %2;18,6D',ib,sum(sgpb(:,:,:,iat)))
C
CC         Alternative integration of sgpb : Int dr r*B(r) Int dr' v(r,r')  r'*J(r')
CC         More integrals than previous, maybe 2x slower than fast form above,
CC         but perhaps more accurate since J is smooth.
CC         Also integrals can be used for sgpp; they are accumulated for both sgpb and sgpp
C
C          do  l = 0, lcuta(iat)
C            do  ig = 1, npwmb
CC             These integrals could be done analytically
C              call rintgp(11,rofi,a,rkpr(1,l,0),jG(1,ig,l),-1d0,nrmx,nr,1,9,1,errmx,JBint)
C              call rintgp(11,rofi,a,rkmr(1,l,0),jG(1,ig,l),-1d0,nrmx,nr,1,9,1,errmx,HBint)
CC             vjG(r) = Int dr' v(r,r') r'*J(r')
C              forall (ir = 1:nr) vjG(ir,ig,l) = rkmr(ir,l,0)*(JBint(1)-JBint(ir)) + rkpr(ir,l,0)*HBint(ir)
CC             Fold in radial integration weights
C              forall (ir = 1:nr) vjG(ir,ig,l) = vjG(ir,ig,l)*rwgt(ir)
C            enddo
C          enddo
C          do  l = 0, lcuta(iat)
C            offr = npb(l+(lcutmx+1)*(iat-1))
C            nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
C            allocate(sig(npwmb,nf))
CC           Replace the following with faster dgemm call
CC            do  n1 = 1, nf
CC              do  ig = 1, npwmb
CC                sig(ig,n1) = ddot(nr,vjG(1,ig,l),1,rprodx(1,offr+n1),1)/(2*l+1)*fpi
CC              enddo
CC            enddo
C            call dgemm('T','N',npwmb,nf,nr,1d0,vjG(1,1,l),nrmx,rprodx(1,offr+1),nrmx,0d0,sig,npwmb)
C            call dscal(npwmb*nf,fpi/(2*l+1),sig,1)
C
C            do  ilm = l*l+1, (l+1)**2
C              do  n1 = 1, nf
C              do  ig = 1, npwmb
C                sgpb(ig,n1,ilm,iat) = dconjg(pjyl(ilm,ig)) * sig(ig,n1)
C              enddo
C              enddo
C            enddo
C            deallocate(sig)
C          enddo
CC#ifdefC DEBUG
CC          call info2(1,0,0,' sumcheck sgpb site %i %2;18,6D',ib,sum(sgpb(:,:,:,iat)))
CC#endif
C
CC     ... Make fouvb = 4*pi/(G^2-E) * integral(J*B) ... note E is negative here
C          do  ig = 1, npwmb
C            ilm = 0
C            do  l = 0, lcuta(iat)
C              offr = npb(l+(lcutmx+1)*(iat-1))
C              nf = npb(l+1+(lcutmx+1)*(iat-1))-offr
C              allocate(sig(nf,1))
C              do  n1 = 1, nf
C                sig(n1,1) = dot3(nr,jG(1,ig,l),rprodx(1,offr+n1),rwgt)
C              enddo
C              do  m = -l, l
C                ilm = ilm+1
C                do  n1 = 1, nf
C                  fouvb(ig,n1,ilm,iat) = fpi/(qpg2(ig)-E)*dconjg(pjyl(ilm,ig))*sig(n1,1)
C                enddo
C              enddo
C              deallocate(sig)
C            enddo
C          enddo
CC#ifdefC DEBUG
CC          call info2(1,0,0,' sumcheck fouvb/10 site %i %2;18,6D',ib,sum(fouvb(:,:,:,iat))/10)
CC#endif
C
CC     --- Add sgpp + fouvp into vcoul (to save memory) ---
C          allocate(sig(npwmb,npwmb),wjj(npwmb,npwmb,0:lcuta(iat)))
CC     ... Wronskian wjj = W{qpg2(ig),qpg2(ig2)} at rmax for l=0:lcuta
C          call wronjje(jobBes,qpg2,qpg2,rmax,npwmb,npwmb,lcuta(iat),wjj)
C          do  l = 0, lcuta(iat)
C            call dgemm('T','N',npwmb,npwmb,nr,1d0,vjG(1,1,l),nrmx,jG(1,1,l),nrmx,0d0,sig,npwmb)
C            call dscal(npwmb*npwmb,fpi/(2*l+1),sig,1)
C
C            do  ig = 1, npwmb
C
CC             if (ib == 2 .and. ig == 2) call snot
C
C              do  ig2 = 1, ig
C                p12 = 0
C                do  ilm = l*l+1, (l+1)**2
CC                 p12 = p12 + dconjg(pjyl(ilm,ig)*phase(ig)) * pjyl(ilm,ig2)*phase(ig2)
C                  p12 = p12 + dconjg(pjyl(ilm,ig)) * pjyl(ilm,ig2)
C                enddo
C
C                wv12 = fpi/(qpg2(ig)-E)*wjj(ig,ig2,l)
C                wv21 = fpi/(qpg2(ig2)-E)*wjj(ig2,ig,l) ! Could exploit wjj(ig,ig2,l) = wjj(ig2,ig,l)
C
CC               Add onsite contribution (prop sig) and Fourier contribution (prop wv12 and wv21) to vcoul
C                vcoul(ntpb+ig,ntpb+ig2) = vcoul(ntpb+ig,ntpb+ig2) + (sig(ig,ig2) + wv12 + wv21) * p12
C              enddo ! ig2
C            enddo ! ig
C          enddo ! l
C          deallocate(sig,wjj)
C
C          deallocate(phase)
C
C        endif  ! mode1 > 0
C      enddo  ! Loop over sites
C      deallocate(rkpr,rkmr)
C      deallocate(JBint,Hbint,vB)
C
C      if (mode1 > 0) then
C
CC#ifdefC DEBUG
CC        call info2(1,0,0,' sumcheck onsite vcoul(PW-PW)/1000 %2;18,6D',sum(vcoul(ntpb:ntmb,ntpb:ntmb))/1000,2)
CC#ifdefC DEBUG2
CC        call zprm0('(1p9e22.12)')
CC        call zprm('vcoul (pp)',2,vcoul(ntpb+1,ntpb+1),ntmb,npwmb,npwmb)
CC#endifC
CC#endif
C
C        deallocate(qpg,yl,qpg2,pjyl)
C        deallocate(jG,vjG)
C
C      endif
C
C      call tcx('onsite-vc')
C
C      end
