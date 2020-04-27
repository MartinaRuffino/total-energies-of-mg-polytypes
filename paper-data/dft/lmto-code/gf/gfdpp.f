      subroutine gfdpp(mode,s_ham,s_pot,offH,ib,nbas,isp,ldp,dp,dpd,dpdd)
C- Change in potential functions from pot. perturbation at one site
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 lncol ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dpf ddpf dddpf
Cio    Passed to:  *
Ci Inputs
Ci   mode  :0 potential perturbation is a constant potential shift
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   ib    :site with potential perturbation
Ci   isp   :current spin channel (1 or 2)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ldp   :leading dimension of dp,dpd,dpdd
Co Outputs
Co   dp    :change in potential function from potential perturbation
Co         :In the noncollinear case dp is dimensioned dp(ldp,2,2)
Co         :and is a 2x2 spinor is returned
Co   dpd   :change in pdot from potential perturbation
Co         :dimensions are same as dp
Co   dpdd  :change in -1/2 P-dotdot/P-dot from potential perturbation
Co         :dimensions are same as dp
Cr Remarks
Cr   This routine produces changes in potential parameters corresponding
Cr   to a perturbation at site ib.  The change is hamiltonian-dependent.
Cr
Cr   ... 2nd generation ASA:
Cr   Potential parameters are (diagonal) P, Pdot, etc.
Cr   mode 1 (perturbation = constant potential shift)
Cr     dP/dV = -Pdot, etc.  Thus:
Cr     dPdot/dV = -P-dot-dot.  For this, use:
Cr       P-dot-dot = -2 * dpf * ddpf, dpf, ddpf gen from mkptfp.
Cr     d(-1/2 P-dotdot/P-dot)/dV  is made in mkptfp, mode 6.
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Mar 00 extended to include iwaves and hwaves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer mode,ib,isp,ldp,offH(n0H,nkap0,1),nbas
      double complex dp(ldp,*),dpd(ldp,*),dpdd(ldp,*)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Local parameters
      integer ldham(16),ldim,lgen3,lhdim,lidim,lncol,nspc,
     .  offpi,nlm0,iblk,nlmi,isum
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lhdim,ldham(3))
      equivalence (nspc,ldham(4))

      lgen3 = s_ham%lgen3
      lncol = s_ham%lncol
      ldham = s_ham%ldham

C --- 3nd generation NMTO --
      if (lgen3 /= 0) then
        call rx('gfdpp not ready for 3rd generation')

C --- 2nd generation ASA ---
      else

C   ... dP for lower, intermediate, higher blocks
        nlm0 = 0
        do  10  iblk = 1, 3
          offpi = offH(iblk,1,ib)
          nlmi  = offH(iblk,1,ib+1) - offpi
          if (nlm0+nlmi > ldp) return
          offpi = offpi
     .          + isum(iblk,offH(1,1,nbas+1),1)-offH(iblk,1,nbas+1)
          call gfdpp2(mode,nlmi,offpi,isp,nspc,lhdim,s_pot%dpf,
     .      s_pot%ddpf,s_pot%dddpf,ldp,dp(1+nlm0,1),dpd(1+nlm0,1),
     .      dpdd(1+nlm0,1))
          nlm0 = nlm0+nlmi
   10   continue
      endif

      end

      subroutine gfdpp2(mode,nlmi,offpi,isp,nspc,lhdim,dpf,ddpf,dddpf,
     .  ldp,dp,dpd,dpdd)
C- Change in 2nd gen. pot. functions from pot. perturbation at one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  0 potential perturbation is a constant potential shift
Ci   nlmi  :dimension of hamiltonian at site
Ci   offpi :offset to first entry in pf,dpf,ddpf for site of pot. shift
Ci   dpf   :1st energy derivative of potential function P-dot (mkptfp.f)
Ci   ddpf  :scaled 2nd energy derivative (see Remarks)
Ci   dddpf :scaled 3rd energy derivative (see Remarks)
Co Outputs
Co   dp    :change in potential function from potential perturbation
Co   dpd   :change in pdot from potential perturbation
Co   dpdd  :change in -1/2 P-dotdot/P-dot from potential perturbation
Cl Local variables
Cl   jsp   :index to spin-diagonal part of gii
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   ksp   :isp in the collinear case, and jsp in the noncollinear
Cl         :Thus, for the spinor parts of gii and gd:
Cl         :ksp = index to dpf etc for current spin
Cl         :(jsp,jsp) = indices to diag part of dp etc for current spin
Cr Remarks
Cr   mode 1 (perturbation is constant potential shift) :
Cr     dP/dV = -Pdot, etc.  Thus:
Cr     dPdot/dV = -P-dot-dot.  For this, use:
Cr       P-dot-dot = -2 * dpf * ddpf, dpf, ddpf gen from mkptfp.
Cr     d(-1/2 P-dotdot/P-dot)/dV  is made in mkptfp, mode 6.
Cr     NB:  modes 2,6 in mkptfp construct these derivatives.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer mode,nlmi,lhdim,offpi,isp,nspc,ldp
      double complex dpf(lhdim,2),ddpf(lhdim,2),dddpf(lhdim,2)
      double complex dp(ldp,nspc,nspc),dpd(ldp,nspc,nspc),
     .               dpdd(ldp,nspc,nspc)
      integer i,jsp,ksp

      call sanrg(.true.,mode,0,0,'gfdpg2','mode')

      do  10  jsp = 1, nspc
        ksp = max(jsp,isp)

C   --- Case perturbation is a constant potential shift ---
        do  12  i = 1, nlmi
          dp(i,jsp,jsp) = -dpf(i+offpi,ksp)
          dpd(i,jsp,jsp) = -2 * dpf(i+offpi,ksp) * ddpf(i+offpi,ksp)
          dpdd(i,jsp,jsp) = dddpf(i+offpi,ksp)
          if (nspc == 2) then
            dp(i,jsp,3-jsp) = 0
            dpd(i,jsp,3-jsp) = 0
            dpdd(i,jsp,3-jsp) = 0
          endif
   12   continue

   10 continue

      end
