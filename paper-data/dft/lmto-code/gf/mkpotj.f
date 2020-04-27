      subroutine mkpotj(nl,nbas,nsp,ipc,lmx,rmax,avw,indxsh,iopt,pp,
     .  zmv,eref,P)
C- Potential functions for free electrons
C-----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   lmx   :lmx(j) = maximum l for atom j (input)
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   iopt  :a compound of digits determining what is generated:
Ci         100s digit
Ci           1 make P for bare representation (alpha=0)
Ci           2 make P for gamma representation (alpha=gamma)
Ci          10s digit determines whether to make P or some derivative:
Ci           0 P <- potential function P
Ci           1 P <- P-dot
Ci           2 P <- -1/2 P-dotdot/P-dot
Ci           3 P <- sqrt(P-dot), choosing abs(sqrt(delta))
Ci           4 P <- sqrt(P-dot(bare)/P-dot(alpha)): then
Ci                  g^alpha_ij = P_i g^bare_ij P_j for i<>j
Ci           5 P <- scaling to convert diagonal part of g^bare
Ci                  to g^alpha
Ci          1s digit
Ci           0 second-order P
Ci           1 third-order P
Ci           2 first-order P
Ci   pp    :potential parameters (atomsr.f)
Ci   zmv   :complex energy inside sphere relative to background pot.
Ci   eref  :reference energy prescribing matching at rmax; see Remarks.
Co Outputs
Co   P:  potential function or some derivative, depending on iopt
Cr Remarks
Cr  *Potential functions connect sphere of constant potential
Cr   where wave function matches smoothly to Hankels, Bessels
Cr   of energy eref at rmx (eref=0 for classical ASA).
Cr
Cr  *Bare Potential function is defined for wave function phi as
Cr     W{K,phi} / W{J,phi}
Cr   For free electrons
Cr     phi = J / <J J>   normalization: <J J> = W{J,J-dot}
Cr  *Energy derivative from Gunnarsson PRB 27, 7144, 1983, Eq. A21
Cr     P-dot = -W{J,J-dot} (avw/2) W^-2{J,phi}
Cr  *For now we are connecting to strux with kappa=0.  So
Cr     P-dot = -W{J,J-dot} avw/2 W^-2 {J,J0}
Cr     -1/2 P-dot-dot/P-dot = W{J,J0-dot} W^-1{J,J0}
Cr  *Transformation to P^alpha:
Cr   P^-1(bare) = P^-1(alpha) => P(alpha) = P(bare) / (1-alpha*P(bare))
C-----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nl,nbas,nsp,ipc(nbas),lmx(nbas),indxsh(nl**2*nbas),iopt
      double precision pp(6,nl,nsp,*),rmax(nbas),avw
      double complex zmv,P(nl**2*nbas,nsp),eref
C Local parameters
      integer lmax,lmr,ibas,l,m,ic,isp,ipr,lgunit,modeP,irep
      double precision rmx,alp
      double complex phold(10)
C     double precision V0,G0,sG0
      double complex wkk(0:20),wkj(0:20),wjk(0:20),wjj0(0:20),wjj(0:20),
     .  wjjd(0:20),pfun,pdot,pfa

C --- Setup and printout  ---
      call getpr(ipr)
c      ipr=55
      irep    = mod(iopt/100,10)
      modeP   = mod(iopt/10,10)

      if (ipr >= 55) then
        call awrit4(' mkpotj: '//
     .    '%?#n==5#delta-##%-1j'//
     .    '%?#n==4#ratio-##%-1j'//
     .    '%?#(n==3|n==4)#sqrt-##%-1j'//
     .    '%?#n>=1#dot-##%-1j'//
     .    '%?#(n==2|n==5)#dot-##'//
     .    'P-%?#n==0#alpha##%-1j%?#n==1#bare##%-1j%?#n==2#gamma##,'//
     .    ' z-v= (%1;6d,%1;6d):',' ',80,lgunit(1),
     .    modeP,irep,zmv,dimag(zmv))
        if (nsp == 1) then
          print '(''  ib   P_l(ib) ...'')'
        else
          print '(''  ib spin  P_l(ib) ...'')'
        endif
      endif

C --- For each site, make P or related function ---
      do  4  isp = 1, nsp
        lmr = 0
        do  3  ibas = 1, nbas
          if (ipc(1) == 0) then
            rmx = rmax(ibas)
            lmax = lmx(ibas)
          else
            ic = ipc(ibas)
            rmx = rmax(ipc(ibas))
            lmax = lmx(ipc(ibas))
          endif
          call wrnhjz(zmv,zmv,rmx,lmax,avw,wkk,wkj,wjk,wjj,111)
          call wrnhjz(eref,zmv,rmx,lmax,avw,wkk,wkj,wjk,wjjd,111)
          call wrnhjz(eref,zmv,rmx,lmax,avw,wkk,wkj,wjk,wjj0,110)

          do  2  l = 0, nl-1
            alp = pp(6,l+1,isp,ic)
            if (irep == 1) alp = 0
            if (irep == 2) alp = pp(5,l+1,isp,ic)
            do  1  m = -l, l
              lmr = lmr + 1
              if (l > lmx(ic)) goto 1
              pfun = wkj(l)/wjj0(l)
              pdot = -wjj(l) * avw/2 / wjj0(l)**2
              pfa  = 1 - alp * pfun
C         ... P
              if (modeP == 0) then
                phold(l+1) = pfun / pfa
C         ... P-dot
              elseif (modeP == 1) then
                phold(l+1) = pdot / pfa**2
C         ... -1/2 P-dot-dot / P-dot
              elseif (modeP == 2) then
                phold(l+1) = wjjd(l) / wjj0(l) - alp*pdot/pfa
C         ... sqrt(P-dot)
              elseif (modeP == 3) then
                phold(l+1) = -sqrt(-wjj(l)) * sqrt(avw/2)/wjj0(l)/pfa
C         ... sqrt(P-dot(bare)/P-dot(alpha))
              elseif (modeP == 4) then
                phold(l+1) = pfa
C         ... delta_bare-alpha (-1/2 P-dot-dot / P-dot)
              elseif (modeP == 5) then
                phold(l+1) = alp*pfa
              else
                call rx('mkpotj: bad iopt')
              endif
              P(indxsh(lmr),isp) = phold(l+1)
    1       continue
    2     continue
          if (ipr >= 55 .and. nsp == 1) then
            print 333,  ibas, (phold(l), l=1,nl)
  333       format(i4,10(1x,2f10.5,2x))
          elseif (ipr >= 55 .and. nsp == 2) then
            print 334,  ibas, isp, (phold(l), l=1,nl)
  334       format(2i4,10(1x,2f10.5,1x))
          endif
    3   continue
    4 continue

      end
