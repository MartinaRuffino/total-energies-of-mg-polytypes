      subroutine asxsgm(opt,aintra,iqin,plat,nl,nbas,offH,indxsh,
     .  nkp,wtkp,lshft,ipq,igstar,ldim,lpdim,eval,nbmx,evec,bas,istab,g,
     .  ag,efermi,wscrk,wscrl,z,dmat,wtr,sigma)
C- ASA screened exchange potential from nonlocal, static linear response
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt  :
Ci          1s digit not used now
Ci         10s digit
Ci           0: use real spherical harmonics
Ci           1: use complex spherical harmonics
Ci   eval,evec
Ci   wscrk     W(q)
Ci   wscrl     W, local part
Ci   wtr       work array to hold W(rotated q)
Ci   z         work array to hold evec(rotated q)
Ci   dmat      work array to hold dmat = sum_nu=occ z_RL1,nu z+_RL2,nu
Co Outputs
Co   sigma
Cr Remarks
Cr   Given q-dependent W = screened coulomb potential, and local part
Cr
Cr      X                        i      i
Cr     V_(q)      = sum(k)      z  (k) z  (k)  [W  (k-q)  - Wl_Rl1 ]
Cr       Rlm1,Rlm2     i=occ     Rlm1   Rlm2     Rl1,Rl2
Cr
Cr  To do:
Cr   replace indxsh with offH in wtrans
C ----------------------------------------------------------------------
      implicit none
      logical lshft(3)
      integer opt,nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer iqin,nbmx,nbas,nl,ldim,nkp,lpdim,istab(nbas,*),
     .  offH(n0H,nkap0,*),indxsh(*),ipq(*),igstar(0:*)
      double precision plat(3,3),bas(3,*),eval(nbmx,nkp),efermi,
     . g(3,3,*),ag(3,*),evec(ldim,ldim,2,nkp),wtkp(1),z(ldim,ldim,2),
     .  dmat(ldim,ldim,2),wtr(lpdim,lpdim,2),wscrk(lpdim,lpdim,2,nkp),
     .  wscrl(lpdim,lpdim,2),sigma(ldim,ldim,2)
      logical aintra
C Local variables
      integer iq0,j1,j2,j3,iq1,iq2,i1,i2,i3,ii1,ii2,ii3,jpq0,ipq1,
     .  jpq1,ipq2,jpq2,nu1,iprint,k,nkxyz(3),is(3),
     .  ifac(3),jj1,jj2,jj3,nuocc
      double precision wght1,rb(3,3),qb(3,3),q1(3),q2(3),qk
C      integer nlx
C      parameter (nlx=4)
C      double precision rmat(nlx**4)
      double precision rmat(nl**4*2)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C --- Setup ---
      call tcn('asxsgm')
C     if (nl > nlx) call rx('asxsgm: increase nlx')
      k = iabs(igstar(0))
      call cmpi123(0,nkxyz(1),nkxyz(2),nkxyz(3),k)
C ... Make is,ifac,qb,wght1
      call pshpr(1)
      call bzmsh0(plat,lshft,0,nkxyz(1),nkxyz(2),nkxyz(3),is,ifac,rb,qb)
      call poppr
      wght1 = 2d0/nkxyz(3)/nkxyz(2)/nkxyz(1)
      call dpzero(sigma,ldim**2*2)

C --- identify iqin ---
      iq0 = 0
      do  10  j3 = 1, nkxyz(3)
      do  10  j2 = 1, nkxyz(2)
      do  10  j1 = 1, nkxyz(1)
        iq0 = iq0+1
        if (iqin /= ipq(iq0)) goto 10
        jpq0 = igstar(iq0)
        if (jpq0 /= 1) goto 10

C --- k-integration for vsxnl(q) ---
        iq1 = 0
        do  20  i3 = 1, nkxyz(3)
        ii3 = i3+j3-1
        if (ii3 > nkxyz(3)) ii3 = ii3-nkxyz(3)
        do  20  i2 = 1, nkxyz(2)
        ii2 = i2+j2-1
        if (ii2 > nkxyz(2)) ii2 = ii2-nkxyz(2)
        do  20  i1 = 1, nkxyz(1)
        ii1 = i1+j1-1
        if (ii1 > nkxyz(1)) ii1 = ii1-nkxyz(1)

C   --- Potential for q1=k; scale by wght1 ---
        iq1 = iq1+1
        ipq1 = ipq(iq1)
        jpq1 = igstar(iq1)
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)
C       call yprm('w',2,wscrk(1,1,1,ipq1),lpdim*lpdim,lpdim,lpdim,lpdim)
        call wtrans(ldim,lpdim,nl,nbas,bas,indxsh,istab(1,jpq1),
     .    g(1,1,jpq1),ag(1,jpq1),q1,aintra,wscrk(1,1,1,ipq1),wtr)
C       call yprm('wtrans',2,wtr,lpdim*lpdim,lpdim,lpdim,lpdim)
        call dscal(lpdim**2*2,wght1,wtr,1)

C   --- Wave function for q2=k+q ---
        iq2 = ii1+nkxyz(1)*((ii2-1)+nkxyz(2)*(ii3-1))
        ipq2 = ipq(iq2)
        jpq2 = igstar(iq2)
        q2(1) = qk(1,ii1,ii2,ii3)
        q2(2) = qk(2,ii1,ii2,ii3)
        q2(3) = qk(3,ii1,ii2,ii3)
        call rotwf(opt,nl,nbas,1,bas,offH,indxsh,istab(1,jpq2),g(1,1,jpq2)
     .    ,ag(1,jpq2),q2,rmat,0,ldim,ldim,ldim,evec(1,1,1,ipq2),z)

C   --- Count number of occupied bands ---
        do  26  nu1 = 1, ldim
        nuocc = nu1-1
   26   if(eval(nu1,ipq1) > efermi) goto 27
   27   continue
C   --- dmat(k+q) = sum_nu=occ z_RL1,nu (z+)_RL2,nu ---
        call yygemm('N4','C',ldim,ldim,nuocc,1d0,
     .    z,z(1,1,2),ldim,z,z(1,1,2),ldim,0d0,dmat,dmat(1,1,2),ldim)

C       call yprm('evec',2,evec(1,1,1,ipq2),ldim*ldim,ldim,ldim,nuocc)
C       call yprm('z',2,z,ldim*ldim,ldim,ldim,nuocc)
C       call yprm('z.z+',2,dmat,ldim*ldim,ldim,ldim,ldim)

C   --- Accumulate sum_nu=occ z_RL1,nu (z+)_RL2,nu W_Rl1,Rl2 into Sigma
        call psxsgm(aintra,nbas,indxsh,ldim,lpdim,wscrl,wght1,dmat,wtr,
     .    sigma)

C      call yprm('sigm',2,sigma,ldim*ldim,ldim,ldim,ldim)
   20 continue
   10 continue

      if (iprint() >= 90) then
C       write(*,*) ' --- Re(sigma) -----'
C       write(*,123) ((sigma(i,j,1),j=1,10),i=1,10)
C       write(*,*) ' --- Im(sigma) -----'
C       write(*,123) ((sigma(i,j,2),j=1,10),i=1,10)
        call yprm('sigma',2,sigma,ldim*ldim,ldim,ldim,ldim)
      endif
C 123 format(10f8.3)
      call tcx('asxsgm')
      end

      subroutine psxsgm(aintra,nbas,indxsh,ldim,lpdim,wscrl,
     .  wght1,dmat,wtr,sigma)
C- Kernel called by asxsgm: accumulate peice of sigma
C ----------------------------------------------------------------------
Ci Inputs
Ci   aintra
Ci   nbas  :size of basis
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lpdim
Ci   wscrl
Ci   wght1
Ci   dmat  :density matrix at k+q
Ci   wtr   :screened coulomb interaction at k
Co Outputs
Co   sigma :dmat(k+q)*wtr(k) is added to sigma
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   07 Nov 07
C ----------------------------------------------------------------------
      implicit none
      integer nbas,ldim,lpdim,indxsh(*)
      double precision dmat(ldim,ldim,2),wtr(lpdim,lpdim,2),
     .  wscrl(lpdim,lpdim,2),sigma(ldim,ldim,2),wght1
      logical aintra
C Local variables
      integer mxorb,nglob,lmxa,
     .  nl1,ib1,off1,nlm1,l1,m1,lmr1,
     .  nl2,ib2,off2,nlm2,l2,m2,lmr2
      double precision xxr,xxi

      mxorb = nglob('mxorb')
      lmxa = nglob('nl')-1

C --- Loop over orbitals ib2, extracting Rl(ib2) for wtr ---
      nl2 = 0
      do  40  ib2 = 1, nbas
        lmr2 = mxorb*(ib2-1)
        do  42  l2 = 0, lmxa
          off2 = indxsh(lmr2+1) - 1
          lmr2 = lmr2 + 2*l2+1
          if (off2 >= ldim) goto 42
          nlm2 = 2*l2+1
          nl2 = nl2+1
C       --- Loop over orbitals ib1, extracting Rl(ib1) for wtr ---
            nl1 = 0
            do  50  ib1 = 1, nbas
              lmr1 = mxorb*(ib1-1)
              do  52  l1 = 0, lmxa
                off1 = indxsh(lmr1+1) - 1
                lmr1 = lmr1 + 2*l1+1
                if (off1 >= ldim) goto 52
                nlm1 = 2*l1+1
                nl1 = nl1+1
                if (aintra) then
                  xxr = wtr(nl1,nl2,1)
                  xxi = wtr(nl1,nl2,2)
                else
                  xxr = wtr(ib1,ib2,1)
                  xxi = wtr(ib1,ib2,2)
                endif
C           ... Accumulate z1_Rlm1 W_Rl1,Rl2 (z+)_Rlm2  into Sigma
                do  58  m2 = 1, nlm2
                do  58  m1 = 1, nlm1
                  sigma(m1+off1,m2+off2,1) = sigma(m1+off1,m2+off2,1) +
     .              dmat(m1+off1,m2+off2,1)*xxr -
     .              dmat(m1+off1,m2+off2,2)*xxi
                  sigma(m1+off1,m2+off2,2) = sigma(m1+off1,m2+off2,2) +
     .              dmat(m1+off1,m2+off2,2)*xxr +
     .              dmat(m1+off1,m2+off2,1)*xxi
   58           continue
   52         continue
   50       continue

C       --- Subtract (z2+) z2 wscrl from diagonal sigma ---
          do  45  m2 = 1, nlm2
            if (aintra) then
              xxr = -wscrl(nl2,nl2,1)*wght1
              xxi = -wscrl(nl2,nl2,2)*wght1
            else
              xxr = -wscrl(ib2,ib2,1)*wght1
              xxi = -wscrl(ib2,ib2,2)*wght1
            endif
            sigma(m2+off2,m2+off2,1) = sigma(m2+off2,m2+off2,1) +
     .        dmat(m2+off2,m2+off2,1)*xxr - dmat(m2+off2,m2+off2,2)*xxi
            sigma(m2+off2,m2+off2,2) = sigma(m2+off2,m2+off2,2) +
     .        dmat(m2+off2,m2+off2,2)*xxr + dmat(m2+off2,m2+off2,1)*xxi
   45     continue
   42   continue
   40 continue
        end
