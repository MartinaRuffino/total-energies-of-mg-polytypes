      subroutine gfwscrbz(aintra,plat,nl,nbas,offH,indxsh,nk1,nk2,nk3,
     .  nkp,lshft,ipq,igstar,ldim,lpdim,nlmx,nlma,
     .  bas,istab,g,ag,wscrk,wscrl,wscrmat)
C- ASA-GF screened Coulomb potential in full BZ rolled to the full
C- extent of H(ldim,ldim)
C ----------------------------------------------------------------------
Ci Inputs
Ci  aintra:logical, if true, the intra-site term of the Coulomb interaction
Ci         given by the second derivative of the total energy witgh respect
Ci         to occupation number is added to the inter-atomic bare(uncreened)
Ci         Coulomb interaction given by the Madelung matrix. It is supposed
Ci         to be true in gf implementation of the screened exchange.
Ci  plat  :lattice vectors.
Ci  nl    :l max
Ci  nbas  :no. of atoms in the basis.
Ci  offH  :Offsets to the Hamiltonian matrix.
Ci  indxsh:permutation ordering orbitals in l+i+h blocks.
Ci  nk1,nk2,nk3:no.of divisions in the BZ mesh.
Ci  nkp   :no.of irreducible q-points
Ci  lshft :shift of the BZ with respect to gamma point.
Ci  ipq   :ipq(i1,i2,i3) point to the irred. q-points into which mesh point
Ci         (i1,i2,i3) is mapped (bzmesh.)
Ci  igstar:info tomao an irred. q-point to its original point (bzmesh).
Ci  ldim  : dimension of the Hamiltonin or Green function.
Ci  lpdim:leading dimension of wscrk.
Ci  nlmx :current index indicating which slice of the wscrmat to be handled.
Ci  nlma :added to nlmx and indicates the slice handled.
Ci  bas  :basis vectors.
Ci  istab:table of site permutationfor each group operation(mksym.f symtbl.f).
Ci  g    :point group operations.
Ci  ag   :translation part of space group.
Ci
Ci
Ci   wscrk     W(q)
Ci   wscrl     W, local part
Ci   wtr       work array to hold W(rotated q)
Ci
Ci
Co Outputs
Co   wscrmat
Cr Remarks
Cr   W(k) is transformed to W(-k) to be used in gfgw (T. Sandu, 27 Sept 04)
Cr
Cr   Given q-dependent(irreducible) W = screened coulomb potential,
Cr   and local part
Cr   one obtains W(k=in full Brillouin) rolled over m-indices
Cr
Cr     W_(k)      = [W  (q*=k)  - Wl_RL1 ]
Cr       RLm1,RLm2    RL1,RL2
Cr
Cr  To do:
Cr   replace indxsh with offH in wtrans
Cu Updates
Cu   20 Sep 04 (T. Sandu) first created
C ----------------------------------------------------------------------
      implicit none
      logical lshft(3)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nlmx,nlma,nbas,nl,ldim,nkp,lpdim,istab(nbas,1),
     .  offH(n0H,nkap0,1),indxsh(1),ipq(1),igstar(0:*)
      double precision plat(3,3),bas(3,*),g(3,3,*),ag(3,*),
     .  wtr(lpdim,lpdim,2),wscrk(lpdim,lpdim,2,nkp),
     .  wscrl(lpdim,lpdim,2),wtr1(ldim,ldim,2)
      logical aintra
C Local variables
      integer iq1,i1,i2,i3,ipq1,jpq1,iprint,k,
     .  nkxyz(3),is(3),ifac(3),jj1,jj2,jj3,nk1,nk2,nk3
      double precision rb(3,3),qb(3,3),q1(3),qk
      double precision wscrmat(2,nk1,nk2,nk3,nlma,ldim)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C --- Setup ---
      call tcn('gfwscrbz')
C     if (nl > nlx) call rx('asxsgm: increase nlx')
C      k = iabs(igstar(0))
C      call cmpi123(0,nkxyz(1),nkxyz(2),nkxyz(3),k)
      nkxyz(3) = nk3
      nkxyz(2) = nk2
      nkxyz(1) = nk1
C ... Make is,ifac,qb
      call pshpr(1)
      call bzmsh0(plat,lshft,0,nkxyz(1),nkxyz(2),nkxyz(3),is,ifac,rb,qb)
      call poppr
C     wght1 = 2d0/nkxyz(3)/nkxyz(2)/nkxyz(1)

C --- k-space expansion for vsxnl(q) ---
        iq1 = 0
        do  20  i3 = 1, nkxyz(3)
        do  20  i2 = 1, nkxyz(2)
        do  20  i1 = 1, nkxyz(1)

C   --- Potential for q1=k ---
        iq1 = iq1+1
        ipq1 = ipq(iq1)
        jpq1 = igstar(iq1)
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)
c        call yprm('gf',2,wscrk(1,1,1,ipq1),lpdim*lpdim,lpdim,lpdim,lpdim)
        call wtrans(ldim,lpdim,nl,nbas,bas,indxsh,istab(1,jpq1),
     .   g(1,1,jpq1),ag(1,jpq1),q1,aintra,wscrk(1,1,1,ipq1),wtr)
c        call yprm('wtransgf',2,wtr,lpdim*lpdim,lpdim,lpdim,lpdim)

C   --- completing (i1,i2,i3)-k-point and Rlm1,Rlm2 dependence of W
c        print*,'i1,i2,i3',i1,i2,i3
        call wsrlm(aintra,nbas,indxsh,ldim,nlma,nlmx,lpdim,wscrl,
     .   i1,i2,i3,nkxyz(1),nkxyz(2),nkxyz(3),wtr,wtr1,wscrmat)

c      print*,'i1 i2 i3',i1,i2,i3
c      call yprm('wscr',2,wtr1,ldim*ldim,ldim,ldim,ldim)
      if (iprint() >= 90) then
        print*,'i1 i2 i3',i1,i2,i3
        call yprm('wscr',2,wtr1,ldim*ldim,ldim,ldim,ldim)
      endif

   20 continue

      call tcx('gfwscrbz')
      end

      subroutine wsrlm(aintra,nbas,indxsh,ldim,nlma,nlmx,lpdim,wscrl,
     .  i1,i2,i3,nk1,nk2,nk3,wtr,wtr1,wwsij)
C- Kernel called by gfwscrbz
C-w static is also transformed from k to -k in order to be handled properly by
C gfgw subroutine
      implicit none
      integer nbas,ldim,lpdim,nlma,indxsh(*)
      integer i1,i2,i3,nk1,nk2,nk3,nlmx
      double precision wtr(lpdim,lpdim,2),
     .  wscrl(lpdim,lpdim,2),wtr1(ldim,ldim,2)
C     .  wtr2(ldim,ldim,2),
      double precision wwsij(2,nk1,nk2,nk3,nlma,ldim)
      logical aintra
C Local variables
      integer mxorb,nglob,lmxa,
     .  nl1,ib1,off1,nlm1,l1,m1,lmr1,
     .  nl2,ib2,off2,nlm2,l2,m2,lmr2
      double precision xxr,xxi
      integer ii1,ii2,ii3

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
C           ... Rolling W_Rl1,Rl2 over_Rlm1,Rlm2 W_Rl1,Rl2
                do  58  m2 = 1, nlm2
                do  58  m1 = 1, nlm1
                  wtr1(m1+off1,m2+off2,1) = xxr
                  wtr1(m1+off1,m2+off2,2) = xxi
C                 wtr2(m1+off1,m2+off2,1)= xxr
C                 wtr2(m1+off1,m2+off2,2)= xxi


   58           continue
   52         continue
   50       continue

C       --- Subtract wscrl from diagonal wscrk ---
c          do  45  m2 = 1, nlm2
c            if (aintra) then
c              xxr = -wscrl(nl2,nl2,1)
c              xxi = -wscrl(nl2,nl2,2)
c            else
c              xxr = -wscrl(ib2,ib2,1)
c              xxi = -wscrl(ib2,ib2,2)
c            endif
c            wtr1(m2+off2,m2+off2,1) = wtr1(m2+off2,m2+off2,1) + xxr
c            wtr1(m2+off2,m2+off2,2) = wtr1(m2+off2,m2+off2,2) + xxi
c   45     continue
   42   continue
   40 continue
c      call yprm('wtr1',2,wtr1,ldim*ldim,ldim,ldim,ldim)
           do 400  m2 = 1, ldim
            do 400  m1 = 1,nlma
             ii3 = nk3 + 2 -i3
             if (i3 == 1) ii3 = i3
             ii2 = nk2 + 2 -i2
             if (i2 == 1) ii2 = i2
             ii1 = nk1 + 2 -i1
             if (i1 == 1) ii1 = i1

             wwsij(1,ii1,ii2,ii3,m1,m2) = wtr1(nlmx+m1,m2,1)
             wwsij(2,ii1,ii2,ii3,m1,m2) = wtr1(nlmx+m1,m2,2)

 400       continue


c         print *,'i1 i2 i3',i1,i2,i3
c        call yprm('wtr2',2,wtr2,ldim*ldim,ldim,ldim,ldim)
        end
