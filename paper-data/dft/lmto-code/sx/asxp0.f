      subroutine asxp0(opt,aintra,iqin,isp,plat,nl,nsp,nbas,offH,indxsh,
     .  nkp,lshft,ipq,igstar,ldim,lpdim,eval,nbmx,evec,zdel,bas,istab,g,
     .  ag,efermi,h1,h2,h3,p00)
C- ASA nonlocal, static linear response P0 for noninteracting particles
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt  :
Ci          1s digit not used now
Ci         10s digit
Ci           0: use real spherical harmonics
Ci           1: use complex spherical harmonics
Ci   aintra:F if to contract over l and m
Ci         :T if to contract over l
Ci   iqin  :make P0 for qp index iqin
Ci   isp   :current spin channel (1 or 2)
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   w(offH):Offsets to hamiltonian matrix (makidx.f)
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   lshft :uniform shift of k-point mesh
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   igstar:contains info needed to map an irreducible qp to its
Ci          original point (bzmesh.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lpdim :dimension of P0 matrix
Ci   eval  :eigenvalues corresponding to evecs.  NB: for spin pol case,
Ci         :eval has spin index, but evecs do not.
Ci   nbmx  :leading dimension of eval
Ci   evec  :eigenvectors for all irreducible qp
Ci   zdel  :energy for which to make P0
Ci   bas   :basis vectors, in units of alat
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   efermi:Fermi energy
Ci   h1,h2,h3 are work arrays of dimension of the hamiltonian
Co Outputs
Co   p00
Cr Remarks
Cr                     i      j       i     j
Cr  P_ab(q) = sum_ijk z (k) z+ (k+q) z+(k) z (k+q) times
Cr                     a      a       b     b
Cr
Cr              (f_i(k) - f_j(k+q))/(e_i(k) - e_j(k+q) + zdel)
Cr
Cr              where z^i_a is coefficient to orbital a of evec i
Cr
Cr  The numerator is grouped as
Cr      i      j       i     j        ( i     i   ) (  j       j     )
Cr     z (k) z+ (k+q) z+(k) z (k+q) = (z (k) z+(k)) (z+ (k+q) z (k+q))
Cr      a      a       b     b        ( a     b   ) (  a       b     )
Cr  To do:
Cr   So far an insulator assumed.
C ----------------------------------------------------------------------
      implicit none
      logical lshft(3)
      integer opt,nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer iqin,nbas,nl,isp,nsp,ldim,nkp,lpdim,offH(n0H,nkap0,*),
     .  indxsh(1),istab(nbas,1),nbmx,ipq(1),igstar(0:*)
      double precision plat(3,3),bas(3,*),eval(nbmx,nsp,nkp),efermi,
     .  g(3,3,*),ag(3,*),evec(ldim,ldim,2,nkp),h1(ldim,ldim,2),
     .  h2(ldim,ldim,2),h3(ldim,ldim,2),p00(lpdim,lpdim,2)
      double complex zdel
      logical aintra
C Local variables
      integer iq0,j1,j2,j3,iq1,iq2,i1,i2,i3,ipq0,ii1,ii2,ii3,k,
     .  jj1,jj2,jj3,jpq0,ipq1,jpq1,ipq2,jpq2,nu1,nu2,nl1,getdig,
     .  nkxyz(3),is(3),ifac(3),i,norbi,ili,li,nlmi,offi,ib,l1,
     .  l2,numax1,numax2,nglob,mxorb,lmr
      double complex z1,z2,ez
      double precision wght1,rb(3,3),qb(3,3),qlat(3,3),q1(3),q2(3),qk
C     integer nlx
C     parameter (nlx=4)
C     double precision rmat(nlx**4*2)
      double precision rmat(nl**4*2)

C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C --- Setup ---
*     COUNT=0
      call tcn('asxp0')
      ipq0 = iqin
C ... Unpack nkxyz
      k = iabs(igstar(0))
      call cmpi123(0,nkxyz(1),nkxyz(2),nkxyz(3),k)
C ... Make is,ifac,qb,qlat,wght1
      call pshpr(1)
      call bzmsh0(plat,lshft,0,nkxyz(1),nkxyz(2),nkxyz(3),is,ifac,rb,qb)
      call poppr
      call dinv33(plat,1,qlat,wght1)
      wght1 = 2d0/nkxyz(3)/nkxyz(2)/nkxyz(1)
C     if (nl > nlx) call rx('asxp0: increase nlx')
      mxorb = nglob('mxorb')

C --- Loop over q-points to get P0(q,zdel) ---
      iq0 = 0
      do  10  j3 = 1, nkxyz(3)
      do  10  j2 = 1, nkxyz(2)
      do  10  j1 = 1, nkxyz(1)
        iq0 = iq0+1
        if (ipq0 /= ipq(iq0)) goto 10
        jpq0 = igstar(iq0)
        if (jpq0 /= 1) goto 10

C   --- k-integration ---
        call dpzero(p00,2*lpdim**2)
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

C   --- Wave function for q2=k+q ---
        iq2 = ii1+nkxyz(1)*((ii2-1)+nkxyz(2)*(ii3-1))
        ipq2 = ipq(iq2)
        jpq2 = igstar(iq2)
        q2(1) = qk(1,ii1,ii2,ii3)
        q2(2) = qk(2,ii1,ii2,ii3)
        q2(3) = qk(3,ii1,ii2,ii3)
C  ...  The following valid for unshifted mesh only ...
C       qx(1) = (Ii1-1)*QB(1,1) + (Ii2-1)*QB(1,2) + (Ii3-1)*QB(1,3)
C       qx(2) = (Ii1-1)*QB(2,1) + (Ii2-1)*QB(2,2) + (Ii3-1)*QB(2,3)
C       qx(3) = (Ii1-1)*QB(3,1) + (Ii2-1)*QB(3,2) + (Ii3-1)*QB(3,3)
C       call daxpy(3,-1d0,q2,1,qx,1)
C       if (ddot(3,qx,1,qx,1) > 1d-12) stop 'bug in asxp0'

*        COUNT = count+1
*        print *, 'count a',count,ipq2
C        if (i1 == 1 .and. i2 == 2 .and. i3 == 1) then
C          print *, 'make h2'
C          call yprm('evec',2,evec(1,10,1,ipq2),ldim*ldim,ldim,9,1)
C        endif
        call rotwf(opt,nl,nbas,1,bas,offH,indxsh,istab(1,jpq2),g(1,1,jpq2)
     .    ,ag(1,jpq2),q2,rmat,0,ldim,ldim,ldim,evec(1,1,1,ipq2),h2)
C        if (i1 == 1 .and. i2 == 2 .and. i3 == 1) then
C          call yprm('h2',2,h2(1,10,1),ldim*ldim,ldim,9,1)
C        endif


C   --- Wave function for q1=k ---
        iq1 = iq1+1
        ipq1 = ipq(iq1)
        jpq1 = igstar(iq1)
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

*        COUNT = count+1
*        print *, 'count b',count
        call rotwf(opt,nl,nbas,1,bas,offH,indxsh,istab(1,jpq1),g(1,1,jpq1)
     .    ,ag(1,jpq1),q1,rmat,0,ldim,ldim,ldim,evec(1,1,1,ipq1),h1)

C   --- Loop over nu1 = occ bands, nu2 = unocc bands ---
C   ... NB: this should be reorganized to vectorize ib loop
        do  26  nu1 = 1, ldim
        numax1 = nu1-1
   26   if(eval(nu1,isp,ipq1) > efermi) goto 27
   27   continue
        do  28  nu2 = 1, ldim
        numax2 = nu2-1
   28   if(eval(nu2,isp,ipq2) > efermi) goto 29
   29   continue
        do  30  nu1 = 1, numax1

C     --- Loop over atoms and accumulate z+ z contracted over m ---
C         to make c^nu2_i = z+^nu1_i z^nu2_i for each orbital i
C         Index nu1 is suppressed: info retained only for current nu1
C         Also c^nu2_i is stored only for orbitals i corresponding to
C         site ib; additionally the i index is contracted over m.
C         C^nu2_i is stored in h3(nu2,i)
C         z^nu1 and z^nu2 are stored in h1 and h2.
C         Loop needs to be optimized and h3 eliminated
          nl1 = 0
          do  40  ib = 1, nbas
            lmr = mxorb*(ib-1)
            norbi = 0
            do  42  li = 0, nl-1
              offi = indxsh(lmr+1) - 1
              lmr = lmr + 2*li+1
              if (offi >= ldim) goto 42
              nl1 = nl1+1
C             Sanity check
C             if (nl1 > lpdim) call rx('bad indexing in asxp0')
              norbi = norbi+1
              nlmi = 2*li+1
              do  43  nu2 = numax2+1, ldim
              h3(nu2,nl1,1) = 0
   43         h3(nu2,nl1,2) = 0
              do  44 nu2 = numax2+1, ldim
                do  45  i = 1, nlmi
                  h3(nu2,nl1,1) = h3(nu2,nl1,1)
     .              + h1(offi+i,nu1,1)*h2(offi+i,nu2,1)
     .              + h1(offi+i,nu1,2)*h2(offi+i,nu2,2)
                  h3(nu2,nl1,2) = h3(nu2,nl1,2)
     .              - h1(offi+i,nu1,1)*h2(offi+i,nu2,2)
     .              + h1(offi+i,nu1,2)*h2(offi+i,nu2,1)
   45           continue
   44         continue
   42       continue

C       ... If contract over l also  ...
            if (.not. aintra) then
              do  46  nu2 = numax2+1, ldim
              h3(nu2,ib,1) = h3(nu2,nl1-norbi+1,1)
              h3(nu2,ib,2) = h3(nu2,nl1-norbi+1,2)
              do  46  ili = 2, norbi
                h3(nu2,ib,1) = h3(nu2,ib,1) + h3(nu2,nl1-norbi+ili,1)
                h3(nu2,ib,2) = h3(nu2,ib,2) + h3(nu2,nl1-norbi+ili,2)
   46         continue
            endif
   40   continue

C   --- qp contr. to integral c+_i c_j / (e^nu1 - e^nu2 + zdel) ---
C        print *, 'nu1,nu2=',nu1,nu2
C        call yprm('h1',2,h1(1,nu1,1),ldim*ldim,ldim,ldim,1)
C        call yprm('h2',2,h2(1,nu2,1),ldim*ldim,ldim,ldim,1)
        do  50  nu2 = numax2+1, ldim
        ez = 2/(eval(nu1,isp,ipq1)-eval(nu2,isp,ipq2)+zdel)
        do  50  l2 = 1, lpdim
          z2 = dcmplx(h3(nu2,l2,1),-h3(nu2,l2,2))*ez*wght1
          do  52  l1 = 1, lpdim
            z1 = z2*dcmplx(h3(nu2,l1,1),h3(nu2,l1,2))

C            if (i1 == 2 .and. i2 == 1 .and. i3 == 1 .or.
C     .          i1 == 1 .and. i2 == 2 .and. i3 == 1) then
C
C            if (l1 == 1 .and. l2 == 2 .and. nu2 == 10) then
C              print 375, i1,i2,i3, j1,j2,j3, nu1, nu2, z1,
C     .          p00(l1,l2,1),p00(l1,l2,2)
C  375         format(3i3,1x,3i3,1x,2i3,1x,2f13.7,1x,2f13.7)
C            endif
C            endif

            p00(l1,l2,1) = p00(l1,l2,1) + dble(z1)
            p00(l1,l2,2) = p00(l1,l2,2) + dimag(z1)
   52     continue
   50   continue

   30   continue
   20 continue
   10 continue

C      if (iprint() >= 50) then
C        call o1r(p00,lpdim,' p0r   ')
C        call o1r(p00(1,1,2),lpdim,' p0i   ')
C      endif

      call tcx('asxp0')
      end

      subroutine p0loc(iq,wght,p0,p0l,lpdim)
      implicit none
      integer iq,lpdim,i
      double precision wght(iq),p0(lpdim,lpdim,2),p0l(lpdim,lpdim,2)

      do  10  i = 1, lpdim
        p0l(i,i,1) = p0l(i,i,1) + abs(wght(iq))/2*p0(i,i,1)
        p0l(i,i,2) = p0l(i,i,2) + abs(wght(iq))/2*p0(i,i,2)
   10 continue
      end

      subroutine p0q0(q,lpdim,p0,p1,p00,p01,p02)
C- Make local part of q->0 limit of P0
      implicit none
      integer lpdim
      double precision p0(lpdim,lpdim,2),p1(lpdim,lpdim,2),q(3),
     .  p00(lpdim,lpdim),p01(lpdim,lpdim),p02(lpdim,lpdim),qq
      integer i,j

      qq = q(1)**2+q(2)**2+q(3)**2
      do  10  i = 1, lpdim
      do  10  j = 1, lpdim
        p00(i,j) = p0(i,j,1)
        p01(i,j) = p1(i,j,2)/dsqrt(qq)
        p02(i,j) = (p1(i,j,1)-p0(i,j,1))/qq
   10 continue

C      if (iprint() >= 50) then
C        call o1r(p00,lpdim,' p00  ')
C        call o1r(p01,lpdim,' p01  ')
C        call o1r(p02,lpdim,' p02  ')
C      endif
      end
