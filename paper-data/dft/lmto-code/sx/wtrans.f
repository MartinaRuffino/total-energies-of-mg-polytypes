      subroutine wtrans(ldim,lpdim,nl,nbas,bas,indxsh,istab,g,ag,q,
     .  aintra,h0,h)
C- Transform the potential h0(q) to h0(qstar)
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lpdim :total number of l's in lower block (or sites if aintra=F)
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   q     :q-point for h0
Ci   aintra:T: h is resolved into l and site; F: h is site-resolved only
Ci   h0    :the potential at q
Co Outputs
Co  h      :the potential transformed from h0.
Cl Local variables
Cr Remarks
Cr   V_Rl,R'l'(q) is transformed into V_Rl,R'l'(g q)
Cu Updates
Cu   03 Oct 01 Revamped to allow downfolding (assumes lower block only)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldim,nl,nbas,lpdim,indxsh(nl**2*nbas),istab(nbas)
      double precision bas(3,nbas),g(3,3),ag(3),tbas(3),q(3),
     .  h(lpdim,lpdim,2),h0(lpdim,lpdim,2)
      logical aintra
C ... Local parameters
      double precision twopi,sp,cosp1,sinp1,xxr,xxi,cosp,sinp
      integer i,i1,i2,ibas1,ibas2,il1,il2,iy1,iy2,j,j1,j2,jbas1,jbas2,
     .  jy1,jy2,l,l1,l2,lmr,mxorb,nglob,offi

      twopi = 8d0*datan(1d0)
      mxorb = nglob('mxorb')

C --- With intraatomic screening -----
      if (aintra) then

C --- For each pair RR' = ibas1,ibas2 do ---
C     iy1 = offset to start of h block, first index
C     jy1 = offset to start of h0 block, first index
      iy1 = 0
      do  10  ibas1 = 1, nbas
        lmr = mxorb*(ibas1-1)
        jbas1 = istab(ibas1)
C       Count l1=number of l's for ibas1,jbas1
        l1 = 0
        do  12  l = 1, nl
          offi = indxsh(lmr+1) - 1
          lmr = lmr + 2*l-1
          if (offi < ldim) then
            l1 = l1+1
          endif
   12   continue
C       Get jy1=offset to start of h0 block
        jy1 = 0
        do  14  j = 1, jbas1-1
        lmr = mxorb*(j-1)
        do  14  l = 1, nl
          offi = indxsh(lmr+1) - 1
          lmr = lmr + 2*l-1
          if (offi < ldim) then
            jy1 = jy1+1
          endif
   14   continue
C       Sanity check
        if (iy1+l1 > lpdim .or. jy1+l1 > lpdim)
     .    call rx('bad indexing in wtrans')

C   ... Phase between ibas1 and jbas1
        do  16  i = 1, 3
          tbas(i) = ag(i) - bas(i,ibas1)
          do  18  j = 1, 3
   18     tbas(i) = tbas(i) + g(i,j)*bas(j,jbas1)
   16   continue
        sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
*       phase1 = zexp(dcmplx(0d0,-sp))
        cosP1 =  dcos(sp)
        sinP1 = -dsin(sp)

C       iy2 = offset to start of h block, second index
C       jy2 = offset to start of h0 block, second index
        iy2 = 0
        do  20  ibas2 = 1, nbas
          jbas2 = istab(ibas2)

C         Count l2=number of l's for ibas2,jbas2
          l2 = 0
          lmr = mxorb*(ibas2-1)
          do  22  l = 1, nl
            offi = indxsh(lmr+1) - 1
            lmr = lmr + 2*l-1
            if (offi < ldim) then
              l2 = l2+1
            endif
   22     continue
C         Get jy2=offset to start of h0 block
          jy2 = 0
          do  24  j = 1, jbas2-1
          lmr = mxorb*(j-1)
          do  24  l = 1, nl
            offi = indxsh(lmr+1) - 1
            lmr = lmr + 2*l-1
            if (offi < ldim) then
              jy2 = jy2+1
            endif
   24     continue
C         Sanity check
          if (iy2+l2 > lpdim .or. jy2+l2 > lpdim)
     .      call rx('bad indexing in wtrans')

C     ... Phase between ibas2 and jbas2
          do  26  i = 1, 3
            tbas(i) = ag(i) - bas(i,ibas2)
            do  28  j = 1, 3
   28       tbas(i) = tbas(i) + g(i,j)*bas(j,jbas2)
   26     continue
          sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
*         phase = phase1*zexp(dcmplx(0d0,sp))
          xxr = dcos(sp)
          xxi = dsin(sp)
          cosP = cosP1*xxr - sinP1*xxi
          sinP = sinP1*xxr + cosP1*xxi

C     --- Loop over l ---
          do  30  il1 = 1, l1
            i1 = iy1+il1
            j1 = jy1+il1
            do  32  il2 = 1, l2
              i2 = iy2+il2
              j2 = jy2+il2
*             h(i1,i2) = phase*h0(j1,j2)
              h(i1,i2,1) = h0(j1,j2,1)*cosP - h0(j1,j2,2)*sinP
              h(i1,i2,2) = h0(j1,j2,2)*cosP + h0(j1,j2,1)*sinP
   32       continue
   30     continue
          iy2 = iy2+l2
   20   continue
        iy1 = iy1+l1
   10 continue

C --- Without intraatomic screening -----
      else
        do  40  ibas1 = 1, nbas
        jbas1 = istab(ibas1)
        do  42  i = 1, 3
          tbas(i) = ag(i) - bas(i,ibas1)
          do  44  j = 1, 3
   44     tbas(i) = tbas(i) + g(i,j)*bas(j,jbas1)
   42   continue
        sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
*       phase1 = zexp(dcmplx(0d0,-sp))
        cosP1 =  dcos(sp)
        sinP1 = -dsin(sp)
        do  46  ibas2 = 1, nbas
          jbas2 = istab(ibas2)
          do  48  i = 1, 3
            tbas(i) = ag(i) - bas(i,ibas2)
            do  50  j = 1, 3
   50       tbas(i) = tbas(i) + g(i,j)*bas(j,jbas2)
   48     continue
          sp = twopi*(tbas(1)*q(1)+tbas(2)*q(2) + tbas(3)*q(3))
*         phase = phase1*zexp(dcmplx(0d0,sp))
          xxr = dcos(sp)
          xxi = dsin(sp)
          cosP = cosP1*xxr - sinP1*xxi
          sinP = sinP1*xxr + cosP1*xxi
*         h(ibas1,ibas2) = phase*h0(jbas1,jbas2)
          h(ibas1,ibas2,1) = cosP*h0(jbas1,jbas2,1)
     .                     - sinP*h0(jbas1,jbas2,2)
          h(ibas1,ibas2,2) = cosP*h0(jbas1,jbas2,2)
     .                     + sinP*h0(jbas1,jbas2,1)

   46   continue
   40 continue
      endif
      end
