      subroutine cffor(nbas,nl,nsp,nsite,ldim,ib,lmx,ipc,iax,npr,
     .  indxsh,s,z,indxH,pv,f)
C- Assembly crystal field force terms
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas,nl,nsp,lmx,ipc;
Ci   nsite: total number of neighbors in all clusters;
Ci   ldim: dimension of basis set;
Ci   ib: band index
Ci   iax: neighbor lists;
Co   npr: pointers for neighbor lists
Ci   indxsh: pointers for constructing permutated H and S from makidx
Ci   s: real-space matrix that is to be summed;
Ci   z: eigenvectors for band ib
Ci   indxH: work array of length nbas
Ci   pv: true if calculating 3PV
Co Outputs
Co   f(ibas): Contribution to force for atom ibas and band ib
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nbas,nl,nsp,nsite,ldim,ib,niax
      parameter (niax=10)
      integer lmx(*),ipc(nbas),iax(0:niax-1,nsite),npr(0:1,nbas),
     .  indxsh(*),indxH(nbas)
      double precision s(nl**2,nl**2,nsite*nsp**2),z(ldim,ldim,2),
     .  f(nbas)
      logical pv

C Local parameters
      integer i,lsp,isite,iat0,iatR,ix0,nlm0,lm0,lm1,j
      integer nlm
      double precision z2s

C Intrinsic functions
      nlm(i) = (1+lmx(i))**2

C --- Form table of indices that mark offset of Rth block ---
      indxH(1) = 0
      do  i = 2, nbas
        indxH(i) = indxH(i-1) + nlm(ipc(i-1))
      enddo

C --- For each RR' pair add contribution to force sums ---
      lsp = ldim / nsp
      do  isite = 1, nsite
        iat0 = iax(0,isite)
        iatR = iax(1,isite)
C ... Skip site with itself
        if (isite == npr(1,iat0)+1) cycle
        ix0 = indxH(iat0)
        nlm0 = nlm(ipc(iat0))
        do  lm0 = 1, nlm0
          i = indxsh(ix0+lm0)
C ... Throw out higher waves here
          if (i > lsp) cycle
          do  lm1 = 1, nlm0
            j = indxsh(ix0+lm1)
C ... Throw out higher waves here
            if (j > lsp) cycle
            z2s = (z(i,ib,1)*z(j,ib,1)
     .        + z(i,ib,2)*z(j,ib,2))*s(lm1,lm0,isite)
            if (nsp == 2) then
              z2s = z2s + (z(i,ib,1)*z(j+lsp,ib,1)
     .          + z(i,ib,2)*z(j+lsp,ib,2))*s(lm1,lm0,isite+2*nsite)
              z2s = z2s + (z(i+lsp,ib,1)*z(j,ib,1)
     .          + z(i+lsp,ib,2)*z(j,ib,2))*s(lm1,lm0,isite+1*nsite)
              z2s = z2s + (z(i+lsp,ib,1)*z(j+lsp,ib,1)
     .          + z(i+lsp,ib,2)*z(j+lsp,ib,2))*s(lm1,lm0,isite+3*nsite)
            endif
            f(iat0) = f(iat0) + z2s
            if (.not. pv) f(iatR) = f(iatR) - z2s
          enddo
        enddo
      enddo

      end
