      subroutine vbare(mode,nsp,aintra,vintra,vmad,vmad1,nbas,nl,
     .  lpdim,indxsh,ipc,ldim,wsr,v,v1)
C- Bare coulomb potential matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 for normal q; mode=1 for limit q-> 0 (see madmtq.f)
Ci   aintra:T if v is l-dependent
Ci   vintra:intra-atomic contribution
Ci   vmad  :Madelung matrix without on-site contributions
Ci   vmad1
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   lpdim :dimensions v,v1
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   wsr   :Wigner-Seitz radius, in a.u. (input; alias rmax)
Co Outputs
Co   v     :Case no vintra: indices i,j to v(i,j) are site indices R
Co         :  v_ij = v_R,R' = 2*vmad_R,R' + 2/wsr delta_R,R'
Co   v     :Case vintra: indices i,j to v(i,j) are compound indices Rl
Co         :  v_ij = 2*vmad_R,R' + (2/wsr + vintra(l,l')) delta_R,R'
Co   v1    :
Cr Remarks
Cu Updates
Cu   06 Oct 04 (T. Sandu) Adapted to handle spin polarized case;
Cu                        vintra independent of spin
Cu   03 Oct 01 Adapted to handle downfolding
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical aintra
      integer nsp,isp1,isp2
      integer mode,nbas,nl,ldim,lpdim,ipc(nbas),indxsh(*)
      double precision vintra(nl,nl,nsp,nsp,nbas),wsr(*),
     .  v(lpdim,lpdim,2),v1(lpdim,lpdim),vmad(nbas,nbas,2),
     .  vmad1(nbas,nbas)
C ... Local parameters
      integer iy,ib,il,jy,jb,jl,lmi,iy0,mxorb,nglob,offi,lmj,offj

      mxorb = nglob('mxorb')

C --- Case v is l-dependent ----
      if (aintra) then
        iy = 0
        do  10  ib = 1, nbas
C       Hang onto offset to start of block corresponding to ib
        iy0 = iy
        lmi = mxorb*(ib-1)
        do  10  il = 1, nl
          offi = indxsh(lmi+1) - 1
          lmi = lmi + 2*il-1
          if (offi < ldim) then
            iy = iy+1

C           Sanity check
            if (iy > lpdim) call rx('bad indexing in vbare')

C        ...Interatomic contribution
            jy = 0
            do  14  jb = 1, nbas
            lmj = mxorb*(jb-1)
            if (mode == 0) then
              do  15  jl = 1, nl
              offj = indxsh(lmj+1) - 1
              lmj = lmj + 2*jl-1
              if (offj < ldim) then
                jy = jy+1
                v(iy,jy,1) = 2d0*vmad(ib,jb,1)
                v(iy,jy,2) = 2d0*vmad(ib,jb,2)
              endif
   15         continue
            else
              do  16  jl = 1, nl
              offj = indxsh(lmj+1) - 1
              lmj = lmj + 2*jl-1
              if (offj < ldim) then
                jy = jy+1
                v(iy,jy,1) = 2d0*vmad(ib,jb,1)
                v1(iy,jy)  = 2d0*vmad1(ib,jb)
              endif
   16         continue
            endif
   14       continue

C        ...Intraatomic contribution: jb=ib
            jy = iy0
            lmj = mxorb*(ib-1)
            do  18  jl = 1, nl
              offj = indxsh(lmj+1) - 1
              lmj = lmj + 2*jl-1
              if (offj < ldim) then
                jy = jy+1
                if (nsp == 1) then
                  v(iy,jy,1) = v(iy,jy,1) +
     .              vintra(il,jl,nsp,nsp,ipc(ib)) + 2d0/wsr(ipc(ib))
                 elseif (nsp == 2) then
                  do 101 isp1 = 1, nsp
                  do 101 isp2 = 1, nsp
                   v(iy,jy,1) = v(iy,jy,1) +
     .               0.25d0*vintra(il,jl,isp1,isp2,ipc(ib))
 101              continue
                  v(iy,jy,1) = v(iy,jy,1) + 2d0/wsr(ipc(ib))
                 else
                   call rx('bad spin index in vbare')
                 endif
              endif
   18       continue
          endif
   10   continue

C --- Case v is l-independent ----
      else
        do  20  ib = 1, nbas
          if (mode == 0) then
            do  22  jb = 1, nbas
            v(ib,jb,1) = 2d0*vmad(ib,jb,1)
   22       v(ib,jb,2) = 2d0*vmad(ib,jb,2)
          else
            do  24  jb = 1, nbas
            v(ib,jb,1) = 2d0*vmad(ib,jb,1)
   24       v1(ib,jb)  = 2d0*vmad1(ib,jb)
          endif
          v(ib,ib,1) = v(ib,ib,1) + 2d0/wsr(ipc(ib))
   20   continue
      endif

C     call yprm('vbare',2,v,lpdim*lpdim,lpdim,lpdim,lpdim)

      end
