      subroutine asldau(nbas,ldim,lmaxu,isp,nsp,nspc,nlibu,idu,
     .  indxsh,hp1,vorb,hk)
C- Add (1+oh)+ V(LDA+U) (1+oh) to the LDA hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0  Copy LDA+U potential to hk
Ci         :1  Add  LDA+U potential to hk
Ci         :2  Copy diagonal part of LDA+U potential to hd in
Ci             hamiltonian order
Ci         :3  Add diagonal part of LDA+U potential to hd in
Ci             hamiltonian order
Ci         :4  Copy LDA+U potential to hk, sans diagonal part
Ci         :5  Add LDA+U potential to hk, sans diagonal part
Ci   nbas  :size of basis
Ci   ldim  :dimension of hamiltonian matrix (makidx.f):
Ci         :orbitals i for which indxsh(i)>ldim are not included in h
Ci   lmaxu :dimensioning parameter for U matrix
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   isp   :current spin channel (1 or 2)
Ci   nlibu :number of U blocks (LDA+U)
Ci   idu   :idu(l+1,ib)>0 => this l has a nonlocal U matrix
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   vorb  :orbital dependent potential matrices
Co   hp1   :1+oh
Co Input/Outputs
Co   hk:   :LDA+U contribution added to hk
Cl Local variables
Cl   nl    :(global maximum l) + 1
Cu Updates
Cu   15 Nov 07 Replaced potential with (1+oh)+ V[U] (1+oh)
Cu   08 Nov 07 (J. Xu) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ldim,lmaxu,isp,nsp,nspc,nlibu
      integer idu(4,nbas)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double precision hk(ldim,nspc,ldim,nspc),hp1(ldim,nspc,ldim,nspc)
      integer indxsh(*)
C ... Dynamically allocated arrays
      complex(8), allocatable :: vh(:),wk(:)
C ... Local parameters
      integer ldimx
      double precision xx

      ldimx = ldim*nspc
      allocate(vh(ldimx**2)); call dpzero(vh,2*ldimx**2)
      allocate(wk(ldimx**2)); call dpzero(wk,2*ldimx**2)
      call asaddu(4,nbas,ldim,lmaxu,isp,nsp,nspc,nlibu,idu,
     .  indxsh,vorb,vh,xx)
C     call yprm('V(U)',2,vh,ldimx*ldimx,ldimx,ldimx,ldimx)
C     call yprm('1+oh',2,hp1,ldimx*ldimx,ldimx,ldimx,ldimx)
      call ygemm('N','N',ldimx,ldimx,ldimx,1d0,
     .  vh,ldimx**2,ldimx,hp1,ldimx**2,ldimx,0d0,
     .  wk,ldimx**2,ldimx)
C     call yprm('V(U)(1+oh)',2,wk,ldimx*ldimx,ldimx,ldimx,ldimx)
      call ygemm('C','N',ldimx,ldimx,ldimx,1d0,
     .  hp1,ldimx**2,ldimx,wk,ldimx**2,ldimx,0d0,
     .  vh,ldimx**2,ldimx)
C     call yprm('(1+oh)+V(U)(1+oh)',2,vh,
C    .  ldimx*ldimx,ldimx,ldimx,ldimx)
      call daxpy(ldimx*ldimx*2,1d0,vh,1,hk,1)
      deallocate(vh,wk)

      end
      subroutine asaddu(mode,nbas,ldim,lmaxu,isp,nsp,nspc,nlibu,idu,
     .  indxsh,vorb,hk,hd)
C- Add LDA+U hamiltonian to the full hamiltonian matrix hk
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0  Copy LDA+U potential to hk
Ci         :1  Add  LDA+U potential to hk
Ci         :2  Copy diagonal part of LDA+U potential to hd in
Ci             hamiltonian order
Ci         :3  Add diagonal part of LDA+U potential to hd in
Ci             hamiltonian order
Ci         :4  Copy LDA+U potential to hk, sans diagonal part
Ci         :5  Add LDA+U potential to hk, sans diagonal part
Ci   nbas  :size of basis
Ci   ldim  :dimension of hamiltonian matrix (makidx.f):
Ci         :orbitals i for which indxsh(i)>ldim are not included in h
Ci   lmaxu :dimensioning parameter for U matrix
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   isp   :current spin channel (1 or 2)
Ci   nlibu :number of U blocks (LDA+U)
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   vorb  :orbital dependent potential matrices
Co Input/Outputs
Co   hk:   :LDA+U contribution added to hk
Cl Local variables
Cl   nl    :(global maximum l) + 1
Cu Updates
Cu   08 Nov 07 (J. Xu) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,ldim,lmaxu,isp,nsp,nspc,nlibu
      integer idu(4,nbas),nglob
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double precision hk(ldim,nspc,ldim,nspc),hd(ldim)
      integer indxsh(*),ldimx
C ... Local parameters
      integer ib,il,lmr,lcount,ih,ih1,ih2,iv1,iv2
      integer mxorb,k1,k2,ll,nl,is,ks,ispc,iblu
      logical lfull

      ldimx = ldim*nspc
      if (mode == 0 .or. mode == 4) call dpzero(hk,ldimx*ldimx*2)
      if (mode == 2) call dpzero(hd,ldim)
      call sanrg(.true.,mode,0,6,'asaddu','mode')
      lfull = mode == 0 .or. mode == 1 .or. mode == 4 .or. mode == 5

C      print *, isp
C      if (lfull) then
C      call yprm('before H',2,hk,ldimx*ldimx,ldimx,ldimx,ldimx)
C      endif

C     is = spin index; ks = index to poke into hk
      is = isp
      ks = 1
      do  ispc = 1, nspc
      if (nspc == 2) then
        is = ispc
        ks = ispc
      endif
      mxorb = nglob('mxorb')
      nl = ll(mxorb)+1
      lmr = 1
      iblu = 0
      do  ib = 1, nbas
          lcount = 0
          do  il = 1, nl
            if (idu(il,ib) /= 0) then
            iblu = iblu+1
              ih = indxsh(lmr)
C         ... Add lda+u for site ib, l block il
              if (ih <= ldim .and. lfull) then
                do  k1 = 1, 2*il-1
                  do  k2 = 1, 2*il-1
                    ih1 = ih + k1 - 1
                    ih2 = ih + k2 - 1
                    iv1 = -il + k1
                    iv2 = -il + k2
                    hk(ih1,ks,ih2,ks) = hk(ih1,ks,ih2,ks) +
     .                dble(vorb(iv1,iv2,is,iblu))
                    hk(ih1,ks,ih2+ldimx,ks) = hk(ih1,ks,ih2+ldimx,ks) +
     .                dimag(vorb(iv1,iv2,is,iblu))
                  enddo
                  if (mode == 4 .or. mode == 5) then
                    hk(ih1,ks,ih1,ks) = hk(ih1,ks,ih1,ks) -
     .                dble(vorb(iv1,iv1,is,iblu))
                    hk(ih1,ks,ih1+ldimx,ks) = hk(ih1,ks,ih1+ldimx,ks) -
     .                dimag(vorb(iv1,iv1,is,iblu))
                  endif
                enddo
              elseif (ih <= ldim) then
                do  k1 = 1, 2*il-1
                  ih1 = ih + k1 - 1
                  iv1 = -il + k1
                  hd(ih1) = dble(vorb(iv1,iv1,is,iblu))
                enddo
              endif
              lcount = lcount + 1
            endif
            lmr = lmr + 2*il - 1
          enddo
      enddo
      enddo

C      if (lfull) then
C        print *, isp
C        call yprm('H, LDA+U',2,hk,ldimx*ldimx,ldimx,ldimx,ldimx)
C      else
C        call yprm('diagonal part of U',1,hd,0,ldim,ldim,1)
C      endif

      end

      subroutine u2pph(mode,nbas,lmaxu,nsp,nlibu,idu,vorb,indxsh,
     .  ldim,lihdim,fac,pph)
C- Add diagonal part of LDA+U potential into pph(1) and pph(2)
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 add diagonal part of LDA+U potential into pph
Ci         :1 return diagonal part of LDA+U potential in pph
Ci         :  In this case, pph is dimensioned pph(ldim,2)
Ci   nbas  :size of basis
Ci   lmaxu :dimensioning parameter for U matrix
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu :number of U blocks (LDA+U)
Ci
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   fac   :add fac*ppv into pph(1) and pph(2)
Co Outputs
Co   pph   :potential parameters in downfolding order (makpph.f)
C          :fac*(diagaonal LDA+U) potential added into pph(1),pph(2)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,lmaxu,nsp,nlibu,indxsh(*),ldim,lihdim
      integer idu(4,nbas)
      double precision fac,pph(5,lihdim,2)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
C ... Local parameters
      double precision ppv(ldim,2)

C     return
      if (nlibu == 0) return

C     call prmx('pph before adding V(LDA+U)',pph,5,5,2*lihdim)
      call asaddu(2,nbas,ldim,lmaxu,1,nsp,1,nlibu,idu,indxsh,vorb,
     .  ppv,ppv)
      call asaddu(2,nbas,ldim,lmaxu,nsp,nsp,1,nlibu,idu,indxsh,
     .  vorb,ppv,ppv(1,nsp))
C     call prmx('V(LDA+U)',ppv,ldim,ldim,2)

      if (mode == 1) then
        call dcopy(ldim*2,ppv,1,pph,1)
C       call prmx('V(LDA+U)',pph,ldim,ldim,2)
        return
      endif

      call daxpy(ldim,fac,ppv,1,pph(1,1,1),5)
      call daxpy(ldim,fac,ppv,1,pph(2,1,1),5)
      if (nsp == 2) then
        call daxpy(ldim,fac,ppv(1,2),1,pph(1,1,2),5)
        call daxpy(ldim,fac,ppv(1,2),1,pph(2,1,2),5)
      endif

C     call prmx('pph after adding V(LDA+U)',pph,5,5,2*lihdim)
      end
