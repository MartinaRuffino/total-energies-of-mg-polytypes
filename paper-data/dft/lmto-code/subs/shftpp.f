      subroutine shftpp(nc,nlsp,pp,vold,vnew,oshft,nshft)
C- Shift or undo shift of pp's by constant potential
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nc    :number of classes
Ci   nlsp  :nl*nsp
Ci   pp    :potential parameters (atomsr.f)
Ci   vold  :see Remarks
Ci   vnew  :see Remarks
Ci   oshft :T, subtract vold from pp; see Remarks
Ci   nshft :T, add      vnew to   pp; see Remarks
Cr Remarks
Cr   Band center Cnu and linearization enu are shifted as follows:
Cr       oshft\nshft   F                 T
Cr             ---------------------------------------
Cr         F  |   do nothing       shift by vnew
Cr         T  | shift by -vold   shift by vnew-vold
Cr             ---------------------------------------
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical oshft,nshft
      integer nc,nlsp
      double precision pp(6,nlsp,nc),vold(nc),vnew(nc)
C ... Local parameters
      integer ic

      if (oshft) then
        call info2(60,0,0,' shiftpp sub old %n;11,6D',nc,vold)
      endif
      if (nshft) then
        call info2(60,0,0,' shiftpp add new %n;11,6D',nc,vnew)
      endif

C --- Shift enu and c for each class ---
      do  ic = 1, nc
        if (oshft) then
          call daxpy(nlsp,-1d0,vold(ic),0,pp(1,1,ic),6)
          call daxpy(nlsp,-1d0,vold(ic),0,pp(2,1,ic),6)
        endif

        if (nshft) then
          call daxpy(nlsp,1d0,vnew(ic),0,pp(1,1,ic),6)
          call daxpy(nlsp,1d0,vnew(ic),0,pp(2,1,ic),6)
        endif
      enddo
      end

      subroutine shftppr(lrel,nc,nl,nsp,pp,pprel,vold,vnew,oshft,nshft)
C- Shift or undo shift of relativistic pp's by constant potential
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lrel  :0 shift pp only
Ci         :1 shift pp only
Ci         :2 shift pp and pprel
Ci   nc    :number of classes
Ci   nlsp  :nl number of l's
Ci   pprel :relativistic potential parameters (atomsr.f)
Ci   vold  :see Remarks
Ci   vnew  :see Remarks
Ci   oshft :T, subtract vold from pp; see Remarks
Ci   nshft :T, add      vnew to   pp; see Remarks
Cr Remarks
Cr   Band center Cnu and linearization enu are shifted as follows:
Cr       oshft\nshft   F                 T
Cr             ---------------------------------------
Cr         F  |   do nothing       shift by vnew
Cr         T  | shift by -vold   shift by vnew-vold
Cr             ---------------------------------------
Cu Updates
Cu   21 Apr 15 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical oshft,nshft
      integer lrel,nc,nl,nsp
      double precision vold(nc),vnew(nc),pp(6,nl*nsp,nc),pprel(5,nl,2*nl,2,2)
C ... Local parameters
      integer ic,il,imu
      integer,parameter :: NULLI=-99999

      if (lrel >= 0 .and. lrel <= 2) call shftpp(nc,nl*nsp,pp,vold,vnew,oshft,nshft)
      if (lrel /= 2) return

C --- Shift enu and c for each class ---
      do  ic = 1, nc
        if (oshft) then
          do  il = 1, nl
            do  imu = 1, 2*il
              if (pprel(5,il,imu,1,1) /= NULLI) then
                pprel(1,il,imu,1,1) = pprel(1,il,imu,1,1) - vold(ic)
                pprel(1,il,imu,2,2) = pprel(1,il,imu,2,2) - vold(ic)
                pprel(5,il,imu,1,1) = pprel(5,il,imu,1,1) - vold(ic)
                pprel(5,il,imu,2,2) = pprel(5,il,imu,2,2) - vold(ic)
              endif
            enddo
          enddo
        endif

        if (nshft) then
          do  il = 1, nl
            do  imu = 1, 2*il
              if (pprel(5,il,imu,1,1) /= NULLI) then
                pprel(1,il,imu,1,1) = pprel(1,il,imu,1,1) + vnew(ic)
                pprel(1,il,imu,2,2) = pprel(1,il,imu,2,2) + vnew(ic)
                pprel(5,il,imu,1,1) = pprel(5,il,imu,1,1) + vnew(ic)
                pprel(5,il,imu,2,2) = pprel(5,il,imu,2,2) + vnew(ic)
              endif
            enddo
          enddo
        endif

      enddo
      end


      subroutine shftpph(nl,nsp,nbas,lihdim,ipc,indxsh,vold,vnew,oshft,nshft,pph)
C- Create a vector of site-dependent potential parameters
C ----------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   lihdim:size of lower + downfolding block
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   vold  :see Remarks
Ci   vnew  :see Remarks
Ci   oshft :T, subtract vold from pp; see Remarks
Ci   nshft :T, add      vnew to   pp; see Remarks
Co Input/Outputs
Co   pph   :potential parameters in downfolding order
Co         :band center shifted (see Remarks)
Cr Remarks
Cr   Band center Cnu and linearization enu are shifted as follows:
Cr       oshft\nshft   F                 T
Cr             ---------------------------------------
Cr         F  |   do nothing       shift by vnew
Cr         T  | shift by -vold   shift by vnew-vold
Cr             ---------------------------------------
Cr
Cr   Mapping of composite RL index to this vector is same as that of
Cr   eigenvectors permutated according to downfolding rules in indxsh
Cr   pph(1) : enu           (shifted)
Cr   pph(2) : calpha        (shifted)
Cr   pph(3) : sqrdel        (untouched)
Cr   pph(4) : palphaz       (untouched)
Cr   pph(5) : oalp          (untouched)
Cu Updates
Cu   23 Sep 13 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nsp,nbas,lihdim,ipc(nbas),indxsh(lihdim)
      logical oshft,nshft
      double precision pph(5,lihdim,nsp),vold(nbas),vnew(nbas)
C ... Local parameters
      integer ibas,ic,l,m,lmr,isp,k,n

C     call yprm('pph',1,pph,0,5,5,lihdim*nsp)

      do  isp = 1, nsp
      lmr = 0
        do  ibas = 1, nbas
        ic = ipc(ibas)
          do  l = 0, nl-1
          k = l+1
          if (indxsh(lmr+1) > lihdim) then
            lmr = lmr + 2*l+1
              cycle
          endif
          do  m = -l, l
            lmr = lmr+1
            n = indxsh(lmr)
            if (oshft) then
              pph(1,n,isp) = pph(1,n,isp) - vold(ibas)
              pph(2,n,isp) = pph(2,n,isp) - vold(ibas)
            endif

            if (nshft) then
              pph(1,n,isp) = pph(1,n,isp) + vnew(ibas)
              pph(2,n,isp) = pph(2,n,isp) + vnew(ibas)
            endif

          enddo
          enddo
        enddo
      enddo

C     call yprm('pph',1,pph,0,5,5,lihdim*nsp)
      end
