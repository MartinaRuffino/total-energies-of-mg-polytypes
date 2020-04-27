      subroutine pmblk(nblk,ibp,ioffo,ioffn,nlao,asrc,mode,
     .  alfa,nlma,nlan,adst,nmax)
C- Permutes subblocks of an array and copiesinto a new array
Ci ioffo: sequence markers to start of subblocks to asrc.  End of
Ci        subblock implied by start of next subblock.  ioffo should be
Ci        always increasing, with dimension 1+no. of subblocks, to mark
Ci        end of last subblock.
Ci ibp:   permutation of subblocks to take from asrc or adst
Ci ioffn: like ioffo, but for adst.
Ci mode:  one's digit:  1, make ioffo offsets relative to ioffo(1)
Ci        ten's digit:  1, make ioffn offsets relative to ioffn(1)
Ci        100's digit:  0, no permutation of indices
Ci                      1, permute sequence of subblocks in asrc by ibp
Ci                      2, permute sequence of subblocks in adst by ibp
Ci        1e3's digit:  1, add into new array, not copy
Ci        1e4's digit:  1, scale adst(*,j) by alfa(j)
Ci                    >=8, subblocks in both rows and cols (see remarks)
Ci nlao,nlan: leading dimensions of asrc and adst
Ci nlma:  number of rows for which to permute columns.
Co adst:  result matrix
Co nmax:  largest column to adst addressed
Cr By default pmblk permutes columns around, doing the same for
Cr each of nlma rows.  if mode 1000's digit >=8, pmblk permutes
Cr both columns and rows.
Cr If size of destination subblock smaller than source, truncate.
Cr If instead larger, accumulate only smaller size.
C ----------------------------------------------------------------
      implicit none
      integer nblk,nlao,nlan,nlma,mode,ioffn(nblk),ioffo(nblk),
     .  ibp(nblk),nmax
      double precision adst(nlan,1),asrc(nlao,1),alfa(1)
      integer ilma,ib,iofo,ilm,iofn,nlm1,ibn,ibo,offo0,offn0,
     . nn,moded(5),jb,jofo,jlm,jofn,nlm2,jbn,jbo

      moded(1) = mod(mode,10)
      moded(2) = mod(mode/10,10)
      moded(3) = mod(mode/100,10)
      moded(4) = mod(mode/1000,10)
      moded(5) = mod(mode/10000,10)

C     if (moded(3) /= 0) call awrit2('pm ibp  %n,3i',' ',80,6,nblk,ibp)
C     call awrit2('pm iofo %n,3i',' ',80,6,nblk+1,ioffo)
C     call awrit2('pm iofn %n,3i',' ',80,6,nblk+1,ioffn)

      offn0 = 0
      offo0 = 0
      if (moded(1) == 1) offo0 = ioffo(1)
      if (moded(2) == 1) offn0 = ioffn(1)
      nmax = -1
      if (moded(5) < 8) then
        do  10  ib = 1, nblk
          ibo = ib
          if (moded(3) == 1) ibo = ibp(ib)
          ibn = ib
          if (moded(3) >= 2) ibn = ibp(ib)
          nn = ioffn(ibn+1)-ioffn(ibn)
          nlm1 = min(nn,ioffo(ibo+1)-ioffo(ibo))
          iofn = ioffn(ibn) - offn0
          iofo = ioffo(ibo) - offo0
          if (moded(4)+moded(5) == 0) then
            do  20  ilm = 1, nlm1
            do  20  ilma = 1, nlma
   20       adst(ilma,ilm+iofn) = asrc(ilma,ilm+iofo)
          elseif (moded(5) == 0) then
            do  25  ilm = 1, nlm1
            do  25  ilma = 1, nlma
   25       adst(ilma,ilm+iofn) = adst(ilma,ilm+iofn)
     .                          + asrc(ilma,ilm+iofo)
          elseif (moded(4) == 0) then
            do  30  ilm = 1, nlm1
            do  30  ilma = 1, nlma
   30       adst(ilma,ilm+iofn) = asrc(ilma,ilm+iofo)*alfa(ilm+iofn)
          else
            do  35  ilm = 1, nlm1
            do  35  ilma = 1, nlma
   35       adst(ilma,ilm+iofn) = adst(ilma,ilm+iofn)
     .                          + asrc(ilma,ilm+iofo)*alfa(ilm+iofn)
          endif
          nmax = max(nn+iofn,nmax)
   10   continue
      else
        do  100  ib = 1, nblk
          ibo = ib
          if (moded(3) == 1) ibo = ibp(ib)
          ibn = ib
          if (moded(3) >= 2) ibn = ibp(ib)
          nn = ioffn(ibn+1)-ioffn(ibn)
          nlm1 = min(nn,ioffo(ibo+1)-ioffo(ibo))
          iofn = ioffn(ibn) - offn0
          iofo = ioffo(ibo) - offo0
          do  110  jb = 1, nblk
            jbo = jb
            if (moded(3) == 1) jbo = ibp(jb)
            jbn = jb
            if (moded(3) >= 2) jbn = ibp(jb)
            nlm2 = min(ioffn(jbn+1)-ioffn(jbn),ioffo(jbo+1)-ioffo(jbo))
            jofn = ioffn(jbn) - offn0
            jofo = ioffo(jbo) - offo0
            if (moded(4) == 0) then
              do  120  jlm = 1, nlm2
              do  120  ilm = 1, nlm1
  120         adst(ilm+iofn,jlm+jofn) = asrc(ilm+iofo,jlm+jofo)
            else
              do  125  jlm = 1, nlm2
              do  125  ilm = 1, nlm1
  125         adst(ilm+iofn,jlm+jofn) = adst(ilm+iofn,jlm+jofn) +
     .                                  asrc(ilm+iofo,jlm+jofo)
            endif
  110     continue
          nmax = max(nn+iofn,nmax)
  100   continue
      endif
      end
