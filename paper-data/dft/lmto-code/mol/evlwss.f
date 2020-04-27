      subroutine evlwss(qval,qsmc,ipan,nhs,evl,del,ewt,nstate,nsp)
C - make occupation number weights
      implicit none
      integer ipan,nsp,nhs,nstate(nsp)
      double precision evl(nhs,1),ewt(nhs,1),del,qval,qsmc
      integer iprint,i1mach,nbpw,obmap,owk,nkp,i,isp
      real w(1)
      common /w/ w

      if (nsp == 1) then
        call evlwsm(nsp,qval,qsmc,ipan,nhs,evl,del,ewt,nstate)
      else
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
        nkp=1
        call defrr(obmap,nhs*nsp*nkp/nbpw+1)
        call dpzero(w(obmap),nhs*nsp*nkp/nbpw+1)
        call defdr(owk, nhs*nsp)
        call ebcpl(0,nhs,nhs,nsp,1,nkp,nbpw,w(obmap),w(owk),evl)
        call evlwsm(nsp,qval,qsmc,ipan,nhs,evl,del,ewt,nstate)
        call ebcpl(1,nhs,nhs,nsp,1,nkp,nbpw,w(obmap),w(owk),evl)
        call ebcpl(1,nhs,nhs,nsp,1,nkp,nbpw,w(obmap),w(owk),ewt)
        call rlse(obmap)
        if (nstate(1) == nhs) then
          nstate(2) = nhs
        else
          do  isp = 1, nsp
            nstate(isp)=0
            do  i = nhs, 1, -1
              if (ewt(i,isp) > 1d-8) then
                nstate(isp) = i
                goto 1
              endif
            enddo
    1       continue
          enddo
        endif
        if (iprint() >= 30) then
          call awrit2(' EVLWSS: nstate=%i,%i',' ',120,i1mach(2),
     .                nstate(1),nstate(2))
        endif
      endif
      end
      subroutine ebcpl(mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,wk,eb)
C- Gather spin-polarized bands into a single group, or redistribute
Ci Inputs
Ci  mode: 0, gather; 1, scatter
Ci  nbpw: a number no larger than the number of bits per integer word
Ci  bmap: an integer work array with dimension at nbmx*nsp*nq/nbpw
Ci  wk:   a work array with dimension nbmx*nsp
Ci  nbmx
      implicit none
      integer mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(1)
      double precision eb(nbmx,nsp,nq),wk(nbmx*nsp)
      integer ib,iq,ib1,ib2,iqb,iget

      if (nsp == 1 .or. nspc == 2) return

C --- Gather bands at each qp ---
      if (mode == 0) then
      iqb = 0
      do  10  iq = 1, nq

C   ... Gather and order +,- bands at this qp into one column
        ib1 = 1
        ib2 = 1
        do  20  ib = 1, nevx*nsp
          iqb = iqb+1
          if (eb(ib1,1,iq) < eb(ib2,2,iq) .and. ib1 <= nevx
     .      .or. ib2 > nevx) then
            wk(ib) = eb(ib1,1,iq)
            ib1 = ib1+1
          else
            wk(ib) = eb(ib2,2,iq)
            ib2 = ib2+1
            call mark1(bmap, nbpw, iqb)
          endif
   20   continue
        call dpcopy(wk,eb(1,1,iq),1,nevx*nsp,1d0)
        if (ib1-1 /= nevx .and. ib2-1 /= nevx) call rx('bug')
   10 continue
      endif

C --- Disperse bands at each qp ---
      if (mode == 1) then
      iqb = 0
      do  110  iq = 1, nq

C   ... Disperse bands into +,- for this qp according to bmap
        ib1 = 1
        ib2 = 1
        call dpcopy(eb(1,1,iq),wk,1,nevx*nsp,1d0)
        do  120  ib = 1, nevx*nsp
          iqb = iqb+1
          if (iget(bmap,nbpw,iqb) == 0) then
            eb(ib1,1,iq) = wk(ib)
            ib1 = ib1+1
          else
            eb(ib2,2,iq) = wk(ib)
            ib2 = ib2+1
          endif
  120   continue
        if (ib1-1 /= nevx .and. ib2-1 /= nevx) call rx('bug')
  110 continue
      endif

      end
      subroutine mark1(bitmap, nbpw, n)
C- put a one in the nth bit of bitmap.
C ----------------------------------------------------------------
Ci Inputs
Ci   bitmap, n
Cr Remarks
Cr    nbpw: a number no larger than the number of bits per integer word
C ----------------------------------------------------------------
      implicit none
      integer bitmap(1), nbpw, n
C Local parameters
      integer nword,nbit,i

      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      i = 2**(nbpw-nbit-1)
      bitmap(nword+1) = bitmap(nword+1) + i*(1-mod(bitmap(nword+1)/i,2))
      end
      integer function iget(bitmap, nbpw, n)
C- Return 0 or 1, depending on the value of the nth bit of bitmap
C ----------------------------------------------------------------
Cr Remarks
Cr   See mark1
C ----------------------------------------------------------------
      implicit none
      integer bitmap(1), nbpw, n
C Local parameters
      integer nword,nbit
      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
      end
