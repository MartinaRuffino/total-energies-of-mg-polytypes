      subroutine evlwsm(nsp,qval,qsmc,ipan,nhs,evl,del,ewt,nstate)
C- Set weights for eigenvalues.
C ----------------------------------------------------------------------
Ci Inputs:  nsp, qval, qsmc, ipan, nhs, evl, del
Ci
Co Outputs: ewt, nstate
Co
Cr Remarks: adaped from mc evlwsm. Each eigenvalue delta function is
Cr          broadened to a gaussian and the resulting density of states
Cr          is integrated. Returns eigenvalue weights for the charge
Cr          density accumulation.
Cr
Cr          In the spin polarised case the eigenvalues are united into
Cr          an ordered set with increasing energy. Their spin cannot be
Cr          identified here, but they are scattered back into two sets
Cr          on exit in evlwss.
Cr
Cr          The user may set -excite=#1,#2 in the command line in which
Cr          case evlwsm will unoccupy the HOMO-#1 and occupy the
Cr          LUMO+#2 eigenvectors. This creates an approximation to an
Cr          excited state. For example -excite=0,0 will create the
Cr          lowest "excited state" while -excite=0,1 takes an electron
Cr          from the HOMO and puts it in the second unoccupied level.
Cr          --LUMO is a synonym for -excite=0,0
Cr
Cr          ewt(*,2-3) is used as workspace
Cu Updates:
Cu         (ATP) generalised the "excited state" to -excite=#1,#2
Cu         (ATP) added --LUMO Feb 2005
Cu         11.01.95  spin polarized
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nsp,ipan,nhs
      double precision qval,qsmc,del
      double precision evl(1),ewt(nhs*nsp,3)
C Local Variables
      integer j1,j2,is,nst0,nstate,nste(2),i,it(2),nqval,nadd,ip
      integer i1mach,iprint,ipr,parg
      double precision wsum,state,dsum
      logical cmdopt, excite
      character*256 strn

      nste(1) = 0
      nste(2) = 0
      if (cmdopt('--LUMO',6,0,strn)) then
        if (nsp /= 2) call rx0(' --LUMO needs nsp=2')
        excite = .true.
      elseif (cmdopt('-excite=',8,0,strn)) then
        ip = 8
        call skipbl(strn,len(strn),ip)
        i = parg(' ',2,strn,ip,len(strn),', ',2,2,it,nste)
        if (i < 0) then
          call rxs2('EVLWSM: failed to parse "',strn(1:ip+5),' ..."')
        else
          if (nsp /= 2) call rx0(' --LUMO needs nsp=2')
          excite = .true.
        endif
      else
        excite = .false.
      endif
      if (excite .and. iprint() > 20) then
        call awrit2(' EVLWSM making excited state HOMO-%i --> LUMO+%i',
     .              ' ',128,i1mach(2),nste(1),nste(2))
      endif

C --- Take care of semicore panel ---
      if (ipan == 2) then
        nqval = qsmc + 0.01d0
        nstate = (qsmc + 0.01d0)/2
        if(nstate*2 /= nqval) call rx('EVLWSM: odd el num for 2nd pan')
        nstate = nsp*nstate
        do is = 1, nstate
          ewt(is,1) = 1d0
        enddo
        do  is = nstate + 1, nhs*nsp
          ewt(is,1) = 0d0
        enddo
        return
      endif

      nst0 = qval + del
      nadd = 0
      call molwts(nsp,qval,nadd,del,nhs,evl,ewt,nstate)

      if (excite) then
C --- make a hole at HOMO-#1 ---
        nadd = -(1 + nste(1))
        ipr = max(10,iprint() - 20)
        ipr = iprint()
        call pshprt(ipr)
        call molwts(nsp,qval,nadd,del,nhs,evl,ewt(1,2),nstate)
        call popprt
        wsum = 0d0
        do i = 1, nsp*nhs
          ewt(i,3) = ewt(i,1) - ewt(i,2)
          wsum = wsum - ewt(i,3)
        enddo
        call daxpy(nsp*nhs,-1d0,ewt(1,3),1,ewt,1)
        if (iprint() >= 20) then
          call awrit3(' EVLWSM: made weights for qval=%d-%d, sum=%d',
     .                ' ',128,i1mach(2),qval+nadd,qval,wsum)
        endif
C --- occupy LMUO+#2 ---
        nadd = 1 + nste(2)
        call pshprt(ipr)
        call molwts(nsp,qval,nadd,del,nhs,evl,ewt(1,2),nstate)
        call molwts(nsp,qval,nadd-1,del,nhs,evl,ewt(1,3),nstate)
        call popprt
        wsum = 0d0
        do i = 1, nsp*nhs
          ewt(i,2) = ewt(i,2) - ewt(i,3)
          wsum = wsum + ewt(i,2)
        enddo
        call daxpy(nsp*nhs,1d0,ewt(1,2),1,ewt,1)
        if (iprint() >= 20) then
          call awrit3(' EVLWSM: made weights for qval=%d-%d, sum=%d',
     .                ' ',128,i1mach(2),qval+nadd,qval+nadd-1,wsum)
        endif

C --- for excited states, set nstate = nhs for making density matrix ---
        state = dsum(nsp*nhs,ewt,1)

C --- Printout ---
        if(iprint() >= 20) then
          call awrit2(' EVLWSM: del=%d state=%d ',' ',120,i1mach(2),
     .                del,state)
          if(iprint() >= 30) then
            j1 = max0(nst0 - 6,1)
            j2 = min0(nst0 + 6,nhs*nsp)
            write(6,10) 'evl:',(evl(i),i = j1,j2)
            write(6,10) 'ewt:',(ewt(i,1),i = j1,j2)
          endif
        endif
   10   format(1x,a4,13f8.4)
      endif

      end

      subroutine molwts(nsp,qval,nadd,del,nhs,evl,ewt,nstate)
C- make weights for evl, qval+nadd electrons
      implicit none
      integer nhs,nsp,nadd,nstate
      double precision qval,del,evl(1),ewt(1)

      integer nqval,nst0,i1,i2,i,irep,j1,j2
      integer i1mach,iprint
      double precision pi,datan,srpi,wsum,ef,sum,dum,dexp,derf,x,
     .                 e1,e2,e0,state,dsum

      pi = 4d0*datan(1d0)
      srpi = dsqrt(pi)

C     nqval = qval + del
      wsum = nsp*qval/2d0
      nst0 = wsum + del

      wsum = wsum + nadd
      i1 = 1
      i2 = nhs*nsp
      e1 = evl(i1)
      e2 = evl(i2)
      ef = evl(nst0+nadd)

C --- Repeat until Fermi energy found ---
      if (iprint() >= 45) print *, 'MOLWTS: seek ''Fermi energy''..'
      do  irep = 1, 30
        if (ef < e1 .or. ef > e2) ef = 0.5d0*(e1 + e2)
        sum = 0d0
        dum = 0d0
        do  i = 1, nhs*nsp
          x = (evl(i) - ef)/del
          ewt(i) = 0.5d0*(1d0 - derf(x))
          sum = sum + ewt(i)
          dum = dum + dexp(-x*x)/(del*srpi)
        enddo
        if (iprint() >= 45) then
          call awrit6(' trial %#2i %#8,8D %#8,8D %#8,8D '//
     .                'sum,wsum=%,8d,%,8d',
     .                ' ',120,i1mach(2),irep,e1,ef,e2,sum,wsum)
        endif
        if (sum > wsum) e2 = ef
        if (sum <= wsum) e1 = ef
        if (dabs(wsum - sum) < 1d-8*nsp) goto 1
        e0 = ef
        ef = 0.5d0*(e1 + e2)
        if (dum > 1d-8) then
          if (iprint() >= 45) then
            call awrit1('  weight at ef: %d',' ',120,i1mach(2),dum)
          endif
          ef = e0 - (sum - wsum)/dum
        endif
      enddo

      call rx('MOLWTS: can''t find fermi energy')
    1 continue

C --- Get "number of states" ---
      state = dsum(nsp*nhs,ewt,1)
      nstate = int(state + del)

C --- Printout ---
      if(iprint() >= 20) then
        call awrit6(' MOLWTS: del=%d nst0=%i nstate=%i'
     .              //' state=%d ef=%d evl=%d',' ',120,i1mach(2),
     .              del,nst0,nstate,state,ef,evl(nstate))
        if(iprint() >= 30) then
          j1 = max0(nst0 - 4,1)
          j2 = min0(nst0 + 4,nhs*nsp)
          write(6,10) 'evl:',(evl(i),i = j1,j2)
          write(6,10) 'ewt:',(ewt(i),i = j1,j2)
        endif
      endif
   10 format(1x,a4,9f8.4)

      end
