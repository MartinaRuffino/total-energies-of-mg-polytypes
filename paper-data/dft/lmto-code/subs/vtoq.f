      subroutine vtoq(nbas,nl,nsp,nclass,ipc,nrclas,emad,frzves,
     .  lmx,clabl,vrmax,rhrmx,rmax,ves,mad,z,pnu,wk,dq,qnu)
C- Adjusts moments to conform to Madelung potential or energy
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas,nl,nsp,nclass,ipc,nrclas,rmax,dq,ves,mad,z,pnu
Ci   emad when frzves = 1
Ci   ves, when frzves = 2
Co Outputs
Co   wk:  inverse M of Madelung matrix
Co   dq:  total charges in each sphere
Co   qnu: adjusted moments
Co   ves  when frzves = 1
Cr Remarks
Cr   Case frzves=1:
Cr     Idea is to readjust charges, restoring to
Cr     given Madelung energy.  Minimizing (DQ)**2 subject to
Cr     constraints dE = 2 * V . DQ and sum DQ = 0 gives to first order
Cr     DQ = DE (V - Vbar) / (V - Vbar) . (V - Vbar)
Cr     DQ is calculated repeatedly until given emad is found
Cr   Case frzves=2:
Cr     is more stringent than frzves=1 because each total sphere
Cr     charge is determined by the input Madelung potential
Cr     (but shifted by a constant to make system neutral)
Cr   Case frzves=3:
Cr     calculates dq as in frzves=2, but qnu untouched
Cr   Case frzves=4:
Cr     calculates wk(j) such that sum wk(j) V(j) = total charge
Cr   Uses old ipc
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      implicit none
C Passed Parameters
      character*(*) clabl
      integer nbas,nl,nclass,nsp,ipc(*),nrclas(*),frzves,lmx(*)
      double precision ves(nclass),mad(nbas,nbas),wk(nclass,nclass),
     .  qnu(3,nl,nsp,nclass),pnu(nl,nsp,nclass),rmax(1),dq(1),z(1),
     .  emad,vrmax(1),rhrmx
C Local Parameters
      integer ib,ic,iclbas,n,ib0,ic0,il,j,isp,iprint,stdo,nglob
      double precision det,qca,amom,qtot,sumqsq,syschg,vconst,sum,
     .  emad0,vmtz(2),vdotv,vbar,trumad,etrms(22),xx
      integer, allocatable :: iwk(:)
      real(8), allocatable :: dwk(:),delq(:),qc(:)

      stdo = nglob('stdo')
      if (iprint() >= 20) write(stdo,*)
C --- Case frzves = 1 ---
      if (frzves == 1) then
      if (emad < 0) stop 'Vtoq: bad emad'
C --- Calculate Madelung potential and energy of this qnu  --
      allocate(delq(nclass),qc(nclass))
      call getq(nsp,nl,lmx,nclass,z,pnu,qnu,0,xx,qc,dq,delq)
      deallocate(delq,qc)

    5 continue
      call pshprt(iprint()-10)
      vmtz(1) = 0
      call madpot(nbas,1,nclass,nrclas,ipc,clabl,dq,0d0,rhrmx,rmax,mad,
     .  xx,xx,0d0,.false.,vrmax,ves,emad0,trumad,vmtz,0,xx,xx,etrms)
      call popprt

C --- Subtract off avg V, make (V - Vbar) . (V - Vbar) and dq --
      vbar = 0
      vdotv = 0
      do  ic = 1, nclass
        vbar = vbar + nrclas(ic)*ves(ic)
      enddo
      vbar = vbar/nbas
      do  ic = 1, nclass
        ves(ic) = ves(ic) - vbar
      enddo
      do  ic = 1, nclass
        vdotv = vdotv + nrclas(ic)*ves(ic)**2
      enddo
      call daxpy(nclass,(emad-emad0)/vdotv,ves,1,dq,1)

C --- Shift charge and repeat until emad-emad = 0 --
      if (dabs(emad0-emad) > 1d-8) goto 5

C --- Case frzves = 2,3,4 ---
      elseif (frzves == 2 .or. frzves == 3 .or. frzves == 4) then

C --- make the matrix M such that Q = M V --
      call dpzero(wk,nclass**2)
      do  ic0 = 1, nclass
        ib0 = iclbas(ic0-1,ipc,nbas)
        if (ib0 == 0) cycle
        do  ib = 1, nbas
          ic = ipc(ib)
          wk(ic0,ic) = wk(ic0,ic) + 2*mad(ib0,ib)
        enddo
        wk(ic0,ic0) = wk(ic0,ic0) + 2/rmax(ic0)
      enddo

      allocate(iwk(nclass),dwk(nclass))
      call dgefa(wk,nclass,nclass,iwk,n)
      if (n /= 0) stop 'Vtoq: madelung matrix singular'
      call dgedi(wk,nclass,nclass,iwk,det,dwk,1)
      deallocate(iwk,dwk)

C --- Find wk(ic) s.t. sum_ic wk(ic) v(ic) = tot Q ---
      if (frzves == 4) then
        do  ic = 1, nclass
          sum = 0
          do  j = 1, nclass
            sum = sum+wk(j,ic)*nrclas(j)
          enddo
          wk(ic,1) = sum
        enddo
        return
      endif

C --- Find the shift in ves to make system neutral ---
      call dmpy(wk,nclass,1,ves,nclass,1,dq,nclass,1,nclass,1,nclass)
      syschg = 0
      vconst = 0
      do  ic = 1, nclass
        sum = 0
        do  j = 1, nclass
          sum = sum+wk(ic,j)
        enddo
        vconst = vconst + nrclas(ic)*sum
        syschg = syschg+nrclas(ic)*dq(ic)
      enddo
      vconst = syschg / vconst
      do  ic = 1, nclass
        ves(ic) = ves(ic) - vconst
      enddo

C --- Find the q corresponding to shifted ves ---
      call dmpy(wk,nclass,1,ves,nclass,1,dq,nclass,1,nclass,1,nclass)

      if (frzves == 3) goto 10
      else
        return
      endif

C --- Adjust moments to fit dq ---
      if (nsp == 2) write(stdo,*) 'need check vtoq for nsp=2'

      if (iprint() > 20) write(stdo,*)
     .  'Class  New DQ     Old DQ   Difference      By l channel ...'

      do  ic = 1, nclass
        call getq(nsp,nl,nl-1,1,z(ic),pnu(1,1,ic),qnu(1,1,1,ic),0,xx,
     .    qca,qtot,amom)

        sumqsq = 0
        do  isp = 1, nsp
          do  il = 1, nl
            sumqsq = sumqsq + qnu(1,il,isp,ic)**2
          enddo
        enddo

        if (iprint() > 20) write(stdo,1) ic, dq(ic), qtot, dq(ic)-qtot,
     .  (((dq(ic)-qtot)*qnu(1,il,isp,ic)**2/sumqsq, isp=1,nsp), il=1,nl)
    1   format(i4,7f11.6:/37x,4f11.6)

        do  isp = 1, nsp
          do  il = 1, nl
          qnu(1,il,isp,ic) = qnu(1,il,isp,ic) +
     .      (dq(ic)-qtot)*qnu(1,il,isp,ic)**2/sumqsq
          enddo
        enddo
      enddo

   10 continue
      if (iprint() > 20 .and. (frzves >= 1.and.frzves <= 3))
     .  write(stdo,2) vconst,syschg
    2 format(' Vtoq:  pot adjusted by',f9.6,
     .       ' to adjust for net charge',f9.6)
      end
