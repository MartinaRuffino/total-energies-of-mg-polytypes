      subroutine vsh2es(mode,nbas,ipc,nrc,vshft,ves)
C- Adds/subtracts vshft from ves to produce true effective potentials
C ----------------------------------------------------------------
C  mode   1s digit   0 add vshft to ves   1 subtract from ves
C        10s digit   1 initialize input ves to zero
C       100s digit   1 reverse role of ves and vshft
C ----------------------------------------------------------------
      implicit none
      integer mode,nbas,ipc(nbas),nrc(nbas)
      double precision vshft(-7:nbas),ves(nbas),scal
      logical mode2
      integer ib,ic,iprint,nclass,stdo,nglob

C ... Initialize
      mode2 = mod(mode/100,10) == 1
      if (mod(mode/10,10) == 1)  then
        nclass = 0
        do  ib = 1, nbas
        ic = ipc(ib)
        nclass = max(nclass,ic)
        if (.not. mode2) ves(ic) = 0
        if (mode2) vshft(ib) = 0
        enddo
      endif

C ... Scale +/- 1
      scal = 1
      if (mod(mode/1,10) == 1) scal = -1

C ... Do the copy
      do  ib = 1, nbas
        ic = ipc(ib)
        if (.not. mode2) ves(ic) = ves(ic) + scal*vshft(ib)/nrc(ic)
        if (mode2) vshft(ib) = vshft(ib) + scal*ves(ic)
      enddo

      if (iprint() >= 30) then
        stdo = nglob('stdo')

        if (.not. mode2) write (stdo,1) 'Class'
        if (mode2) write (stdo,1) ' Site'
        if (mode2) then
          do  ic = 1, nclass
            write (stdo,2) ic,ves(ic)
          enddo
        else
          do  ib = 1, nbas
            write (stdo,2) ib,vshft(ib)
          enddo
        endif
      endif
    1 format(/1x,a,3x,'Shifted V')
    2 format(i4,f14.6)

      end
