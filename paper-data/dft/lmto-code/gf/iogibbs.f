      subroutine rdgibbs(nclass,nth,maxth,lweiss,temp,theta,wtq,wtg)
C- Reads or computes Gibbs weights and records them in wtg
C ----------------------------------------------------------------------
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   nth   :number of DLM angles
Ci   nth   :array with number of DLM angles for all classes (0 if no DLM)
Ci   maxth
Ci   lweiss
Ci   temp
Ci   theta
Ci   wtq
Ci   wtg
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lweiss
      integer nclass,maxth,nth(*)
      double precision wtq(*),wtg(*)
      double precision temp,theta(*)
C ... Local parameters
      double precision th(maxth),wt(maxth),znorm
      integer ic,i,ioff,ifi,fopn

      if (.not.lweiss) then
        ifi = fopn('gibbs')
        rewind ifi
        ioff = 1
        do 10 ic = 1, nclass
C ...     Normal sites have weight 1
          if (nth(ic) <= 1) then
            wtg(ic) = 1d0
            goto 10
          endif
C ...     Dummy sites have weight 0
          wtg(ic) = 0d0
          read(ifi,930,err=8,end=8) (wt(i),i=1,nth(ic))
          goto 9
C ...     If weights not found, set unnormalized weights to to 1
  8       do i = 1,nth(ic)
            wt(i) = 1d0
          enddo
  9       write(6,931) ic, (wt(i),i=1,nth(ic))
          call dpscop(wt,wtg,nth(ic),1,nclass+ioff,1d0)
          ioff = ioff + nth(ic)
 10     continue
        call fclose(ifi)
      else
        ioff = 1
        do 20 ic = 1,nclass
          if (nth(ic) <= 1) then
            wtg(ic) = 1d0
            goto 20
          endif
          wtg(ic) = 0d0
          call dpscop(theta,th,nth(ic),ioff,1,1d0)
          do i = 1,nth(ic)
            wt(i) = dexp(-temp*dcos(th(i)))
          enddo
          write(6,931) ic, (wt(i),i=1,nth(ic))
          call dpscop(wt,wtg,nth(ic),1,ioff+nclass,1d0)
          ioff = ioff + nth(ic)
 20     continue
      endif
C --- Normalizing the weights and recording them
      ioff = nclass
      do 30 ic = 1, nclass
        if (nth(ic) <= 1) goto 30
        znorm = 0d0
        do i = 1,nth(ic)
          wtg(ioff+i) = wtg(ioff+i) * wtq(ioff-nclass+i)
          znorm = znorm + wtg(ioff+i)
        enddo
        do i = 1,nth(ic)
          wtg(ioff+i) = wtg(ioff+i)/znorm
        enddo
        write(6,932) ic,(wtg(ioff+i),i=1,nth(ic))
        ioff = ioff + nth(ic)
  30  continue
      write(6,933) (wtg(i),i=1,ioff)
 930  format(100F13.8)
 931  format(' Gibbs weights for class',i3,': '/100(8F10.6/))
 932  format(' Total weights for class',i3,': '/100(8F10.6/))
 933  format(' Weights for all classes:'/100(8F10.6/))
      end

