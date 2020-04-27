      subroutine iodosw(ifi,dos,sigmrq,nemx,ndos,nchan,emin,
     .  emax,nspin,eferm)
C- Write weighted density-of-states
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :>0 => read (not allowed here), <0 => write
Ci   nemx  :leading dimension of dos array (must be >ndos)
Ci         :Not referenced unless mode>=2
Ci   ndos  :number of energy mesh points
Cio Inputs/Outputs
Cio  The following header is read from or written to disk
Cio  emin  :dos ranges from (emin,emax)
Cio  emax  :dos ranges from (emin,emax)
Cio  nchan :number dos channels per spin
Cio  nspin :number of spin channels
Cio  eferm :Fermi level
Cio  del   :Sampling delta
Cio  After the header, this is read from or written to disk
Cio  dos   :density of states ... not referenced unless mode>=2
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   18 Mar 10 Adapted from dosio
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,ndos,nemx,nchan,nspin
      double precision eferm,emax,emin
      double precision dos(nemx,*),sigmrq(nemx,*)
C ... Local parameters
      integer ie,ild
      double precision de,omg,wt(nchan*nspin)

C --- Write branch ---
      if (ifi < 0) then
        write(-ifi,763) ndos, nchan*nspin+1, eferm
        de = (emax-emin)/(ndos-1)
        do  ie = 1, ndos
          do  ild = 1, nchan*nspin
            if (sigmrq(ie,ild) /= 0) then
              wt(ild) = dos(ie,ild)/dsqrt(sigmrq(ie,ild))
            else
              wt(ild) = 0
            endif
          enddo
          omg = emin + de*(ie-1)
          write(-ifi,764) omg,(wt(ild),ild=1,nchan*nspin)
        enddo
  763   format('% rows',i6,'  cols',i5,'  eferm=',f11.6)
  764   format(f12.5,6f13.6:/(12x,6f13.6))
      endif

C --- Read branch ---
      if (ifi > 0) then
        call rx('iodosw has no read branch')
      endif
      end
