      subroutine ioqpp(lio,s_ctrl,s_pot)
C- File I/O for phi-phi, phi-dot, dot-dot products
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nbas nl nspin
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:qpp
Cio    Passed to: *
Ci Inputs:
Ci   lio: true for write, false for read
Ci          <0 write
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   07 Aug 10 Make MPI compatible
Cu   08 Nov 07 (J. Xu) qpp is complex
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lio
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
C ... Local parameters
      integer i1,i2,nbas,nl,nsp,ifi,fopna,rdm,ipr
      integer procid,mpipid

      call getpr(ipr)
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      i1 = nl**2*(nl**2+1)  ! /2 * 2 for complex
      i2 = 4*nsp*nbas

C ... MPI: only master does file I/O
      procid = mpipid(1)
      if (procid == 0) then

      ifi = fopna('qpp',-1,4+8)
      if (lio) then
        call ywrm(1,'lmasa',1,ifi,' ',s_pot%qpp,1,i1,i1,i2)
        if (ipr >= 30) print *, 'IOQPP:  wrote qpp to disk'
      else
        if (ipr >= 30) print *, 'IOQPP:  reading qpp from disk ...'
        call pshpr(1)
        if (rdm(ifi,2,i1*i2,' ',s_pot%qpp,i1,i2) < 0) then
          if (ipr >= 0)
     .      print *,'IOQPP:  (warning) failed to read qpp file'
          call dvset(s_pot%qpp,1,1,-1d0)
        endif
        call poppr
      endif
      call fclose(ifi)
      endif

      if (.not. lio) then
        call mpibc1(s_pot%qpp,i2*i2,4,.false.,'ioqpp','qpp')
      endif

      end
