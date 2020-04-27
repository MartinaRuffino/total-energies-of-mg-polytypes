      subroutine rlxstp(s_ctrl,s_site,natrlx,nvar,indrlx,xyzfrz,pdim)
C- Set up variables for relaxation
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nbas nitmv mdprm defm ltb lfrce
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: relax
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs:
Co Outputs:
Co   indrlx(1,i) points to the ith relaxing component and
Co   indrlx(2,i) points to the corresponding site
Co   natrlx      # of relaxing degrees of freedom for atoms
Co   nvar:       # of variables to relax (natrlx)
Co   xyzfrz(i)   T means all the ith components are frozen (T on input)
Co   pdim:       dimension of the work array p, needed in relax
Cr Remarks:
Cr   At version 6.15 first attempt to restore volume and shear relaxations
Cr   At version 6.10 volume and shear relaxations have been removed;
Cr   hence nvar=natrlx.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   21 Mar 06 mdprm(1)>100 signifies shear relaxation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical xyzfrz(3)
      integer indrlx(2,*),nvar,natrlx,nitrlx,pdim
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
C ... Local parameters
      double precision mdprm(7),defm(6)
C     double precision autime,autmp
      integer nbas,i,j,k,iprint,i1mach,ltb,ifrlx(3),lrlx
      logical force,md,bittst
C     character*1 comp(3)
C     data comp /'x','y','z'/
C     data autime/0.048377d0/
C     data autmp/0.6333328d-5/

      nbas = s_ctrl%nbas
      nitrlx = s_ctrl%nitmv
      mdprm = s_ctrl%mdprm
      defm = s_ctrl%defm
      ltb = s_ctrl%ltb
      nvar = 0
      force = bittst(ltb,16) .or. s_ctrl%lfrce > 0
      if (.not. force .or. nint(mdprm(1)) == 0) return
      md = nint(mdprm(1)) <= 3
      lrlx = mod(nint(mdprm(1)),100)

C --- Set relaxation variables ---
      j = 0
      if (md) then
        xyzfrz(1) = .false.
        xyzfrz(2) = .false.
        xyzfrz(3) = .false.
        return
      elseif (force .and. mdprm(1) >= 100) then
        do  i = 1, 6
          if (defm(i) == 1) then
            j = j+1
            indrlx(1,j) = i
          endif
        enddo
      elseif (force) then
        do  i = 1, nbas
          ifrlx = s_site(i)%relax
          do  k = 1, 3
            if (ifrlx(k) == 1) then
              j = j + 1
              indrlx(1,j) = k
              indrlx(2,j) = i
              xyzfrz(k) = .false.
            endif
          enddo
      enddo
      endif
      natrlx = j
      nvar = natrlx
      if (nvar == 0) return
      pdim = 0
      if (.not. md) then
        if (lrlx == 4) pdim = nvar*7
        if (lrlx == 5) pdim = nvar*(7+nvar)
        if (lrlx == 6) pdim = nvar*(12+2*nvar)
      endif

C --- Printout ---
      if (iprint() >= 30) then
        if (lrlx == 4) then
          call info(0,1,0,
     .        ' RLXSTP: Molecular statics (conjugate gradients) ..',0,0)
        elseif (lrlx == 5) then
          call info(0,1,0,
     .        ' RLXSTP: Molecular statics (Fletcher-Powell) ..',0,0)
        else
          call info(0,1,0,
     .        ' RLXSTP: Molecular statics (Broyden) ..',0,0)
        endif
        call info2(0,0,0,
     .      '         relaxing %i variables, %i iterations',nvar,nitrlx)
        call awrit4('         x-tol=%d, g-tol=%d, step=%d (pdim=%i)',
     .              ' ',120,i1mach(2),mdprm(3),mdprm(4),mdprm(5),pdim)
      endif

      end
