      subroutine mkbfld(nl,nbas,lihdim,indxsh,bsite,nbf,eula,neul,bdots)
C- Generate B . sigma
C ----------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   bsite :magnetic field, site representation
Ci   nbf   :1 if bfield is l-independent, nl if l-dependent
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Co Outputs
Co   bdots :(magnetic field . Pauli matrix), in downfolding order
Cr Remarks
Cr   Routine makes B . (sigma = Pauli matrices):
Cr   This is the potential acting on electrons (spin 1/2), since:
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   Electron spin operator is S = 1/2 sigma or (hbar/2 sigma)
Cr   Coupling is g0 S B => E = (g0.ms.muB) B
Cr   Splitting between ms=1/2 and ms=-1/2 is  g0 muB B = 2 muB B.
Cr
Cr   Mapping of composite RL index to this vector is same as that of
Cr   eigenvectors permutated according to downfolding rules in indxsh.
Cu Updates
Cu   07 Jan 11 Fix erroneous 1/2 scaling of B
Cu             and use same sign convention for B as for Bxc.
Cu   27 Mar 04 bfield can be l- or lm- dependent
Cu   12 Feb 03 First created
Cb Bugs
Cb   makpph assumes that all m's corresponding to a given l
Cb   are folded down in the same way.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nbas,nbf,neul,lihdim,indxsh(lihdim)
      double precision bsite(nbas,nbf,3),eula(nbas,neul,3)
      double complex bdots(2,2,lihdim)
C ... Local parameters
      logical lwrite
      integer ib,l,m,lmr,ilmr,ipr,stdo,lgunit,ilm,nlm
      double precision B(3),bloc(3),rotm(3,3),facB,alpha,beta,gamma
C     facB covers convention for B:
C     facB= 1 => +B induces positive -M
C     facB=-1 => +B induces positive +M
      parameter (facB=1d0)

      call getpr(ipr)
      stdo = lgunit(1)

C ... Number of independent local B field channels
      if (nbf == nl*nl .or. neul == nl*nl) then
        nlm = nl*nl
      elseif (nbf > 1 .or. neul > 1) then
        nlm = nl
      else
        nlm = 1
      endif

C ... Printout header
      call info2(30,1,0,
     .  '%?#n==1#%17f#%20f#(global)%9fApplied B-field%11f(local)',nlm,0)
      call info5(30,0,0,
     .  '   ib%?#n==0#%4f##%?#n==0#  l    ##%?#n==0# ilm   ##'//
     .  'bx%8fby%8fbz%20fbx%8fby%8fbz',nlm-1,nlm-nl,nlm-nl**2,0,0)

      lmr = 0
      do  ib = 1, nbas

C   ... Magnetic field for site ib in global coordinates
        B(1) = bsite(ib,1,1)
        B(2) = bsite(ib,1,2)
        B(3) = bsite(ib,1,3)
        alpha = eula(ib,1,1)
        beta  = eula(ib,1,2)
        gamma = eula(ib,1,3)

C   ... Rotate magnetic field to local coordinates
        call eua2rm(alpha,beta,gamma,rotm)
C       call prmx('rotm',rotm,3,3,3)
        call dgemm('N','N',3,1,3,1d0,rotm,3,B,3,0d0,bloc,3)
C       call prmx('bfield, global coordinates',B,3,3,1)
C       call prmx('bfield, local  coordinates',bloc,3,3,1)

        lwrite = .true.
        ilm = 0
        do  l = 0, nl-1

          if (indxsh(lmr+1) > lihdim) then
            lmr = lmr + 2*l+1
            ilm = ilm + 2*l+1
            goto 2
          endif

          do   m = -l, l

            ilm = ilm+1
            lmr = lmr+1
            ilmr = indxsh(lmr)

C       ... If Euler angles are l or lm dependent
            if (neul == nl) then
              alpha = eula(ib,l+1,1)
              beta  = eula(ib,l+1,2)
              gamma = eula(ib,l+1,3)
            elseif (neul == nl*nl) then
              alpha = eula(ib,ilm,1)
              beta  = eula(ib,ilm,2)
              gamma = eula(ib,ilm,3)
            elseif (neul /= 1) then
              call rxi('mkbfld: bad value neul=',neul)
            endif

C       ... If b-field is l or lm dependent
            if (nbf == nl) then
              B(1) = bsite(ib,l+1,1)
              B(2) = bsite(ib,l+1,2)
              B(3) = bsite(ib,l+1,3)
            elseif (nbf == nl*nl) then
              B(1) = bsite(ib,ilm,1)
              B(2) = bsite(ib,ilm,2)
              B(3) = bsite(ib,ilm,3)
            elseif (nbf /= 1) then
              call rxi('mkbfld: bad value nbf=',nbf)
            endif

C       ... l-dependent fields : remake bloc
            if (neul /= 1 .or. nbf /= 1) then
              call eua2rm(alpha,beta,gamma,rotm)
              call dgemm('N','N',3,1,3,1d0,rotm,3,B,3,0d0,bloc,3)
            endif

            bdots(1,1,ilmr) =  bloc(3)                   *facB
            bdots(1,2,ilmr) =  dcmplx(bloc(1),-bloc(2))  *facB
            bdots(2,1,ilmr) =  dcmplx(bloc(1), bloc(2))  *facB
            bdots(2,2,ilmr) = -bloc(3)                   *facB

C       ... Printout
            if (nlm == nl**2 .and. ipr >= 30) then
              write(stdo,'(i5,i3,3f10.6,3x,8x,1x,3f10.6)')
     .            ib, ilm, B, bloc
            elseif (nlm == nl .and. ipr >= 30) then
              if (lwrite)
     .        write(stdo,'(i5,i3,3f10.6,3x,8x,1x,3f10.6)') ib,l,B,bloc
              lwrite = .false.
            elseif (ipr >= 30) then
              if (lwrite)
     .          write(stdo,'(i5,3f10.6,3x,8x,1x,3f10.6)') ib, B, bloc
              lwrite = .false.
            endif

C           debugging
C            if (l /= 2) then
C            bdots(1,1,ilmr) = 0
C            bdots(1,2,ilmr) = 0
C            bdots(2,1,ilmr) = 0
C            bdots(2,2,ilmr) = 0
C            endif

          enddo
          if (nlm >= nl) lwrite = .true.
    2     continue
        enddo
      enddo

C     call zprm('B.sigma',2,bdots,4,4,lihdim)

      end

      subroutine prbfield(nl,nbas,offH,bsite,nbf)
C- Prints out B-fields
C ----------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   bsite :magnetic field, site representation
Ci   nbf   :1 if bfield is l-independent, nl if l-dependent
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Co Outputs
Co   bsite is printed to stdout
C ----------------------------------------------------------------
      implicit none
C ... Passed variables
      integer nkap0,n0H,nl,nbas,nbf
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,*)
      double precision bsite(nbas,nbf,3)
C ... Local parameters
      logical lwrite
      integer ib,l,ipr,stdo,lgunit,ilm,nlm
      double precision B(3),bloc(3),facB
C     facB covers convention for B:
C     facB= 1 => +B induces positive -M
C     facB=-1 => +B induces positive +M
      parameter (facB=1d0)

      call getpr(ipr)
      stdo = lgunit(1)

C ... Number of independent local B field channels
      if (nbf == nl*nl) then
        nlm = nl*nl
      elseif (nbf > 1) then
        nlm = nl
      else
        nlm = 1
      endif
      if (nlm /= 1) call rx('mkbsz not ready for l-dependent field')

C ... Printout header
      call info2(30,1,0,
     .  '%?#n==1#%17f#%20f#(global)%9fApplied B-field%11f(local)',nlm,0)
      call info5(30,0,0,
     .  '   ib%?#n==0#%4f##%?#n==0#  l    ##%?#n==0# ilm   ##'//
     .  'bx%8fby%8fbz%20fbx%8fby%8fbz',nlm-1,nlm-nl,nlm-nl**2,0,0)

      do  ib = 1, nbas

C   ... Magnetic field for site ib in global coordinates
        B(1) = bsite(ib,1,1)
        B(2) = bsite(ib,1,2)
        B(3) = bsite(ib,1,3)
        bloc = B

C  ... Printout
        lwrite = .true.
        if (nlm == nl**2 .and. ipr >= 30) then
C         ! Oops ilm is not defined
          write(stdo,'(i5,i3,3f10.6,3x,8x,1x,3f10.6)')
     .      ib, ilm, B, bloc
        elseif (nlm == nl .and. ipr >= 30) then
C         ! Oops l is not defined
          if (lwrite)
     .      write(stdo,'(i5,i3,3f10.6,3x,8x,1x,3f10.6)') ib,l,B,bloc
          lwrite = .false.
        elseif (ipr >= 30) then
          if (lwrite)
     .      write(stdo,'(i5,3f10.6,3x,8x,1x,3f10.6)') ib, B, bloc
          lwrite = .false.
        endif

        if (nlm >= nl) lwrite = .true.
      enddo

C     call zprm('B.sigma',2,bdots,4,4,lihdim)

      end
