      subroutine ioden(sopts,s_lat,s_site,s_spec,s_rhat,smrho)
C- File I/O charge density on a uniform mesh in a plane or full 3d mesh
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat nabc ng vol alat qlat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos kv gv
Cio    Passed to:  rhgcmp
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  rhgcmp rhogkl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z lmxl rg lfoca rfoca qc ctail etail stc lmxb p pz
Ci                 rmt a nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  rhgcmp corprm rhogkl
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2 rhoc
Cio    Passed to:  rhgcmp rhogkl
Ci Inputs
Ci   sopts :character string specifying plane and other options.
Ci         :Specifications and options are separated by a delimiter,
Ci         :which is the first character in sopts.
Ci         :
Ci         :if option is 3d then the density is written to
Ci         :disk on the full 3d grid.   Otherwise:
Ci         :The density is written to disk for a uniform of mesh of
Ci         :points in one plane.  This information is specified by three
Ci         :groups of numbers: the origin, a first direction vector with
Ci         :its number of points, and a second direction vector with its
Ci         :number of points.
Ci         :
Ci         :At present, these points must be taken from the points on
Ci         :the smooth mesh, smrho.  In this case, all three groups of
Ci         :information are sets of integers.  For example,
Ci         :specify the origin by three numbers:
Ci         :    o=#1,#2,#3
Ci         :The point (#1,#2,#3) corresponds to the Cartesian coordinate
Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
Ci         :where (n1,n2,n3) are the number of divisions in the
Ci         :mesh along the three lattice vectors (p1,p2,p3).
Ci         :o=0,0,0 corresponds to the origin.
Ci         :
Ci         :Specify the direction vectors by
Ci         :    l1=#1,#2,#3[,#4]
Ci         :    l2=#1,#2,#3[,#4]
Ci         :
Ci         :l1 and l2 specify the first and second direction vectors,
Ci         :respectively.  #1,#2,#3 select the
Ci         :increments in mesh points along each of the three lattice
Ci         :vectors that define the direction vector.  Thus in Cartesian
Ci         :coordinates a direction vector is
Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
Ci         :where as before (n1,n2,n3) are the number of divisions in
Ci         :the mesh along the three lattice vectors (p1,p2,p3).
Ci         :The last number (#4) specifies how many points to take
Ci         :in that direction.
Ci         :
Ci         :Other options:
Ci         :  3d        output of 3d grid plus headers in xsf format ready for xcrysden
Ci         :  fn=name   specifies the file name for file I/O
Ci         :
Ci         :  core=#    specifies how local rho is to be included
Ci         :            #=0 include core densities - nuclear charge
Ci         :            #=1 include core densities
Ci         :            #=2 (default) exclude core densities
Ci         :            #=-1 no local densities to be included
Ci         :            #=-2 true local density, no smoothed part
Ci         :            #=-3 istl-local sm densities, no true part
Ci         :
Ci         :Example: use '~' as the delimiter, and suppose
Ci         :n1=n2=48 and n3=120, the specification
Ci         :  ~fn=myrho~o=0,0,60~l1=1,1,0,49~l2=0,0,1,121
Ci         :writes 'myrho.ext' a mesh (49,121) points.
Ci         :The origin (first point) lies at (p3/2)
Ci         :The first vector points along (p1+p2), and has that length;
Ci         :the second vector points along p3, and has that length.
Ci   smrho :smooth density or potential on uniform mesh, depending on context
Co Outputs
Co   density mode: sum local gaussian densities + smrho is written to disk
Co   potential mode: vext is stored in smrho (name smrho is misleading!)
Cl Local variables
Cl   modrhg:controls what part of core density is added
Cl         : 2 exclude core densities (default)
Cl         : 1 include core densities, no nuclear contributions
Cl         : 0 include core density, including smoothed nuclear charge
Cl         :-1 smooth density only --- add no local densities
Cl         :-2 local density, with no smoothed part
Cl         :-3 interstitial and local smoothed densities
Cr Remarks
Cr    sopts specifies which plane(s) are written to disk
Cr Bugs
Cr   Routines create smoothed approximation to density, not true density
Cu Updates
Cu   30 Jun 16 New option to generate external potential
Cu   26 Apr 16 New option to write density at 1 kp
Cu   05 Sep 15 New option to write spin density
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   July 7 05 (Walter Lambrecht) New option g3d
Cu   25 Aug 04 New modes -2, -3
Cu   24 May 03 Corrections to errors in messages
Cu   23 Oct 01 Local densities are more accurately represented
Cu             in G_kL expansion:  k=0..kmax.  Added core= option.
Cu   25 Apr 01 Simplified the plane specification
Cu   02 Mar 01 Spin polarized
Cu   09 Feb 01 Added local gaussian densities to file I/O
Cu             and extended how a plane may be specified
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) sopts
      double complex smrho(*)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: z(:)
      complex(8), allocatable :: psrho(:,:,:,:)
      complex(8), allocatable :: cn(:,:)
      complex(8), allocatable :: wk(:)
C ... Local parameters
      logical lspin
      integer ngabc(3),n1,n2,n3,k1,k2,k3,kmax
      integer kkk,lgunit,ng,nglob,nsp,nwk,stdo,modrhg,nbas,i
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision vol,xx0,xx,plat(3,3),alat,xv(1)
      integer ib,is,isw

C ... External calls
      external awrit2,dcopy,defcc,defrr,fclose,fftz3,fftz30,gvgetf,
     .         gvputf,icopy,ioden2,ivset,mkils0,mkilst,nwordg,poppr,
     .         pshpr,rhgcmp,rhomom,rx,upack,upack2

C ... Unpack and setup
      nsp = nglob('nsp')
      nbas = nglob('nbas')
      stdo = lgunit(1)
      plat = s_lat%plat
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol
      alat = s_lat%alat
      allocate(z(nbas))
      do  ib = 1, nbas
        is = s_site(ib)%spec
        z(ib) = s_spec(is)%z
      enddo
      call fftz30(n1,n2,n3,k1,k2,k3)
      kkk = k1*k2*k3
      kmax = 3
      modrhg = 2
C     Call only parses argument list to returns modrhg and lspin
      call ioden2(0,s_lat,sopts,nsp,nbas,z,k1,k2,k3,
     .  xv,xv,nwk,modrhg,lspin,xv,0,[0],xv)
      call sanrg(.true.,modrhg,-3,2,'ioden','core option')

C ... Overwrite smrho+, smrho- with smrho, smrho+ - smrho-
      if (nsp == 2) then
        call dsumdf(kkk*2,1d0,smrho,0,1,smrho(1+kkk),0,1)
      endif

C ... Put n0(G) into psrho; translate into cn
      allocate(psrho(k1,k2,k3,nsp))
      allocate(cn(ng,nsp))
      call dcopy(kkk*2*nsp,smrho,1,psrho,1)
      call fftz3(psrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
      call gvgetf(ng,nsp,s_lat%kv,k1,k2,k3,psrho,cn)
      if (lspin) call dswap(ng*2,cn(1,1),1,cn(1,2),1)
      xx0 = cn(1,1)

C ... Add sum of local gaussian densities to mesh density
      if (modrhg >= 0 .or. modrhg <= -2) then
        if (modrhg ==  0) i = 131
        if (modrhg ==  1) i =  31
        if (modrhg ==  2) i =   1
        if (modrhg == -2) i =   2
        if (modrhg == -3) i =   3
        if (i == 2) call dscal(ng*nsp*2,0d0,cn,1)
        if (i == 3) call dscal(ng*nsp*2,-1d0,cn,1)
        if (lspin) i = i+100000
        call rhgcmp(i,1,nbas,s_site,s_spec,s_lat,s_rhat,kmax,ng,cn)
        if (i == 3) call dscal(ng*nsp*2,-1d0,cn,1)
      endif
      xx = cn(1,1)

C ... FFT (n0 + gaussians) (G) to real-space mesh
      if (modrhg >= 0) then
      call gvputf(ng,nsp,s_lat%kv,k1,k2,k3,cn,psrho)
      call awrit3('%N ioden : local '//
     .  '%?#n#spin ##'//
     .  'densities + envelope density, Qloc=%,6;6d  Q=%,6;6d',
     .  ' ',80,stdo,isw(lspin),(xx-xx0)*vol,xx*vol)
      if (modrhg < 2) then
        call awrit1('%9fLocal densities include core'//
     .    '%?#n==0#+nuclear## contributions.',' ',80,stdo,modrhg)
      endif
      elseif (modrhg == -2) then
        call info2(0,0,0,'%N ioden : local'//
     .    ' densities (true-smooth terms),  Qloc=%,6;6d',xx*vol,0)
      elseif (modrhg == -3) then
        call info2(0,0,0,'%N ioden : smooth (envelope - local)'//
     .    ' density,  Qs=%,6;6d  Qs-Qs(loc)=%,6;6d',
     .    xx0*vol,xx*vol)
      else
        call info2(0,0,0,'%N ioden : smooth density only'//
     .    ' (no local densities added) Qsm=%,6;6d',xx0*vol,0)
      endif

      call fftz3(psrho,n1,n2,n3,k1,k2,k3,nsp,0,1)

C ... File I/O
      nwk = 12*max(k1,k2,k3)
C     call defcc(owk,nwk**2*nsp)
      allocate(wk(nwk**2*nsp))
      call ioden2(2,s_lat,sopts,nsp,nbas,z,k1,k2,k3,
     .  psrho,wk,nwk,modrhg,lspin,[0d0],0,[0],xx)
      deallocate(wk,z,psrho,cn)

C ... Restore smrho+, smrho-
      if (nsp == 2) then
        call dsumdf(kkk*2,.5d0,smrho,0,1,smrho(1+kkk),0,1)
      endif
      end

      subroutine ioden2(mode,s_lat,sopts,nsp,nbas,z,k1,k2,k3,
     .  smrho,wk,nwk,modrhg,lspin,qp,nblst,iblst,eps)
C- Manipulation of mesh density/potential (kernel called by ioden or mkpot)
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc plat alat pos ng kv nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv ng kv ips0 bgv
Cio    Passed to:  symsmr
Ci Inputs
Ci   mode  :0 just return these arguments parsed from sopts:
Ci         :  modrhg, lspin
Ci         :1 just return these arguments parsed from sopts:
Ci         :  qp, nblst, iblst
Ci         :2 copy and save density in appropriate plane
Ci         :10 create external potential ('potential' mode). Output stored in smrho
Ci   sopts :string containing options; see ioden above.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :number of atoms in unit cell
Ci   z     :atomic number of atoms in cell
Ci   k1..3 :dimensions smrho
Ci   wk    :2D work array
Co         :Not touched in 'potential' mode
Ci   nwk   :dimensions wk
Co         :Not touched in 'potential' mode
Cio Inputs/Outputs
Cio  smrho :density mode: smooth density on uniform mesh (input)
Cio        :potential mode: potential (output)
Co   modrhg:controls what part of core density is added
Co         :not touched if mode=1 or ('potential') mode=10
Co   lspin :(all but potential mode) returned F, or T if tag 'spin' is encountered
Co         :(potential mode) returned F, or or T if tag 'pot0' is encountered
Co   nblst :number of states in lst=...
Co         :If on input, nblst is < 0,
Co         :input |nblst| sets an upper bound to the maximum size of list
Co         :not touched unless mode=1
Co   iblst :enumeration of states in lst=...
Co         :not touched unless mode=1
Cl Local variables
Cl   lic   :lic(1..3,1..2) = step length in each vector axes 1..2
Cl         :lic(1..3,3)    = (vext only) step length for vector normal to plane
Cl         :lic(4,1..2)    = number of points along axes 1..2
Cl         :lic(5..7,1)    = starting element for each lattice vector
Cl   vtype :0 => no potential specified
Cl         :1 => Specify at a lattice point (not implemented)
Cl         :2 => Specify at a G vector
Cl         :3 => constant potential in a sequence of planes
Cl   lzavg :T = subtract off average vext so <vext>=0 (potential mode)
Cs Command-line switches
Cs   fn=   :(density mode) fn=nam specifies that density is written to file nam
Cs   spin  :(density mode) If present, return lspin=T
Cs   3d    :(density mode) If present, write full 3D density
Cs   ro=.. :ro=x1,x2,x3 specifies the origin (actually the closest mesh point)
Cs   o=..  :o=i1,i2,i3  specifies the origin, as multiples of plat/ngabc
Cs         :Default in the absence of specification: (1,1,1) mesh point
Cs   q=..  :q=q1,q2,q3  specifies a single q for which density is made
Cs   lst=..:(density mode, use with q=), list-of-integers, e.g, lst=7,8. specifying QP levels
Cs         :Density from a specific set of states.  See Integer-list-syntax.html for list syntax.
Cs   l1=.. :l1=#1,#2,#3[,#4]  direction vector of 1st plane, as increments in mesh
Cs         :                  points along each of the three lattice vectors.
Cs         :                  #4 how many points to take on this line
Cs         :Default in the absence of specification: l1=1,0,0,ngabc(1)+1
Cs   l2=.. :l2=#1,#2,#3[,#4]  similar to l1=...  but for the second vector
Cs         :Default in the absence of specification: l2=0,1,0,ngabc(2)+1
Cs   core= :(density mode) specifies how local and core charges are to be added.
Cs         :See modrhg in ioden
Cs   vext=#:(potential mode) specifies potential. See Remarks.
Cs   v=#   :(potential mode) size of potential
Cs   step=# or igv=#,#,#:(potential mode) spcifies shape of potential; see Remarks
Cs   atpos :a special-purpose tag that connects site positions with
Cs         :points on the interstitial mesh.  Information is printed and program exits.
Cr Remarks
Cr   ioden is called in one of two modes:
Cr   *density mode, in one of two formats:
Cr    The projection of smrho onto a plane is written to disk
Cr    or if option is g3d the full 3d grid is written to disk
Cr
Cr    ioden2 has a special option to specify the contribution to density
Cr    from a single q-point.  ioden2 is called in a special setup mode
Cr    in which returns the k-point and possible list of QP states (mode=1)
Cr    The calling routine should generate the partial density for the
Cr    k-point and list of states (qp, nblst, iblst) returned by mode=1.
Cr
Cr   *potential mode:
Cr    An external potential vext is generated and written to smrho.
Cr    Specify it through --vext[~...]~v=#; # is the size of the potential.
Cr    At present the potential is specified by a rule which can assume one of two forms
Cr      a single Fourier component in reciprocal space
Cr      a set of planes in real space
Cr      In either case it is symmetrized by whatever symmetry operations are present.
Cr    Generate vext with one of the following:
Cr      ~step=#1       (real space) Make vext in a sequence of #1 adjacent planes
Cr      ~stepx=#1      (real space) Make vext in a sequence of #1 adjacent planes
Cr                     but all planes are shifted by a constant to make the average potential 0.
Cr      ~igv=#1,#2,#3  (Fourier component)  vext generated by a single Fourier component.
Cr                     (#1,#2,#3) is a fractional G vector specifed by the density mesh.
Cr    Other options:
Cr      ~eps=#         return # in argument eps (to tell the caller how to analyze the result)
Cr      ~eps           return NULLI in argument eps
Cr
Cr   For both density and potential modes:
Cr  *Planes are specified through switches l1= and l2=
Cr  *Origin is specifed through switches o= or ro=
Cr  *Plane is specified by two vectors, through switches l1= and l2=.
Cr   creates a default set of vectors, choosing the other lattice
Cr   vectors in the system.  Thus:
Cr   p1=# => lic(1..7,1) = (0 1 0 nb 1 1 1) and  lic(1..7,2) = (0 0 1 nc **)
Cr   p2=# => lic(1..7,1) = (0 0 1 nc 1 1 1) and  lic(1..7,2) = (1 0 0 na **)
Cr   p3=# => lic(1..7,1) = (1 0 0 na 1 1 1) and  lic(1..7,2) = (0 1 0 nb **)
Cr
Cr   Examples for density mode (delimiter ~)
Cr    *Given a density on a (48,48,120) mesh:
Cr        ~fn=myrho~o=0,0,60~l1=1,1,0,49~l2=0,0,1,121
Cr     The origin (first point) lies at (p3/2).
Cr     The first vector points along (p1+p2), and has that length.
Cr     The second vector points along p3, and has that length.
Cr    *Given a density on a (40,40,40) mesh in an fcc lattice, standard vectors
Cr      ~q=0,0,0.001~core=-1~ro=0,0,.25~lst=7~l1=-1,1,1~l2=1,-1,1
Cr     generates the smoothed part of the density from the 7th band at Gamma,
Cr     in a plane normal to the z axis, passing through (0,0,1/4).
Cr   Examples that generate potentials (delimiter ~)
Cr   *Create up a heaviside potential .01 Ry normal to z, passing through (0,0,1/4).
Cr    Potential given on a (40,40,40) mesh, standard vectors for an fcc lattice.
Cr     ~step=2~ro=0,0,.25~v=.01~l1=-1,1,1~l2=1,-1,1
Cu Updates
Cu   29 Jun 18 Bug fix (mode 10): ioden2 wrongly did not symmetrize potential
Cu   25 Jan 17 vext has new eps and pot0 options
Cu   26 Apr 16 New option in mode, and switches q= lst= ro=
Cu   09 Feb 01 Revised how a plane is defined.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) sopts
      logical lspin
      integer mode,k1,k2,k3,nwk,nsp,modrhg,nbas,nblst,iblst(nblst)
      double complex smrho(k1,k2,k3,nsp)
      double precision eps,wk(nwk,nwk,nsp),qp(3),z(nbas)
C ... For structures
!       include 'structures.h'
      type(str_lat)::   s_lat
C ... Pointers and dynamically allocated arrays
      integer, allocatable :: countp(:,:)
      real(8), allocatable :: xpos(:,:),xpost(:,:)
      real(8), pointer :: pos(:,:)
C ... Local parameters
      logical lok,lg3d,lvext,latpos,lzavg
      integer stdo,j2,j1,fopn,ifi,i,j,k,isp,i1,i2,i3,nblstl,mod0,mod1,vtype,nstep,iporigin
      integer lic(7,3),iv(3),jv(3),kv(4),igv(3),ngabc(3),nrpt(4)
      double precision alat,pspace,vv,vext
      double precision plat(3,3),qlat(3,3),qpl(3),pout(3),origin(3),vecs(3,4)
      character*120 dc*1, fn, prfmt*40, bndlst
      integer,parameter :: NULLI=-99999
      procedure(logical) :: isanrg
      procedure(integer) a2vec,lgunit,isw,iprint
      procedure(real(8)) ddot,dlength

      stdo = lgunit(1)
      dc = sopts(1:1)
      fn = 'smrho'
C     prfmt = '(8f10.6)'
      prfmt = '(1p8e14.6)'
      lg3d = .false.; latpos = .false.; lzavg = .false.
C     Default: origin at (1,1,1), (na+1,nb+1) points along 1st,2nd axes
      call iinit(lic,21)
C     call ivset(lic,5,7,1)
      lic(1,1) = 1
      nrpt(1) = ngabc(1)
      lic(4,1) = nrpt(1)+1
      lic(2,2) = 1
      nrpt(2) = ngabc(2)
      lic(4,2) = nrpt(2)+1
      qpl(1) = NULLI
      bndlst = ' '
      mod0 = mod(mode,10); mod1 = mod(mode/10,10)
      lvext = mod(mode/10,10) /= 0; vtype = 0 ! must be set by switch
      nstep = 0
      ngabc(1:3) = s_lat%nabc
      plat = s_lat%plat
      alat = s_lat%alat
      call dinv33(plat,1,qlat,vv)
      pos => s_lat%pos
      if (lvext) call dpzero(smrho,2*size(smrho))
      lspin = .false.
      vext = NULLI

      if (dc /= ' ') then
C ... Return here to resume parsing for arguments
      j2 = 0
   10 continue
      j2 = j2+1
      if (j2 <= len(sopts)) then
        if (sopts(j2:j2) == dc) goto 10
        j1 = j2
        call nwordg(sopts,0,dc//' ',1,j1,j2)
      else
        j1 = j2+1
      endif
      if (j2 >= j1) then
        if (.false.) then
C   ... switch fn
        elseif (.not. lvext .and. sopts(j1:j1+2) == 'fn=')  then
          if (j1+3 > j2) goto 99
          fn = sopts(j1+3:j2)
C   ... switch spin
        elseif (.not. lvext .and. sopts(j1:j1+3) == 'spin') then
          if (mod0 == 1) goto 10
          lspin = .true.
C   ... switch 3d
        elseif (.not. lvext .and. sopts(j1:j1+1) == '3d') then
          if (mod0 == 1) goto 10
          lg3d = .true.
C   ... switch ro=#,#,#
        elseif (sopts(j1:j1+2) == 'ro=')  then
          if (j1+2 > j2) goto 99
          if (mod0 == 1) goto 10
          i = j1+2
          if (a2vec(sopts,j2,i,4,', '//dc,3,2,3,kv,vecs) /= 3) goto 99
          call dgemm('T','N',3,1,3,1d0,qlat,3,vecs,3,0d0,vecs(1,2),3)
          do  i = 1, 3
            iv(i) = nint(ngabc(i)*vecs(i,2))
          enddo
          lic(5,1) = mod(iv(1)+ngabc(1),ngabc(1))
          lic(6,1) = mod(iv(2)+ngabc(2),ngabc(2))
          lic(7,1) = mod(iv(3)+ngabc(3),ngabc(3))
C   ... switch o=#,#,#
        elseif (sopts(j1:j1+1) == 'o=')  then
          if (j1+2 > j2) goto 99
          if (mod0 == 1) goto 10
          i = j1+1
          if (a2vec(sopts,j2,i,2,', '//dc,3,2,3,kv,iv) /= 3) goto 99
          lic(5,1) = mod(iv(1)+ngabc(1),ngabc(1))
          lic(6,1) = mod(iv(2)+ngabc(2),ngabc(2))
          lic(7,1) = mod(iv(3)+ngabc(3),ngabc(3))
C   ... switch q=#,#,#
        elseif (.not. lvext .and. sopts(j1:j1+1) == 'q=')  then
          if (j1+2 > j2) goto 99
          i = j1+1
          if (a2vec(sopts,j2,i,4,', '//dc,3,2,3,kv,qpl) /= 3) goto 99
          if (mod0 == 1) qp = qpl
C   ... switch lst=
        elseif (.not. lvext .and. sopts(j1:j1+3) == 'lst=')  then
          if (j1+4 > j2) call rx('ioden: bad list, lst=..')
          call mkils0(sopts(j1+4:j2),k,i)
          nblstl = k
          if (mod0 /= 1) goto 10
          if (nblst < 0 .and. k > iabs(nblst)) call rx('ioden: increase size of passed iblst')
          call mkilst(sopts(j1+4:j2),nblst,iblst)
C   ... switch l[12]=#,#,#,#
        elseif (sopts(j1:j1+2) == 'l1=' .or.
     .          sopts(j1:j1+2) == 'l2=')  then
          if (mod0 == 1) goto 10
          if (j1+3 > j2) goto 99
          i = 0
          call chrps2(sopts(j1+1:),'12',2,0,i,iv)
C         this check should never be necessary
C         call sanrg(.true.,iv,1,2,' ','iv in ioden ... bug ..')
          i = j1+2
          j = iv(1)
          k = a2vec(sopts,j2,i,2,', '//dc,3,2,4,kv,lic(1,j))
          if (k /= 3 .and. k /= 4) goto 99
          call ioden4(ngabc,lic(1,j),nrpt(j)) ! number of repeats needed to recover periodicity
          if (k == 3) lic(4,j) = nrpt(j)+1  ! And copy 1+nrpt to lic(4,j) if not set previously
C   ... switch core=#
        elseif (.not. lvext .and. sopts(j1:j1+4) == 'core=') then
          if (mod0 == 1) goto 10
          if (j1+5 > j2) goto 99
          i = j1+4
          if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,modrhg) /= 1) goto 99
C   ... switch pot0
        elseif (lvext .and. sopts(j1:j1+3) == 'pot0') then
          lspin = .true.
C   ... switch v=#
        elseif (lvext .and. sopts(j1:j1+1) == 'v=') then
          if (j1+2 > j2) goto 99
          i = j1+1
          if (a2vec(sopts,j2,i,4,' '//dc,2,1,1,kv,vext) /= 1) goto 99
C   ... switch eps=
        elseif (lvext .and. sopts(j1:j1+3) == 'eps=') then
          i = j1+3
          if (a2vec(sopts,j2,i,4,' '//dc,2,1,1,kv,eps) /= 1) goto 99
C   ... switch eps
        elseif (lvext .and. sopts(j1:j1+2) == 'eps') then
          eps = NULLI
C   ... switch stepx=#
        elseif (lvext .and. sopts(j1:j1+5) == 'stepx=') then
          lzavg = .true.
          vtype = 3
          if (j1+5+1 > j2) goto 99
          i = j1+5
          if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,nstep) /= 1) goto 99
C   ... switch step=#
        elseif (lvext .and. sopts(j1:j1+4) == 'step=') then
          vtype = 3
          if (j1+4+1 > j2) goto 99
          i = j1+4
          if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,nstep) /= 1) goto 99
C   ... switch igv=#,#,#
        elseif (lvext .and. sopts(j1:j1+3) == 'igv=') then
          vtype = 2
          if (j1+4 > j2) goto 99
          i = j1+3
          if (a2vec(sopts,j2,i,2,', '//dc,3,2,3,kv,igv) /= 3) goto 99
C   ... switch atpos : locate points on smooth mesh corresponding to site positions
        elseif (sopts(j1:j1+4) == 'atpos') then
          latpos = .true.
        else
          call rxs('ioden: unrecognised option ... ',sopts(j1:j2))
        endif
        goto 10
      endif
      endif
      if ((mod0==0 .or. mod0==1) .and. mod1/=1) return
      if (lvext .and. vtype == 0 .and. .not. latpos)
     .  call rx(' must specify potential through igv=#,#,# or through step=#')

C --- Printout, locate points on smooth mesh corresponding to site positions
      if (latpos) then
        call info2(10,1,0,
     .    ' Site positions as (possibly fractional) indices of %sx%3i smooth mesh',
     .    ngabc,0)
        allocate(xpos(3,nbas),xpost(nbas,3))
        do  i = 1, nbas
          call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,i),3,0d0,xpos(1,i),3)
          xpost(i,1:3) = ngabc(1:3)*xpos(1:3,i)
          do  j = 1, 3
            do while (xpost(i,j) < -1d-8)
              xpost(i,j) = xpost(i,j) + ngabc(j)
            enddo
          enddo
        enddo
        call arrprt(' Site  p1    p2    p3',
     .    '%,4i%;6,1D%;6,1D%;6,1D','Iddd',nbas,0,3,
     .    0,'  | ',[0],xpost(1,1),xpost(1,2),xpost(1,3),vv,vv,vv,vv)

        call info0(10,1,0,
     .    ' Smallest mesh density increments (p1,p2,p3) parallel to specified coordinates:'//
     .    '%N  vector%25frs lattice    repeat in   recip lattice')
        do  i = 1, 4
          vecs(:,1) = 0; vecs(min(i,3),1) = 1
C         Normal vector
          if (i == 4) then
            do  k = 1, 2
              vecs(:,k) = 0
              do  j = 1, 3
                vv = dble((lic(4,k)-1)*lic(j,k))/dble(ngabc(j))
                call dpadd(vecs(1,k),plat(1,j),1,3,vv)
              enddo
            enddo
            call cross(vecs(1,1),vecs(1,2),vecs(1,3)) ! Normal vector
            vecs(:,1) = vecs(:,3)
            call dscal(3,1/dlength(3,vecs,1),vecs,1) ! Normalize it
          endif
          call ioden3(1,plat,ngabc,6,vecs,iv,vecs(1,2))
          call ioden3(1,qlat,ngabc,6,vecs,jv,vecs(1,2))

C         Find smallest multiple of mesh (p1,p2,p3) that recovers a lattice vector
          k = maxval(ngabc)*max(4*9*5*7,iv(1)*iv(2)*iv(3))
          do  j = 1, k
            do  j2 = 1, 3
              if (ngabc(j2)*(j*iv(j2)/ngabc(j2)) /= j*iv(j2)) goto 20  ! Not an integer multiple of plat(j2)
            enddo
            j2 = j
            exit
   20       continue
            j2 = 0
          enddo
          if (sum(abs(iv(1:3))) /= 0) then
            call info5(10,0,0,' '//
     .        '%s,(%3;6d)'//
     .        '%35p%s,(%3i)'//
     .        '%50p%,3i'//
     .        '%?#n<5#%60p%s,(%3i)##',
     .        vecs,iv,j2,i,jv)
          endif
        enddo

        call fexit(0,0,'xx',0)

C --- Printout for 3D density ---
      elseif (lg3d) then        ! In xsf format ; see www.xcrysden.org
        ifi = fopn(fn)
        rewind ifi
        call awrit3('%9fWriting smooth density to file '//trim(fn)//
     .  ' : full 3d grid (%i,%i,%i).',' ',80,stdo,k1,k2,k3)
        do isp = 1, nsp
          write(ifi,'("CRYSTAL")')
          write(ifi,'("PRIMVEC")')
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .                           ,i2=1,3)
          write(ifi,'("PRIMCOORD")')
          write(ifi,'(2i5)') nbas,1
          do i = 1, nbas
            write(ifi,'(i4,2x,3f10.5)') int(z(i)),
     .           (pos(i2,i)*alat*0.529177208,i2=1,3)
          enddo
          write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
          write(ifi,'("charge_density_spin_",i1)') isp
          write(ifi,'("BEGIN_DATAGRID_3D_isp_",i1)') isp
          write(ifi,'(3i4)') k1,k2,k3
          write(ifi,'(3f10.5)') 0.,0.,0.
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .                           ,i2=1,3)
          write(ifi,'(1p8e14.6)')
     .      (((dble(smrho(i1,i2,i3,isp)),i1=1,k1),i2=1,k2),i3=1,k3)
          write(ifi,'("END_DATAGRID_3D_isp_",i1)') isp
          write(ifi,'("END_BLOCK_DATAGRID_3D")')
        enddo
        call fclose(ifi)
        return

C --- Potential for one Fourier component ---
      elseif (vtype == 2) then
        if (vext == NULLI) call rx('You must specify potential through v=#')
        vecs(1:3,1) = matmul(qlat,igv)
        call info5(10,1,0,
     .    ' Add external V_G = %g at igv = %s,(%3i)G -> (%3,6;6d)2*pi/alat',
     .    vext,igv,vecs(1,1),0,0)
C       Find index to this G
        do  i = 1, s_lat%ng
          if (dlength(3,vecs(1:3,1)-s_lat%gv(i,:),1) < 1d-6) then
            jv(1:3) = s_lat%kv(i,:)
            do  isp = 1, nsp
              smrho(jv(1),jv(2),jv(3),1:nsp) = vext
            enddo
            exit
          endif
          if (i == s_lat%ng) call rx('missing G vector missing from table')
        enddo
C       FT to real space mesh
        call fftz3(smrho,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,nsp,0,1)

C --- Density in plane or potential in nstep planes ---
      else
      if (lvext) then         ! Make external potential
        call info5(10,1,0,' Add external potential vext = %g in %i planes passing through origin at %s,(%3i)',
     .    vext,nstep,lic(5:7,1),0,0)
        if (vext == NULLI) call rx('You must specify potential through v=#')
      else
        call info2(10,0,0,' Writing smooth density to file '//fn//
     .    '%a : origin at %s,(%3i).',lic(5:7,1),0)
      endif
      call ioden6(0,lic(5,1),ngabc,plat,origin)
      call info2(10,0,0,' Origin, in cartesian coordinates %3:1,6;6d',origin,1)
      if (qpl(1) /= NULLI) then
        call awrit2(' Density generated from q =%3:1;6d for %i bands',' ',120,stdo,qpl,nblstl)
      endif

C     Printout info about vectors making up plane
      call dpzero(vecs,6)
      do  k = 1, 2
        do  i = 1, 3
          vv = dble((lic(4,k)-1)*lic(i,k))/dble(ngabc(i))
          call dpadd(vecs(1,k),plat(1,i),1,3,vv)
        enddo
        vv = dsqrt(ddot(3,vecs(1,k),1,vecs(1,k),1))
        call info5(10,0,0,' v%i: (%i+1 pts) = %s,(%3i)(p1,p2,p3) = (%3,6;6d) l=%,6;6d',k,lic(4,k)-1,lic(1,k),vecs(1,k),vv)
      enddo
      vv = acos(ddot(3,vecs(1,1),1,vecs(1,2),1)/dlength(3,vecs(1,1),1)/dlength(3,vecs(1,2),1))
      call info(10,0,0,' Angle between vectors : %d radians = %;1d deg',vv,180/(4*datan(1d0))*vv)

C      recip lattice of microcells
C      forall (i = 1:3) pmic(1:3,i) = plat(1:3,i)/ngabc(i)
C      call mkqlat(pmic,qmic,vv)

C     vecs(:,4) = vector normal to plane defined by lines 1 and 2
      call cross(vecs(1,1),vecs(1,2),vecs(1,4))  ! Normal vector
      call dscal(3,1/dlength(3,vecs(1,4),1),vecs(1,4),1)  ! Normalize it
C     Find (n1,n2,n3) that connects to closest neighboring plane and with the smallest length
      j2 = max(sum(abs(lic(1:3,1:2))),3)
      call ioden3(1,plat,ngabc,j2,vecs(1,4),lic(1,3),vecs(1,3)) ! smallest projection along normal
C     call ioden4(ngabc,lic(1,3),nrpt(4)) ! number of repeats needed to recover periodicity
      vv = ddot(3,vecs(1,3),1,vecs(1,4),1) ! Length of NN connecting vector normal to plane
      call ioden3(0,plat,ngabc,j2,vecs(1,4),lic(1,3),vecs(1,3))  ! Allow nonorthogonal => maybe pick up nearer planes
      pspace = ddot(3,vecs(1,3),1,vecs(1,4),1) ! Projection onto vector normal to plane
      call ioden4(ngabc,lic(1,3),nrpt(3)) ! number of repeats needed to recover periodicity
      nrpt(4) = vv/pspace*nrpt(3)
      vv = ddot(3,vecs(1,4),1,origin,1)/pspace ! Projection onto interplanar vector
      iporigin = nint(vv)
      if (abs(vv-iporigin) > 1d-8)
     .  call rx('bug in ioden: interplanar spacing not an integer multiple of p1,p2,p3')
      call info2(10,1,0,' Unit normal vector %s,(%3,6;6d).  Origin lines in plane %i',vecs(1,4),iporigin)
      call info5(10,0,0,' NN planes connected by %s,(%3i)(p1,p2,p3) = (%3,6;6d)',
     .  lic(1,3),vecs(1,3),0,0,0)
      call info5(10,0,0,' Planes separated by %g, %i planes to shortest period,'//
     .  ' %i planes in normal direction',pspace,nrpt(3),nrpt(4),0,0)

      if (lvext) then           ! Make external potential
        lic(4,3) = nstep
        if (nrpt(4) /= 0) then
          call info5(10,0,0,' Put vext in %i planes in interval [%i,%i] (%,1;2d %%)',
     .    lic(4,3),iporigin,iporigin+nstep-1,dble(100*nstep)/nrpt(4),0)
        else
          call info5(10,0,0,' Put vext in %i planes in interval [%i,%i]',
     .      lic(4,3),iporigin,iporigin+nstep-1,4,5)
        endif
      endif

C     Sanity checks
      lok = .true.
      lok = lok .and. lic(4,1) > 0
      lok = lok .and. lic(4,2) > 0
      if (.not. lok) call fexit2(-1,1,' Exit -1 ioden: number of '//
     .  'points along axes (%i,%i) are not > 0',lic(4,1),lic(4,2))
      if (.not. lvext) then
        if (nwk<lic(4,1) .or. nwk<lic(4,2)) call rx('increase nwk')
      endif

C --- Copy into vext into smrho for each point in planes 0 .. nstep-1 ---
      if (lvext) then
        allocate(countp(0:nstep-1,2)); countp = 0
C       forall (i = 1:3) lic(4,i) = min(lic(4,i),nrpt(i))  ! Limit increments to the number of repeats to periodicity
        do  i3 = 1, ngabc(3)
          do  i2 = 1, ngabc(2)
            do  i1 = 1, ngabc(1)

              call ioden6(1,[i1-1,i2-1,i3-1],ngabc,plat,pout) ! Candidate point r
              pout = pout - origin ! subtract origin o => v = r-o
              vv = ddot(3,vecs(1,4),1,pout,1)/pspace ! Projection onto interplanar vector
              if (vv < -1d-8 .or. nint(vv) >= nstep) cycle
              if (abs(vv-nint(vv)) > 1d-8)
     .          call rx('bug in ioden: interplanar spacing not an integer multiple of p1,p2,p3')
              countp(nint(vv),1) = countp(nint(vv),1) + 1

C              print 345, i1,i2,i3, pout+origin, pout, vv
C  345         format('i123 ',3i4,'  pos',3f12.6,'  rel to origin',3f12.6,'  proj',f12.6)

            do  isp = 1, nsp
              smrho(i1,i2,i3,isp) = vext
            enddo

          enddo
        enddo
      enddo

C     Find average vext; count number of points with nonzero vext
      i = 0; vv = 0
      do  i3 = 1, ngabc(3)
        do  i2 = 1, ngabc(2)
          do  i1 = 1, ngabc(1)
            vv = vv + smrho(i1,i2,i3,1)
            if (smrho(i1,i2,i3,1) /= 0) i = i+1
          enddo
        enddo
      enddo
      j = sum(countp); k = ngabc(1)*ngabc(2)*ngabc(3); vv = vv/k
      call info5(10,0,0,' vext added to %i points out of %i (%,1;2d %%)  <v>=%g'//
     .  '%?#n# (remove avg => <v>=0)##',j,k,dble(100*j)/k,vv,isw(lzavg))
      j = isw(isanrg(i,j,j,' ioden2 (warning)!','nonzero values of vext',.false.))
      if (lzavg) smrho = smrho - vv

C      if (iprint() >= 20) then
C        forall (j = 0:nstep-1) countp(j,2) = j+iporigin
C        call arrprt(' Plane  n','%,5i%,5i','ii',nstep,stdo,8,0,' | ',
C     .    countp(0,2),countp,0,0,0,0,0,0)
C      endif
      deallocate(countp)

C     call zprm3('vext(r)',0,smrho,k1,k2,k3)

C --- Copy points to wk from rho in plane ---
      else
        call icopy(3,lic(5,1),1,kv,1)   ! Initial point, to be shifted on 2nd axis
        do  i2 = 1, lic(4,2)            ! For each point on 2nd axis, do
          call icopy(3,kv,1,iv,1)       ! Initial point + 2nd axis shift
          do  i1 = 1, lic(4,1)          ! For each point on 1st axis, do
            jv = iv
            call ioden6(1,jv,ngabc,plat,pout) ! Shorten to equivalent point in unit cell
            do  isp = 1, nsp
              wk(i1,i2,isp) = dble(smrho(jv(1)+1,jv(2)+1,jv(3)+1,isp))
            enddo
            iv(1:3) = iv(1:3) + lic(1:3,1) ! Add 1st axis to increment position
          enddo
          kv(1:3) = kv(1:3) + lic(1:3,2) ! Add 2nd axis to increment position
        enddo

C   ... Write density to data file smrho
        ifi = fopn(fn)
        rewind ifi
        do  isp = 1, nsp
          call ywrm(0,' spin 1',1,ifi,prfmt,wk(1,1,isp),0,nwk,lic(4,1),lic(4,2))
        enddo
        call fclose(ifi)
        return
      endif                     ! lvext branch
      endif                     ! Density in plane or potential in nstep planes

C --- Cleanup (potential case) ---
      if (lvext) then         ! Symmetrize potential
        call symsmr(s_lat,nsp,0,k1,k2,k3,smrho)
C       Take real part of potential
        do  i3 = 1, ngabc(3)
          do  i2 = 1, ngabc(2)
            do  i1 = 1, ngabc(1)
              smrho(i1,i2,i3,1:nsp) = dble(smrho(i1,i2,i3,1:nsp))
            enddo
          enddo
        enddo
C       call zprm3('vext(r) before symmetrization',0,smrho,k1,k2,k3)
        call symsmr(s_lat,nsp,0,k1,k2,k3,smrho)

C       Debugging
C       call fftz3(smrho,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,nsp,0,-1)
C       call zprm3('vext(G) after symmetrization',0,smrho,k1,k2,k3)

        return
      endif

C --- Error handling ---
   99 continue
      call rxs('ioden: failed to parse option ... ',sopts(j1:j2))

      end
      subroutine ioden3(mode,plat,ngabc,nplat,vnorm,ivdotn,vdotn)
C- Find combination of plat that best approximates given normal vector
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 find smallest vector that connecting nearest neighboring plane
Ci         :1 find smallest vector parallel to given vnorm
Ci   plat  :primitive lattice vectors, in units of alat
Ci   ngabc :number of divisions in smooth mesh
Ci   nplat :number of multiples of plat to make excursions over
Ci   vnorm :normal vector (normal to plane, mode 0)
Co Outputs
Co  ivdotn :smallest vector, multiples of the mesh vectors p1,p2,p3
Co   vdotn :smallest vector
Cu Updates
Cu   12 Jun 16 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nplat,ngabc(3),ivdotn(3)
      double precision plat(3,3),vnorm(3),vdotn(3)
C ... Local parameters
      logical candidate
      integer i1,i2,i3,i
      double precision vdoti(3),vn,vcut(2)
      procedure(real(8)) :: ddot,dlength

      ivdotn = 0
      vcut(:) = 9d9

      do  i1 = -nplat, nplat
        do  i2 = -nplat, nplat
          do  i3 = -nplat, nplat
            do  i = 1, 3
              vdoti(i) = (plat(i,1)/ngabc(1))*i1 + (plat(i,2)/ngabc(2))*i2 + (plat(i,3)/ngabc(3))*i3
            enddo
            vn = ddot(3,vdoti,1,vnorm,1) ! Projection onto vector normal to plane
C           Vector is a candidate if its projection is smaller than or equal to any prior vector
            if (mode == 0) then
              candidate = vn > 1d-6 .and. vn < vcut(1)+1d-6
C           Vector is a candidate if it has 100% projection onto vnorm
            else
              if (dlength(3,vdoti,1) == 0.0_8) then
                candidate = vn == 0.0_8
              else
                candidate = abs(vn/dlength(3,vdoti,1)-1) < 1d-6
              end if
            endif
            if (candidate) then
              vcut(1) = vn ! spacing between plane
C             Update vector if its length is smaller than prior candidate
              vn = dlength(3,vdoti,1)
C             print *, sngl(vn), i1,i2,i3, sngl(vdoti)
              if (vn < vcut(2)) then
                vcut(2) = vn         ! Length of connecting vector
                ivdotn(1:3) = [i1, i2, i3]
                vdotn = vdoti
              endif
            endif
          enddo
        enddo
      enddo
      end

      subroutine ioden4(ngabc,n123,nrpt)
C- Find number of repeats for line connecting points in plat to recover periodicity
C ----------------------------------------------------------------------
Ci Inputs
Ci   ngabc :number of divisions (N1,N2,N3) of each of the three lattice vectors
Ci   n123  :(n1,n2,n3) in Remarks
Co Outputs
Co   nrpt  :smallest multiple of (n1,n2,n3) that makes (n1/N1,n2/N2,n3/N3) all integers
Cr Remarks
Cr  Vector (A,A') given by (n1/N1) p1 + (n2/N2) p2 + (n3/N3) p3
Cr  Find nrpt = smallest integer k/i that renders k/i * (n1/N1,n2/N2,n3/N3) all integers
Cr        *  *  *  *  *A'
Cr       *--*--*--*--*
Cr      *  *  *  *  *
Cr     *  *  *  *  *
Cr   A*  *  *  *  *
Cr    --*--*--*--*
Cr
Cu Updates
Cu   29 Jun 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngabc(3),n123(3),nrpt
C ... Local parameters
      integer i,k,kv(3)
      procedure(integer) comfac

      k = ngabc(1) * ngabc(2) * ngabc(3)
      kv(1) = n123(1) * k / ngabc(1)
      kv(2) = n123(2) * k / ngabc(2)
      kv(3) = n123(3) * k / ngabc(3)
      i = comfac([2,3,5,7,11,13,17,19],8,kv,3)
      nrpt = k / i
      end

      subroutine ioden6(mode,iv,ngabc,plat,pos)
C- Return position in Cartesian coordinates, and possibly shortened fractional coordinates
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 make pos only
Ci         :1 Shorten pos (and also iv), bringing position to within
Ci         :  first quadrant of unit cell at origin
Ci   ngabc :number of divisions of the lattice vectors
Ci   plat  :primitive lattice vectors, in units of alat
Co Input/Outputs
Cio  iv    :position as fractions of plat, fraction of jth vector is iv(j)/ngabc(j)
Cio        :iv will by modified is pos is shortened.
Co Outputs
Co   pos   :position given by iv, in Cartesian coordinates, possibly shortened
Cu Updates
Cu   30 Jun 16
C ----------------------------------------------------------------------
      implicit none
      integer mode,iv(3),ngabc(3)
      double precision pos(3),plat(3,3)
      integer i
      double precision pout(3)

      pos = 0
      call dpadd(pos,plat(1,1),1,3,dble(iv(1))/dble(ngabc(1)))
      call dpadd(pos,plat(1,2),1,3,dble(iv(2))/dble(ngabc(2)))
      call dpadd(pos,plat(1,3),1,3,dble(iv(3))/dble(ngabc(3)))

      if (mod(mode,10) == 1) then
        call shorps(-1,plat,[1,1,1],pos,pout)
        if (abs(pout(1))+abs(pout(2))+abs(pout(3)) > 1d-8) then
          forall (i = 1:3) iv(i) = iv(i) + nint(pout(i))*ngabc(i)
          if (iv(1) < 0 .or. iv(1) < 0 .or. iv(3) < 0) call rx('bug in ioden6')
          pos = 0
          call dpadd(pos,plat(1,1),1,3,dble(iv(1))/dble(ngabc(1)))
          call dpadd(pos,plat(1,2),1,3,dble(iv(2))/dble(ngabc(2)))
          call dpadd(pos,plat(1,3),1,3,dble(iv(3))/dble(ngabc(3)))
        endif
      endif

      end
