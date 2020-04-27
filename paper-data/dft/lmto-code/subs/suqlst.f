      subroutine suqlst(s_lat,sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,qp,onesp)
C- Set up a list of q-points in various modes for energy bands
C  See also entry suqlsw(nband,jsp,nsp,evl)
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :character string describing special options that can
Ci         :selects the mode for which qp are read.
Ci         :*default mode: the qp file consists of a list of
Ci         :               lines and the number of qp in each line.
Ci         :*list mode   : the qp file consists of a list of qp
Ci         :*contour mode: input for specifying qp on a uniform 2D
Ci                       : mesh, for contour plots.
Ci         :See Remarks for the syntax of each mode.
Ci
Ci         :Options are separated by delimiters; the first character is
Ci         :the delimiter.  The following list the strings declaring
Ci         :options assuming the delimiter is '~'. (Space as the first
Ci         :implies that there are no options)
Ci         :
Ci         :suqlst generates qp in one of 4 modes.
Ci         :mode 1 qp are generated along a user-supplied list of symmetry lines
Ci         :       This is the default.
Ci         :mode 2 qp are generated from a user-supplied list
Ci         :mode 3 qp are generated on a regular mesh in a plane
Ci         :mode 4 qp are generated in a shell centered at a user-specified origin
Ci         :modes 1-3 read qp specifications from a file (see Remarks)
Ci         :
Ci         :Switches
Ci         :~fn=fnam     specify qp by file 'fnam' (defaults to 'qp')
Ci         :~con         contour plot mode (mode 3)
Ci         :             Bands are written for a uniform mesh of qp
Ci         :             in a parallelepiped defined by (q1,q2) and an offset.
Ci         :             See remarks for specification of (q1,q2).
Ci         :             Fixed row: vary q2 while q1 is fixed
Ci         :             Fixed col: vary q1 while q2 is fixed
Ci         :             When contours are plotted with fplot -con mode,
Ci         :             abscissa corresponds to q1, ordinate to q2.
Ci         :~qp          input file specifies a list of qp (mode 2)
Ci         :~box=#[,#...],shape,q0=#1,#2,#3  (mode 4)
Ci         :             Points generated in a shell around q0; no qp file is read.
Ci         :             # indicates the radius of the shell.
Ci         :             Up to 9 radii may be supplied
Ci         :             'shape' must be one of
Ci         :               n=6[+o]  n=8[+o]  n=12[+o]  n=20[+o]  n=32[+o]
Ci         :             for cube faces, corners, icosahedra vertices or vertices+faces
Ci         :             If +o is appended, the origin q0 is added
Ci         :             q0=#1,#2,#3 specifies q at center of shell
Ci         :~mq          q-points are given as multiples of reciprocal lattice vectors
Ci         :             Applies to symmetry line and qp-list modes only
Ci         :~rot=string  rotate q-point according to rotation specified
Ci         :             by string (syntax defined in a2rotm)
Ci         :~long        write bands with extra digits precision
Ci         :             (has no effect for symmetry line mode)
Ci         :~spin1       generate bands only for first spin
Ci         :~spin2       generate bands only for second spin
Ci         :~nband=#     Write out no more than # bands
Ci         :~lst=list    write only those bands specified in a list.
Ci         :             For syntax of list, see slatsm/mkilst.f
Ci         :~evn=#       keep track of smallest, largest eval for
Ci                       #th band, and printout at close.
Ci         :~ef=#        change efermi to #.
Ci         :~bin         Write into a binary file
Ci         :~col=lst     Assign a list of elements in the eigenvector
Ci                       to make up the first color weight.
Ci         :~col2=lst    Assign a list of elements in the eigenvector
Ci                       to make up a second color weight.
Ci         :~col3=lst    Assign a list of elements in the eigenvector
Ci                       to make up a third color weight.
Ci         :~colst=lst   (noncollinear magnetism only)
Ci                       Assign four color weights ac
Ci         :Example: --band~long~qp~lst=2:5
Ci   iop   :options passed by the caller
Ci         :1s digit suppresses actions for parallel mode
Ci         : mode 1:
Ci         : 1: suppress writing line header info file
Ci         : 2: Return total number of qp to be generated in all lines
Ci         :    and also suppress writing line header info file
Ci         :    Thus iop=2 => all modes return total number of qp to be generated
Ci         : 5: same as as mode 1, but do not open or write anything to bands file
Ci         : 6: same as as mode 2, but do not open or write anything to bands file
Ci         : mode 3:
Ci         : nonzero: do not allocate evsav until iq>nq
Ci         : 7: Make initializations and set arguments (eg nfbn,ifblst,onesp,efermi)
Ci         :    read by parsing sopts. Do not return nq or open a disk file.
Ci         :    Use in conjunction with input nq=0.
Ci         :10s digit
Ci         :1 -> number of bands generated may depend on qp
Ci         :     In that case, the number of bands is stored for each qp.
Ci   nband :(suqlst) maximum number of energy bands to write
Ci         :         If nband==0 => set nfbn=0 and nblst=1; do not set ifblst,iblst
Ci         :(suqlsw) actual number of energy bands to write
Ci
Ci   efermi:Fermi energy (written to bnds file)
Ci
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :NB used only in file write to indicate how many times
Ci         :nsp is ALTERED to nsp=1 if spin1 option is set
Ci
Ci   evl   :eigenvalues for current qp (used only in suqlsw)
Cio Inputs/Outputs (see Remarks)
Cio  nq    :Input nq=0 :if zero, flags suqlst to set up bands mode:
Cio                    :some initializations; then returns ready
Cio                    :for first block of qp
Cio        :           :if nonzero, setup for next qp
Cio        :Output nq  :if zero, suqlst has no more qp to calculate.
Cio                    :if nonzero, nq is the number of k-points in the
Cio                    :current block.  Caller then calls suqlst nq
Cio                    :times (one for each k-point), which will return
Cio                    :qp for each k-point.  After exactly nq calls,
Cio                    :suqlst will start another block, if it exists.
Cio                    :See Remarks for schematic of calling sequence.
Co Outputs
Co  nfbn(1):(color weights) number of elements in ifblst, first color
Co         :Note: sign of nfbn(1) is used to flag spin texture
Co         :nfbn(1) is returned as a negative value if colst used instead of col
Co  nfbn(2):(color weights) number of elements in ifblst, second color
Co  ifblst :(color weights) list of orbital indices for color weights, which is:
Co         :the sum of contributions from eigenvector components in ifblst
Co  qp     :if iq>nq, qp is not set.  Call to suqlst starts a new block (see Remarks).
Co         :Otherwise, qp holds k-point at which to generate bands
Co  onesp  :if spin1 flag is encountered, onesp is set to 1 and nsp is set to 1
Co         :if spin2 flag is encountered, onesp is set to 2 and nsp is set to 1
Cl Local variables
Cl   iq    :current qp in this block
Cl   mode  :1 symmetry-line mode
Cl         :2 list-of-qp mode
Cl         :3 contour mode
Cl         :4 box mode
Cl   q1    :starting qp for symmetry mode; only meaningful in that mode
Cl   q2    :ending   qp for symmetry mode; only meaningful in that mode
Cl   nevn  :(optional) band index; routine monitors largest, smallest value
Cl         :for that index.
Cl   evmnn :smallest value found for qp of specified band index
Cl   evmxn :largest  value found for qp of specified band index
Cl   ifiq  :file logical unit for input qp file
Cl   ifib  :file logical unit for output bands file
Cl   ql    :local copy of current qp
Cl   nblst :0 => number of bands will be q-dependent
Cl         :>0 number of bands to write
Cl   iblst :list of indices to bands to print
Cl         :iblst(1)=0 => list is 1,2,3,...
Cl   qpl   :list of k points
Cl   smapq :ascii string representing algebraic expressions for mapping
Cl         :file q into another q.
Cl         :smapq holds expressions to map (qx,qy,qz) -> kx,ky,kz
Cr Remarks
Cr   suqlst is designed to be called to generate qp in groups or
Cr   blocks.  The calling sequence is:
Cr
Cr     nq = 0  <- flags that first call to suqlst, to set up mode
Cr     do  iblock = 1, forever
Cr      *This call generates nq, the number of points in this block
Cr       call suqlst(s_lat,nband,efermi,nsp,evl,nq,qp,onesp)
Cr       if (nq == 0) stop
Cr       do  iq = 1, nq
Cr        *This call generates qp for current block
Cr         call suqlst(s_lat,ndimh,ef0,nsp,xx,nkp,qp,onesp) <- returns qp
Cr         do  isp = 1, nsp
Cr         call suqlsw(ndimh,qp,evl(1,isp)) <- saves evl for this qp
Cr                                             (call is optional)
Cr         enddo
Cr       enddo
Cr     enddo
Cr
Cr   The following modes are implemented:
Cr     mode=1 reads qp from syml file, and generates qp along each
Cr            specified symmetry line.  Structure of qp file:
Cr            file has one line for each symmetry line as follows:
Cr               nq      q1x   q1y   q1z      q2x   q2y   q2z
Cr               ...
Cr            the line entries have meanings:
Cr            --# qp-   ---starting qp---    --- ending qp ---
Cr            Any line with nq=0 implies no more lines.
Cr
Cr     mode=2 reads qp from specified file and generates qp for each
Cr            specified qp.  File consists sets of qpx,qpy,qpz for each
Cr            qp sought.  Files can be in the 'standard' format defined
Cr            in getqp.f, which has this structure:
Cr
Cr              nkp=#  (some additional information may be here)
Cr               1  q1x  q1y  q1z  [w1]
Cr               2  q2x   q2y   q2z
Cr               ...
Cr            Alternatively you can supply information in a generic format.
Cr            It typically consists of lines like this:
Cr               q1x   q1y   q1z
Cr               q2x   q2y   q2z
Cr               ...
Cr            suqlst uses rdm to read the qp file in this format, which
Cr            permits algebraic expressions and (optional) header info.
Cr
Cr     mode=3 generates qp for a uniform mesh in a plane (contour plot)
Cr            The file supplies information describing a rectangle in
Cr            the Brillouin zone.  It consists of a single line,
Cr            which contains the following:
Cr             v1    range  n     v2    range  n   height  list-of-bands
Cr
Cr            v1 and v2 are two vectors specifying the plane of the
Cr            contour.  range and n (one each for v1 and v2) are the
Cr            starting and final amplitudes of those vectors, and the
Cr            the number of points within the vector.  list-of-bands
Cr            is a list of integers which specify which bands are to
Cr            be written to the output file.  'height' is the 'z' axis.
Cr            For example,
Cr             v1    range  n     v2    range  n   height  list-of-bands
Cr            1 0 0  -1 1   51   0 1 0  -1 1   51   0.00    4,5
Cr            creates a file of 51x51 points, with the four corners
Cr            (-1,-1,0),  (1,-1,0),  (-1,1,0),  (1,1,0)
Cr            going through the gamma-point. Two bands (4,5) are stored.
Cr
Cr   Color weights.  The weights are determined by the Mulliken
Cr   decomposition of the norm.
Cr   Orthogonal basis : D_in = (z_in)+ z_in where
Cr     i = orbital index and n = band index
Cr   Nonorthogonal basis : D_in = (z~_in)+ z_in
Cr     Here z~ is contravariant form of z.
Cr     Overlap matrix is S = (z z+)^-1
Cr     z~+ = z+ S = z+ (z+)^-1 z^-1 = z^-1
Cr   Algebraic transformation of q.
Cr     At times, q actually sought may differ from nominal q,
Cr     e.g. in a photoemission experiment, k|| surface is modified
Cr     You can supply algrebraic expressions involving q,qx,qy,qz
Cr     that generate a transformed q.  This is done by inserting
Cr     a line at the top of the q-points input file such as:
Cr        % kx=qx+.1 ky=qx+qy
Cr     The first character must be a '%', followed by one or more strings
Cr        kx=expr1  ky=expr2  kz=expr3
Cr     Any of the x,y,z components of q for which a string is supplied
Cr     is transformed.
Cu Updates
Cu   06 Jun 17 Color weights for spin texture
Cu   13 Nov 16 New iop=5,6
Cu   03 Oct 15 New 'box' mode
Cu   27 Mar 15 New mq switch
Cu   02 Jul 12 Option to transform q by rotation or expr in q-points file
Cu   30 Mar 11 Bug fix, qp list mode, binary save with list
Cu   01 Jul 10 Bug fix, head info qp list mode, spin polarized case
Cu   17 Feb 10 New switch modifiers to 'qp' --- can do
Cu             qp[,inc=expr][,merge=expr][,save=expr]
Cu   23 Jan 10 Substitute newer mkilsd for mkilst
Cu   12 Aug 09 qp mode can read io file format as defined in getqp.f
Cu             better treatment of qp-dependent number of bands
Cu   08 Jul 08 Extend to case where number of bands can be q dependent
Cu             modes 1,2: suqlsw writes out number of bands with qp
Cu   09 Jul 07 configured to with MPIK mode
Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1
Cu   02 Jul 06 Color mode extended to two colors
Cu   02 Jul 06 New color mode (one color weight only)
Cu   14 Feb 05 contour mode saves both spins in spin-polarized case
Cu   20 Oct 03 suqlst works properly in contour mode
Cu   28 Aug 01 added ef switch
Cu   23 Jan 01 first written
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Passed parameters
      character*(*) sopts
      integer nq,nband,nsp,nfbn(3),ifblst(nband,3),onesp,iop
      double precision efermi,qp(3),evl(nband)
C ... Dynamically allocated arrays
      real(8), pointer :: qpl(:,:),qptmp(:)
C ... Local variables
      character strn*512, strn2*512, dc*1, fn*120, prfmt*40, rots*120
      character qcnst*80, savef*80, mergef*80, smapq*160
      logical:: fmlong,lbin,lmapqp,lmapq(3)=.false.,lqlat=.false.
      integer i,ifib,ifiq,ip,iq,j,j1,j2,jsp,k,mode,nblst,nbox,
     .  nbsave,nevn,nq2,nqall,nqx,nqy,nrad,nspx,op1,op2,stdo,iqx,iqy
      integer iblst(200),iv(12)
      integer,parameter :: NULLI =-99999
      real(8),parameter :: pi = 4d0*datan(1d0)
      double precision xx,xx1,xx2,sclp,plndst,xxv(12),dlength
      double precision x1,x2,y1,y2,evmxn,evmnn,eferm,rotm(3,3),
     .  q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3),qm1(3),qm2(3),box(100)
      real(8),allocatable :: evsav(:,:,:,:),evwk(:,:)
      procedure(logical) :: rdstrn
      procedure(integer) :: iprint,nglob,parg,a2vec,fopna,fopnn,fopno,mkilsd
      procedure(real(8)) :: angle
C ... MPI
      integer procid,master,mpipid
      common /suqlsd/
     .  x1,x2,y1,y2,vectx,vecty,vect3,
     .  q1,q2,ql,qpl,evmxn,evmnn,ifiq,ifib,mode,iq,nevn,nblst,
     .  iblst,fmlong,nqx,nqy
      save evsav,lbin,nbsave,nspx,eferm,smapq,rots
      data lbin /.false./ rots / ' '/
C ... External calls
      external a2rotm,awrit2,awrit3,awrit4,awrit5,clrsyv,cross,dcopy,dfclos,
     .  dgemm,dpscop,dscal,fclose,fexit2,fpiint,getqp,getsyv,ilst2a,info0,info2,
     .  info5,isanrg,ivset,lodsyv,mkils0,mkilss,mkilst,numsyv,nwordg,rx0,rx2,
     .  rxi,rxs,rxx,strip,suqlsc,suqlsf,suqlsz,words,ywrm,zcopy,zgetrf,zgetri

      procid = mpipid(1)
      master = 0
      stdo = nglob('stdo')
      op1 = mod(mod(iop,10),4)
      op2 = mod(iop/10,10)
      nqall = 0

C --- First call ... setup and parse options ---
      if (nq == 0) then

C   ... Defaults
        fmlong = .false.
        mode = 1
        fn = 'qp'
        qcnst = ' '; savef = ' '; mergef = ' '; smapq = ' '
        nblst = 0; if (op2 == 0) nblst = max(nband,1)  ! fixed number of bands
        nevn = 0; evmxn = -99d9; evmnn = 99d9
        nfbn = 0; if (nband>0) iblst(1) = 0
        if (len(sopts) == 0) goto 15
        dc = sopts(1:1)
        if (dc /= ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 0
   10   continue
        j2 = j2+1
        if (j2 > len(sopts)) goto 15
        if (sopts(j2:j2) == dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 < j1) goto 15
        if (.false.) then
        elseif (sopts(j1:j1+2) == 'fn=')  then
          if (j1+3 > j2) call rx('suqlst: bad file name')
          fn = sopts(j1+3:j2)
        elseif (sopts(j1:j1+3) == 'box=') then
          if (j1+4 > j2) goto 996
          mode = 4; nbox=NULLI; q1(1) = NULLI
          i = j1+3
C          j = a2vec(sopts,j2,i,4,',~',2,1,1,iv,box)
C          if (j /= 1) goto 996
          i = j1+3
          j = a2vec(sopts,j2,i,4,',~',2,2,10,iv,box)
          if (j == -1) goto 996
          if (j > 9) goto 996
          nrad = -j-1
          j1 = i+1
   12     continue
C          if (j1 > j2 .and. q1(1) == NULLI) goto 996 ! q0 must be assigned before exit
C          if (j1 > j2 .and. nbox  == NULLI) goto 996 ! nbox must also be assigned
          if (j1 <= j2)  then
            k = index(sopts(j1:j2),',') + j1-2
            if (k < j1) k = j2
            if (.false.) then
            elseif (sopts(j1:j1+2) == 'q0=') then
              i = j1+2
              j = a2vec(sopts,j2,i,4,', ',2,3,3,iv,q1)
              if (j /= 3) goto 996
              k = i
C            elseif (sopts(j1:j1+1) == 'n=') then
C              i = j1+1
C              j = a2vec(sopts,j2,i,2,', ',2,1,2,iv,nbox)
C              if (j /= 1) goto 996
C              k = i
            elseif (sopts(j1:j1+4) == 'n=6+o') then
              nbox = -6
              k = j1+5
            elseif (sopts(j1:j1+4) == 'n=8+o') then
              nbox = -8
              k = j1+5
            elseif (sopts(j1:j1+5) == 'n=12+o') then
              nbox = -12
              k = j1+6
            elseif (sopts(j1:j1+5) == 'n=20+o') then
              nbox = -20
              k = j1+6
            elseif (sopts(j1:j1+5) == 'n=32+o') then
              nbox = -32
              k = j1+6
            elseif (sopts(j1:j1+2) == 'n=6') then
              nbox = 6
              k = j1+3
            elseif (sopts(j1:j1+2) == 'n=8') then
              nbox = 8
              k = j1+3
            elseif (sopts(j1:j1+3) == 'n=12') then
              nbox = 12
              k = j1+4
            elseif (sopts(j1:j1+3) == 'n=20') then
              nbox = 20
              k = j1+4
            elseif (sopts(j1:j1+3) == 'n=32') then
              nbox = 32
              k = j1+4
C            elseif (sopts(j1:j1+5) == 'face+o') then
C              nbox = -6
C              k = j1+6
C            elseif (sopts(j1:j1+3) == 'face') then
C              nbox = 6
C              k = j1+4
C            elseif (sopts(j1:j1+5) == 'corn+o') then
C              nbox = -8
C              k = j1+6
C            elseif (sopts(j1:j1+3) == 'corn') then
C              nbox = 8
C              k = j1+4
C            elseif (sopts(j1:j1+5) == 'icos+o') then
C              nbox = -12
C              k = j1+6
C            elseif (sopts(j1:j1+3) == 'icos') then
C              nbox = 12
C              k = j1+4
            else
              goto 996
            endif
            j1 = k+1
            goto 12
          endif
        elseif (sopts(j1:j1+2) == 'qp,')  then
          mode = 2
          j1 = j1+3
   11     continue
          if (j1 <= j2)  then
            k = index(sopts(j1:j2),',') + j1-2
            if (k < j1) k = j2
            if (.false.) then
            elseif (sopts(j1:j1+3) == 'inc=') then
              if (j1+4 > k) goto 997
              qcnst = sopts(j1+4:k)
            elseif (sopts(j1:j1+4) == 'save=') then
              if (j1+5 > k) goto 997
              savef = sopts(j1+5:k)
            elseif (sopts(j1:j1+5) == 'merge=') then
              if (j1+6 > k) goto 997
              mergef = sopts(j1+6:k)
            elseif (sopts(j1:j1) == ',') then
              j1 = j1+1
              goto 11
            else
              goto 997
            endif
            j1 = k+1
            goto 11
          endif
        elseif (sopts(j1:j2) == 'qp')  then
          mode = 2
        elseif (sopts(j1:j2) == 'spin1')  then
          onesp = 1
          nsp = 1
        elseif (sopts(j1:j2) == 'spin2')  then
          onesp = 2
          nsp = 1
        elseif (sopts(j1:j2) == 'con')  then
          mode = 3
        elseif (sopts(j1:j2) == 'bin')  then
          lbin = .true.
        elseif (sopts(j1:j2) == 'mq')  then
          lqlat = .true.
        elseif (sopts(j1:j2) == 'long')  then
          fmlong = .true.
        elseif (sopts(j1:j1+3) == 'rot=')  then
          rots = sopts(j1+4:j2)
        elseif (sopts(j1:j1+3) == 'col=' .or. sopts(j1:j1+5) == 'colst=')  then
          if (nband<=0) goto 10
          i = 4; if (sopts(j1:j1+5) == 'colst=')  i = 6
          if (j1+i > j2) call rx('suqlst: bad list, col=..')
          nfbn(1) = mkilsd(sopts(j1+i:j2),-1,ifblst)
          if (nfbn(1) < 0) call rx('suqlst: bad list, col=..')
          if (nfbn(1) > nband) call rxi('suqlst: orbital list, col=.. exceeds nband = ',nband)
          nfbn(1) = mkilsd(sopts(j1+i:j2),nband,ifblst)
          if (i == 6) nfbn(1) = - nfbn(1)
        elseif (sopts(j1:j1+4) == 'colst') then
          if (nband<=0) goto 10
          nfbn(1) = -nband/2
          nfbn(2) = 0
          forall (i = 1:nband/2) ifblst(i,1) = i
        elseif (sopts(j1:j1+4) == 'col2=') then
          if (nband<=0) goto 10
          if (j1+4 > j2) call rx('suqlst: bad list, col2=..')
          if (nfbn(1) < 0) call rx('suqlst: cannot use col2 with colst')
          nfbn(2) = mkilsd(sopts(j1+5:j2),-1,ifblst(1,2))
          if (nfbn(2) < 0) call rx('suqlst: bad list, col=..')
          if (nfbn(2) > nband) call rxi
     .      ('suqlst: orbital list, col2=.., exceeds nband = ',nband)
          nfbn(2) = mkilsd(sopts(j1+5:j2),nband,ifblst(1,2))
        elseif (sopts(j1:j1+4) == 'col3=') then
          if (nband<=0) goto 10
          if (j1+4 > j2) call rx('suqlst: bad list, col3=..')
          if (nfbn(2) <= 0) call rx('suqlst: cannot use col3 without col2')
          nfbn(3) = mkilsd(sopts(j1+5:j2),-1,ifblst(1,3))
          if (nfbn(3) < 0) call rx('suqlst: bad list, col=..')
          if (nfbn(3) > nband) call rxi
     .      ('suqlst: orbital list, col3=.., exceeds nband = ',nband)
          nfbn(3) = mkilsd(sopts(j1+5:j2),nband,ifblst(1,3))
        elseif (sopts(j1:j1+5) == 'nband=') then
          if (nband<=0) goto 10
          j = 0
          i = parg('nband=',2,sopts(j1:),j,len(sopts(j1:)),dc//' ',1,1,i,nblst)
          if (i <= 0) call rxs('suqlst: failed to parse string for number of bands ',sopts(j1:))
        elseif (sopts(j1:j1+3) == 'lst=')  then
          if (nband<=0) goto 10
          if (j1+4 > j2) call rx('suqlst: bad list, lst=..')
          call mkils0(sopts(j1+4:j2),nblst,iblst)
          if (nblst > 200) call rx('suqlst: increase size of iblst')
          call mkilst(sopts(j1+4:j2),nblst,iblst)
        elseif (sopts(j1:j1+2) == 'ef=')  then
          j = 0
          i = parg('ef=',4,sopts(j1:),j,len(sopts(j1:)),dc//' ',1,1,i,efermi)
          if (i <= 0) call rxs('suqlst: failed to parse string for fermi level:  ',sopts(j1:))
          call info2(20,0,0,' suqlst:  user specified Fermi level to be Ef=%,5;5d',efermi,0)
        elseif (sopts(j1:j1+3) == 'evn=')  then
          if (j1+4 > j2) call rx('suqlst: bad list')
          i = j1+3; xxv(1) = 0
          j = a2vec(sopts,j2,i,2,dc//' ',2,3,1,iv,nevn)
          if (j /= 1 .or. nevn > nband) call rx('suqlst: bad value for evn')
        else
          call rxs('suqlst: failed to parse argument, ',sopts(j1:j2))
        endif
        goto 10
        endif
   15   continue

        if (mode == 4) then
          call info2(20,0,0,' suqlst: generate bands in cluster around qp=%3;11,6D (box mode)',q1,0)
        elseif (nfbn(1) == 0) then
          call info2(20,0,0,' suqlst: read qp for bands, mode %i',mode,0)
        elseif (nfbn(2) == 0) then
          call ilst2a(ifblst,iabs(nfbn(1)),strn)
          call info5(20,0,0,
     .      ' suqlst:  generate bands with %?#(n<0)#spin texture#color# weights, mode %i.'//
     .      ' %N%10f%i components: '//strn//'%a',nfbn(1),mode,iabs(nfbn(1)),4,5)
        else
          j = 2 ; if (nfbn(3) /= 0) j = 3
          call info2(20,0,0,' suqlst:  generate bands with %i color weights, mode %i.',j,mode)
          do  i = 1, j
            call ilst2a(ifblst(1,i),nfbn(i),strn)
            call info2(20,0,0,'%10f%i components, color %i:  '//strn//'%a',nfbn(i),i)
          enddo
        endif

        if (mod(iop,10) == 7) return

        if (procid == master) then
C     ... Open qp file
          if (mode /= 4) then
            ifiq = fopno(fn)
            rewind ifiq
          endif

          if (mod(iop,10) > 4) goto 25 ! Do not open or write anything to bnds file

C     ... Open bands file
          if (lbin) then
            if (mode == 1) then
              call rx('suqlst: binary write not implemented with mode 1')
            endif
            if (mode == 2 .or. mode == 4) then
              ifib = fopna('tmp2',-1,4+2)
            else
              ifib = fopna('bbnds',-1,4+2)
            endif
          else
            ifib = fopnn('BNDS')
          endif
          rewind ifib

C     ... Write header
          eferm = efermi
          if (mode == 1) then
            i = nblst
            if (nblst == 0) i = nband
            if (nfbn(1) == 0) then
C           Use separate format statment to circumvent gfortran bug
              write(ifib,335) i,efermi,0
  335         format(i5,f10.5,i6)
            elseif (nfbn(2) == 0) then
              call ilst2a(ifblst,iabs(nfbn(1)),strn)
              call strip(strn,j1,j2)
              if (nfbn(1) < 0) then
                write(ifib,338) i,efermi,4,strn(j1:j2)
  338           format(i5,f10.5,i6,'  spin-texture col= ',a)
              else
                write(ifib,336) i,efermi,1,strn(j1:j2)
              endif
  336         format(i5,f10.5,i6:'  col= ',a:'  col2= ',a)
            else
              j = 2 ; if (nfbn(3) /= 0) j = 3
              write(ifib,336,advance='no') i,efermi,j
              do  i = 1, j
                call ilst2a(ifblst(1,i),nfbn(i),strn2)
                call awrit1('%x col%i='//strn2,strn,len(strn),0,i)
                write(ifib,'(a)',advance='no') trim(strn)
              enddo
              write(ifib,'("")')
            endif
          endif

   25     continue
        endif

C   ... Other initializations
        iq = 0
        nq = -1
      endif

C --- Setup for a new block of k-points, depending on mode ---
C     This branch occurs on completion of the last qp of the current block
C     which is marked by iq>nq.
C     First pass: iq = 0 and nq = -1 => iq>nq
C     At the completion of this block: nq must be computed and:
C     (mode=1) q1,q2 set up.
C              Note: if 1s digit of iop is set in this mode,
C              this branch returns sum of all qp in all lines.
C              No setup for q1,q2; no
C     (mode=2) qpl allocated and loaded
C     (mode=3) nblst = number of bands to save
C     (mode=4) qpl allocated and loaded
      if (iq > nq) then
        iq = 1

C   ... Bands along specified symmetry lines
        if (mode == 1) then
  725     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 999
          strn2(1:2) = adjustl(strn)
          if (strn2(1:1) == '#' .or. strn2(1:1) == ' ') goto 725
          if (strn(1:1) == '%') then
            j1 = 2
            call nwordg(strn,1,' ',1,j1,j2)
            if (strn(j1:j2) /= 'map') goto 725
            do  i = 1, 3
              lmapq(i)=index(strn,' k'//char(ichar('x')+i-1)//'=')>0
            enddo
            smapq = strn(j2+1:)
            if (lmapq(1) .or. lmapq(2) .or. lmapq(3)) then
              call info0(20,0,0,' suqlst: map qp using the following:')
              call info0(20,0,0,trim(smapq))
            endif
            goto 725
          endif

C         Exit if first entry in line is zero
          i = 0
          xxv(1) = 0
          i = a2vec(strn,len(strn),i,4,', ',2,3,1,iv,xxv)
          if (i == 1 .and. xxv(1) == 0) goto 999
          i = 0
          i = a2vec(strn,len(strn),i,4,', ',2,3,7,iv,xxv)
          if (i /= 7 .and. iprint()>=10) then
            write(stdo,'(/'' suqlst (warning) skipping line:''/''  '',a)') strn
            goto 725
          endif
          nq = xxv(1)
C         1 qp is nonsensical for a line
          if (nq == 1) nq = 0
C         No qp: exit
          if (nq <= 0) goto 999
C         setup q1,q2
          if (lqlat) then ! Convert q1,q2 (multiples of qlat) to Cartesian coordinates
            call dgemm('N','T',1,3,3,1d0,xxv(2),1,s_lat%qlat,3,0d0,q1,1)
            call dgemm('N','T',1,3,3,1d0,xxv(5),1,s_lat%qlat,3,0d0,q2,1)
          else
            call dcopy(3,xxv(2),1,q1,1)
            call dcopy(3,xxv(5),1,q2,1)
          endif
          if (iprint()>=20) write(stdo,"(1x)")
          if (iprint()>10) write(stdo,785) nq,q1,q2
  785     format(' suqlst:  nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
          if (rots /= ' ') then
            call info0(20,0,0,' Shift q by rotation : '//trim(rots))
          endif
          if (lmapqp(smapq,lmapq,q1,qm1) .or.
     .        lmapqp(smapq,lmapq,q2,qm2)) then
            if (iprint()>=10) write(stdo,
     .        "(8x,' mapped to q1=',3f7.4,'   q2=',3f7.4)") qm1,qm2
          endif
C         Write line header information to disk
          if (op1 == 0) then
            write(ifib,337) nq*nsp
  337       format(2i5)
C         Accumulate all qp and cycle until all lines are read
          elseif (op1 == 2) then
            nqall = nqall + nq
            goto 725
          endif

C   ... (mode==2) Bands for a list of qp specified in file ifiq
C       (mode==4) Bands in a sphere around q1
        elseif (mode == 2 .or. mode == 4) then
C         Only one block for this mode.
C         Flag that prior block already completed: nq>0
          if (nq > 0 .and. lbin) then
            rewind ifib
            read(ifib) i,j
            if (i/=nq*nspx) call rx('suqlst: tmp file was corrupted')
            allocate(evsav(nq*nspx,j,1,1))
            do  i = 1, nq*nspx
              read(ifib) evsav(i,:,1,1)
            enddo
C           Close and remove temporary file tmp2
            call dfclos(ifib)
            ifib = fopna('bbnds',-1,4+2)
            rewind ifib
            call ywrm(1,' ',1,ifib,' ',evsav,0,nq*nspx,nq*nspx,j)
            deallocate(evsav)
            call fclose(ifib)
          endif
          if (nq > 0) goto 999

C     ... Read qp from file in one go
          if (mode == 2) then
            nq = 0
            nq2 = 0
            call suqlsf(0,ifiq,0,nq,iv,xxv)
            if (nq < 0) goto 998
C           Case new file is to be merged with qp list
            if (mergef /= ' ') then
              call suqlsf(0,fopno(mergef),0,nq2,iv,xxv)
              if (nq2 < 0) then
                fn = mergef
                goto 998
              endif
            endif
            allocate(qpl(3,(nq+nq2)))
C           Read qp from first qp file
            call suqlsf(1,ifiq,0,nq,iv,qpl)
            if (nq < 0) goto 998
C           If constraint given, reduce qp to those meeting constraint
            if (qcnst /= ' ') then
              allocate(evwk(3,nq))
              call dcopy(nq*3,qpl,1,evwk,1)
              call suqlsc(iprint(),qcnst,nq,evwk,nq,qpl)
              deallocate(evwk)
              call ivset(iv,1,7,0)
            endif

C           Case new file is to be merged with qp list
            if (mergef /= ' ') then
              call suqlsf(1,fopno(mergef),nq,nq+nq2,iv,qpl)
              call ivset(iv,1,7,0)
              nq = nq+nq2
            endif
C           Case qp list is to be saved on disk: save and quit.
            if (savef /= ' ') then
              ifiq = fopnn(savef)
              rewind ifiq
              call getqp(1,-ifiq,nq,iv,iv(4),iv(7),qpl,xx,xx)
              call rx0('wrote qp list to file: '//trim(savef))
            endif

C           Case qp are in multiples of qlat
            if (lqlat) then     ! Convert qp (multiples of qlat) to Cartesian coordinates
              allocate(qptmp(3*nq))
              call dcopy(3*nq,qpl,1,qptmp,1)
              call dgemm('N','N',3,nq,3,1d0,s_lat%qlat,3,qptmp,3,0d0,qpl,3)
              deallocate(qptmp)
            endif

C         mode 4 : qp in a sphere around q1
          elseif (mode == 4) then
            if (q1(1) == NULLI) goto 996 ! q0 must be assigned
            if (nbox  == NULLI) goto 996 ! nbox must also be assigned
            nq = iabs(nbox)
            if (nq/=6 .and. nq/=8 .and. nq/=12 .and. nq/=20 .and. nq/=32)
     .      call rx('suqlst: --box requires 6, 8, 12, 20, or 32 qp')
            nq = iabs(nbox)*nrad
            if (nbox < 0) nq = nq+1
            allocate(qpl(3,nq),evwk(3,iabs(nbox)))
            call fpiint(-iabs(nbox),0,i,evwk,qpl) ! Store them into evwk
            qpl(:,1) = q1
            k = 0
            if (nbox < 0) k = k+1
            do  j = 1, nrad
              do  i = 1, iabs(nbox)
                k = k+1
                qpl(1:3,k) = q1(1:3) + box(j)*evwk(1:3,i)
              enddo
            enddo
            mode = 2            ! Subsequent flow identical to mode=2
            fmlong = .true.     ! Always write many digits
            if (k /= nq) call rx('bug in suqlst')
          endif

          i = nblst
          if (nblst == 0) i = nband
          nspx = nsp
          if (onesp /= 0) nspx = 1
          call info5(30,0,0,' generating %i band%-1j%?#(n>1)#s## at %i qp'//
     .      '%?#(n==2)#, %i spin(s)##',i,nq,nglob('nsp'),nspx,0)
          if (i == 0) call rx('suqlst mode 4 must have nband>0')

          if (mod(iop,10) <= 4) then
            if (nfbn(1) /= 0) i = i*2
            nbsave = i
            if (lbin) then
              write(ifib) nq*nspx,i+3,1
            elseif (nblst > 0) then
              call awrit4('%% rows %i cols %i  efermi=%;6d  nsp=%i',' ',
     .          80,ifib,nq*nspx,i+3,efermi,nspx)
            else
              call awrit4('%% rows %i cols %i  nlb  efermi=%;6d  nsp=%i',
     .          ' ',80,ifib,nq*nspx,i+4,efermi,nspx)
            endif
          endif

C       Bands on a uniform mesh in a specified plane (contour plot)
        elseif (mode == 3 .and. nq == -1) then
  825     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 998
          if (strn(1:1) == '#') goto 825
          call words(strn,i)
          if (i > 14) then
            call word(strn,15,j1,j2)
            if (strn(j1:j1) /= '#') then
              call rxi('suqlst con mode: expected 14 arguments from input file but read',i)
            endif
          endif
C         Parse argument 14: list of bands
          iblst(1) = -1
          nblst = 1
          call word(strn,14,j1,j2)
          call mkilss(11,strn(j1:j2),nblst,iblst)
          if (nblst <= 0) call rx('suqlst: no bands in list')
C         Parse the first 12 arguments
          j2 = j1-1
          j1 = 1
          ip = 0
          ip = a2vec(strn,len(strn),ip,4,' ',1,-2,-12,iv,xxv)
          if (ip /= 12) call rx(' suqlst: failed to parse '//strn(j1:j2))

          if (lqlat) then ! Convert qstart,qend (multiples of qlat) to Cartesian coordinates
            call dgemm('N','T',1,3,3,1d0,xxv(1),1,s_lat%qlat,3,0d0,vectx,1)
            call dgemm('N','T',1,3,3,1d0,xxv(7),1,s_lat%qlat,3,0d0,vecty,1)
          else
            call dcopy(3,xxv(1),1,vectx,1)
            call dcopy(3,xxv(7),1,vecty,1)
          endif
          x1 = xxv(4); x2 = xxv(5); nqx = nint(xxv(6))
          y1 = xxv(10);y2 = xxv(11);nqy = nint(xxv(12))

C         Fortran read:
C         backspace ifiq; read(ifiq,*) vectx,x1,x2,nqx,vecty,y1,y2,nqy
          nq = nqx*nqy

C         Parse argument 13: height (1 expression) or origin (3 expressions)
          call word(strn,13,j1,j2)
          ip = 0
          ip = a2vec(strn(j1:j2),j2-j1+1,ip,4,', ',2,-3,-3,iv,vect3)
          if (ip /= 1 .and. ip /= 3) call rx(' suqlst: string "'//strn(j1:j2)//
     .      '" in file '//trim(fn)//' must consist of 1 or three expressions')

          if (ip == 1) then
            plndst = vect3(1)
C           01 Mar 2011 No longer orthonormalize vectx, vecty
C           call dscal(3,1/dsqrt(sclp(vectx,vectx)),vectx,1)
C           call dscal(3,1/dsqrt(sclp(vecty,vecty)),vecty,1)
C           Subtract from vecty projection onto vectx
C           call daxpy(3,-sclp(vectx,vecty),vectx,1,vecty,1)
            call cross(vectx,vecty,vect3)
            call dscal(3,plndst/dsqrt(sclp(vect3,vect3)),vect3,1)
            if (iprint() >= 10) write(stdo,716) plndst, vect3
  716       format(' Using h=',f9.6,'  origin =',3f10.6)
C          elseif (ip /= 3) then
          endif
          if (iprint() >= 10) then
            write(stdo,717) vectx,x1,x2,nqx,vect3,vecty,y1,y2,nqy,nq
  717       format(' qx=',3f9.6,'  x1,x2=',2f9.6,'  nx=',i3,'   o=',3f10.6
     .            /' qy=',3f9.6,'  y1,y2=',2f9.6,'  ny=',i3,'  np=',i5)
            do  k = 1, 3
              xxv(k)   = x1*vectx(k) + y1*vecty(k) + vect3(k)
              xxv(3+k) = x1*vectx(k) + y2*vecty(k) + vect3(k)
            enddo
            write(stdo,719) 1,1,(xxv(k), k=1,3), 1,nqy,(xxv(k), k=4,6)
  719       format(' q(',i4,',',i4,')',3f10.6,'  q(',i4,',',i4,')',3f10.6)
            do  k = 1, 3
              xxv(6+k) = x2*vectx(k) + y1*vecty(k) + vect3(k)
              xxv(9+k) = x2*vectx(k) + y2*vecty(k) + vect3(k)
            enddo
            xx1 = dlength(3,xxv(7:9)-xxv(1:3),1)
            xx2 = dlength(3,xxv(4:6)-xxv(1:3),1)
            write(stdo,719) nqx,1,(xxv(k), k=7,9),nqx,nqy,(xxv(k),k=10,12)
            write(stdo,720) xx1,xx2
  720       format(8x,'width',f10.6,4x,'height',f10.6)
            write(strn,'('' save %i bands: %'',i2,'':1i'')') nblst
            if (strn(17:18) == '% ') strn(17:18) = ' %'
            call awrit2(strn,strn,80,stdo,nblst,iblst)
            xx1 = angle(3,vectx,vecty)
            if (dabs(xx1) > 1d-6) call info2(20,0,0,
     .        ' suqlst (warning): input plane vectors not orthogonal: angle %;3d deg',
     .        180/pi*acos(xx1),2)
          endif
          if (op1 == 0) then
          allocate(evsav(nqx,nqy,nblst,nsp))
          endif

C       Contour plot, cleanup.
C       Note: cleanup handled by suqlsw when last qp is called.
        elseif (mode == 3) then
C          if (op1 /= 0) then
C            allocate(evsav(nqx,nqy,nblst,nsp))
C            return
C          endif
          call rx('suqlst: caller should never reach this branch')
        endif

C --- Generate qp for this iq, depending on mode ---
      else
        if (mode == 1) then
          xx = dble(iq-1)/dble(nq-1)
          qp(1) = xx*q2(1) + (1-xx)*q1(1)
          qp(2) = xx*q2(2) + (1-xx)*q1(2)
          qp(3) = xx*q2(3) + (1-xx)*q1(3)
        elseif (mode == 2) then
          call dpscop(qpl,qp,3,iq*3-2,1,1d0)
        elseif (mode == 3) then
C         Inner Loop:  excursions in y; outer loop: excursions in x
          j = mod(iq-1,nqy)
          i = (iq-1-j)/nqy
          if (nqx <= 1) then
            xx1 = x1
          else
            xx1 =i*(x2-x1)/(nqx-1)+x1
          endif
          if (nqy <= 1) then
            xx2 = y1
          else
            xx2 =j*(y2-y1)/(nqy-1)+y1
          endif
          do  k = 1, 3
            qp(k) = xx1*vectx(k) + xx2*vecty(k) + vect3(k)
          enddo
          if (j == 0 .and. iprint()>=20) write(stdo,718) i+1,nqx,qp
  718     format(' line',i3,' of',i3,'   q(1)=',3f10.6)
        else
          call rx('suqlst: bad mode')
        endif

C   ... Rotate or otherwise transform q
        if (rots /= ' ') then
          call a2rotm(rots,.false.,0,rotm)
C         call a2rotm(rots,.true.,0,rotm)
C         In-line multiply avoids bug in DEC fort compiler
          do  i = 1, 3
            qm1(i) = rotm(i,1)*qp(1) + rotm(i,2)*qp(2) + rotm(i,3)*qp(3)
          enddo
          call dcopy(3,qm1,1,qp,1)
        endif
        if (lmapqp(smapq,lmapq,qp,qp)) then
        endif
        iq = iq+1
C       Hold onto local copy of qp
        call dcopy(3,qp,1,ql,1)
      endif
      return

C --- No more qp blocks: cleanup ---
  999 continue
      nq = 0
      if (nevn /= 0) call awrit3(' eval no. %i:  minimum eval'//
     .  ' = %;8F  maximum eval = %;8F',' ',80,stdo,nevn,evmnn,evmxn)

      if (mode == 1) then
        if (op1 == 0) then
          write(ifib,337) 0
        elseif (op1 == 2) then
          nq = nqall
        endif
      elseif (mode == 3) then
        call rx('suqlst mode==3 not ready')
C       call xxxbnd(w(oev),nblst,nqx,nqy,ifib)
C       return
      endif
      return

C --- Error exit ---
  998 call rxs('suqlst: failed to read file contents, file ',fn)
  997 call rxs('suqlst: failed to parse options modifying "qp" : ',sopts(j1:j2))
  996 call info0(1,0,0,' suqlst: box syntax is:  box=#[,#..],shape,q0=#,#,#'//
     .  '%N%9fwith shape one of n=6[+o]  n=8[+o]  n=12[+o]  n=20[+o]  n=32+o')
      call rx('failed to parse arguments to box mode')

      entry suqlsm(iop)
C- Return qlist mode
      iop = mode
      return

      entry suqlsxy(iqx,iqy)
C- Return nqx, nqy (2D contour mode)
      iqx = nqx
      iqy = nqy
      return

      entry suqlsw(nband,jsp,nsp,evl)
C- Write or store the energy bands to file for this qp

C ... line mode
      if (mode == 1) then
        if (nblst == 0) then
          prfmt = '(3f10.5,i6/(10f8.4))'
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        elseif (nblst > 0 .and. iblst(1) == 0) then
          prfmt = '(3f10.5/(10f8.4))'
          if (nband < nblst) then
            call info2(30,0,0,' suqlsw (warning) writing %i bands '//
     .        'but generated only %i',nblst,nband)
            allocate(evwk(nblst,1))
            evwk = evl(nband)
            evwk(1:nband,1) = evl(1:nband)
            write(ifib,prfmt) ql, (evwk(i,1),i=1,nblst)
            deallocate(evwk)
          else
            write(ifib,prfmt) ql, (evl(i),i=1,nblst)
          endif
        elseif (nblst > 0) then
          prfmt = '(3f10.5/(10f8.4))'
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        else
          call rx('bug in suqlsw')
        endif

C ... qp list mode
      elseif (mode == 2) then
        if (fmlong .and. nblst == 0) then
          prfmt = '(3f15.10,i6/(5f15.10))'
        elseif (nblst == 0) then
          prfmt = '(3f10.6,i6/(8f10.6))'
        elseif (fmlong .and. nblst > 0) then
          prfmt = '(3f15.10/(5f15.10))'
        elseif (nblst > 0) then
          prfmt = '(3f10.6/(8f10.6))'
C         if (nblst <= 5) prfmt = '(8f10.6)'
        else
          call rx('bug in suqlsw')
        endif
        if (lbin) then
          allocate(evwk(nbsave+3,1))
          evwk = 9999d0
          k = min(nband,nbsave)
          evwk(1:3,1) = ql
          if (nblst == 0 .or. iblst(1) == 0) then
            evwk(3+1:3+k,1) = evl(1:k)
          else
            do  i = 1, k
              evwk(3+i,1) = evl(iblst(i))
            enddo
          endif
          write(ifib) evwk
          deallocate(evwk)
        else


          if (nblst == 0) then ! General case, q-dependent nband
            write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
            if (nbsave > nband) then
              do  i = nband+1, nbsave-1
                write(ifib,"(' 9999')",advance='no')
              enddo
              write(ifib,"(' 9999')")
            endif
          elseif (nblst > 0 .and. iblst(1) == 0) then !Fixed nblst
            if (nband < nblst) then
              call info2(30,0,0,' suqlsw (warning) writing %i bands '//
     .          'but generated only %i',nblst,nband)
              allocate(evwk(nblst,1))
!             evwk(1,1) = evl(nband)
              evwk(1:nband,1) = evl(1:nband)
              write(ifib,prfmt) ql, evwk(1:nblst,1)
              deallocate(evwk)
            else
              write(ifib,prfmt) ql, (evl(i),i=1,nblst)
            endif
          else  ! Write list of specific bands
            write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
          endif
        endif

C ... Contour plot mode
      elseif (mode == 3) then
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
        stdo = nglob('stdo')
        if (iprint() >= 60) write(stdo,345) i+1,j+1,ql
  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 == nqx .and. j+1 == nqy .and. jsp == nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
          do  k  = 1, nblst
          if (lbin) then
            write(ifib) nqx,nqy,1
            write(ifib) evsav(:,:,k,j1)
          else
            if (nsp == 1)
     .        call awrit4('%% rows %i cols %i  band=%i efermi=%;6d',
     .        ' ',80,ifib,nqx,nqy,iblst(k),eferm)
            if (nsp == 2)
     .        call awrit5(
     .        '%% rows %i cols %i  spin %i band=%i efermi=%;6d',
     .        ' ',80,ifib,nqx,nqy,j1,iblst(k),eferm)
            do  i = 1, nqx
              write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
            enddo
          endif
          enddo

          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

C ... Keep running tab on smallest,largest eval
      if (nevn /= 0) then
        evmxn = max(evmxn,evl(nevn))
        evmnn = min(evmnn,evl(nevn))
      endif
      return

      end

      subroutine suqlstst(s_lat,sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,qp,onesp,spintexture)
C- Initializes suqlst to up a list of q-points, and sets spin texture variables
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Passed parameters
      logical :: spintexture
      character*(*) sopts
      integer nq,nband,nsp,nfbn(3),ifblst(nband,3),onesp,iop
      double precision efermi,qp(3),evl(nband)

      call suqlst(s_lat,sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,qp,onesp)
      if (nfbn(1) < 0) then
        spintexture = .true.
        nfbn(1) = -nfbn(1)
      endif

      end

      subroutine suqlse(lio,nband,jsp,nsp,ndimhx,ifbn,nfbn,ifblst,ndlst,ndhamx,evec,wt)
C- Write to file the projection of eigenvector subblock for this qp
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio   :0  Generate weights and write to disk
Ci         :1  Generate weights, copy them to wt, and return
Ci         :2  Copy weights from wt, and write to disk
Ci   nband :number of energy bands to write
Ci   jsp   :current spin index (not used now)
Ci   nsp   :number of spins (not used now)
Ci   ndimhx:dimensions evec
Ci   ifbn  :index to color list (nfbn,ifblst)
Ci         :0 Spin texture (3 weights generated)
Ci         :>0 Single weight; use list ifblst(ifbn)
Ci   nfbn  :number of elements for color weights projection
Ci   ifblst:list of elements for  color weights projection
Ci   ndlst :leading dimension of ifblst
Ci   evec  :eigenvectors
Ci   ndhamx :dimensions wt
Cio Inputs/Outputs
Cio   wt   :lio = 0 => wt not used
Ci         :lio = 1 => wt is output
Ci         :lio = 2 => wt is input
Cl Local variables
Cl         :
Cr Remarks
Cr   Inefficient, but it works
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   23 Sep 09 Packaged Mulliken decomposition into suqlsz
Cu   08 Jul 08 New argument ndlst so ifblst can be dimensioned
Cu             independently from nband
Cu   05 Jun 06 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lio,jsp,nsp,nband,ifbn,nfbn(2),ndlst,ndhamx,ifblst(ndlst,4)
      integer ndimhx
      real(8) :: wt(ndhamx,4)
      double complex evec(ndimhx,ndimhx)
C ... Dynamically allocated arrays
      real(8), pointer :: qpl(:,:)
      real(8),allocatable :: evsav(:,:,:,:),wtl(:,:)
C ... Local parameters
      logical ltmp
      procedure(logical) :: isanrg
      character prfmt*40
      logical fmlong
      integer i,j,j1,k,stdo,iwt,nwt
      procedure(integer) :: iprint,nglob
C     Common with suqlst
      integer iblst(200),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
      double precision x1,x2,y1,y2,evmxn,evmnn,
     .  q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
      common /suqlsd/
     .  x1,x2,y1,y2,vectx,vecty,vect3,
     .  q1,q2,ql,qpl,evmxn,evmnn,ifiq,ifib,mode,iq,nevn,nblst,
     .  iblst,fmlong,nqx,nqy

      ltmp = isanrg(lio,0,2,' suqlse: ','lio',.true.)
      nwt = 1; if (ifbn == 0) nwt = 4
      allocate(wtl(ndimhx,4))
      if (lio == 2) then
        do  iwt = 1, nwt
          call dcopy(nband,wt(1,iwt),1,wtl(1,iwt),1)
        enddo
      else
        call suqlsz(ndimhx,ifbn,nfbn,ifblst,ndlst,evec,wtl)
C       print 444, ' gen color weight',wtl(1,:)
      endif
C  444 format(a,4f12.6)

      if (lio == 1) then
        do  iwt = 1, nwt
          call dcopy(nband,wtl(1,iwt),1,wt(1,iwt),1)
        enddo
        return
      endif

      do  iwt = 1, nwt

      if (mode == 1) then
        prfmt = '(3f10.5/(10f8.4))'
        if (nblst == 0) then ! General case, q-dependent nband
          write(ifib,prfmt) ql, (wtl(i,iwt),i=1,nband)
        elseif (nblst > 0 .and. iblst(1) == 0) then !Fixed list of bands
          write(ifib,prfmt) ql, (wtl(i,iwt),i=1,nblst)
        else ! Given list of bands
          write(ifib,prfmt) ql, (wtl(iblst(i),iwt),i=1,nblst)
        endif

      elseif (mode == 2) then
        if (fmlong) then
          prfmt = '(5f15.10)'
        else
          prfmt = '(8f10.6)'
C         if (nblst > 0 .and. nblst <= 5) prfmt = '(8f10.6)'
        endif
        if (nblst == 0) then ! General case, q-dependent nband
          write(ifib,prfmt) (wtl(i,iwt),i=1,nband)
        elseif (nblst > 0 .and. iblst(1) == 0) then !Fixed list
          write(ifib,prfmt) (wtl(i,iwt),i=1,nblst)
        else ! Given iblst
          write(ifib,prfmt) (wtl(iblst(i),iwt),i=1,nblst)
        endif

      elseif (mode == 3) then
        call rx('need copy weights into different place, mode=3')
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
        stdo = nglob('stdo')
        if (iprint() >= 60) write(stdo,345) i+1,j+1,ql
  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = wtl(iblst(k),iwt)
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 == nqx .and. j+1 == nqy .and. jsp == nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
          do  k  = 1, nblst
            if (nsp == 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
            if (nsp == 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)
            do  i = 1, nqx
              write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
            enddo
          enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif
      enddo ! iwt

      deallocate(wtl)

      end

      subroutine suqlsr(mode,ifi,nsp,nbf,nc,lde,nb1,nb2,ls1,ls2,nq,ef,qp,eb)
C- Count qp, optionally read them and energy bands from file
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit controls what is read.
Ci         :Any combination of the following is allowed:
Ci         :1 Find, and return:
Ci         :  nq = number of qp in file
Ci         :  nbf = maximum number of bands
Ci         :  nc  = column dimension of bnds file (qp mode)
Ci         :  ef = Fermi level, or NULL if missing from file
Ci         :2 read qp from file, return in array qp
Ci         :4 read bands in file; return in array eb
Ci         :10s digit flags the format of the qp file (see suqlst)
Ci         :0 bands in lines (default format is suqlst)
Ci         :1 same as 0
Ci         :2 qp list format (qp format, suqlst)
Ci   ifi   :read from logical unit ifi
Ci   nsp   :number of spins in band file
Ci   lde   :leading dimension of eb
Ci         :(only used when 4s bit of mode set)
Ci   nb1,nb2: read bands nb1..nb2 into eb
Ci         :(only used when 4s bit of mode set)
Ci   ls1   :(nsp=2): read first spin only
Ci         :(only used when 4s bit of mode set)
Ci   ls2   :(nsp=2): read second spin only
Ci         :(only used when 4s bit of mode set)
Cio Inputs/Outputs
Cio  nbf   :number of bands in band file
Cio        :nbf is input  if 1s bit of mode is zero
Cio        :nbf is output if 1s bit of mode is set
Cio  nc    :column dimension of band file (qp list mode)
Cio  nq    :number of k-points
Cio        :nq is input  if 1s bit of mode is zero
Ci         :   NB: the file nq cannot exceed input nq
Cio        :nq is output if 1s bit of mode is set
Co Outputs
Co   qp    :k-points, returned if 4s bit of mode is set
Co   eb    :energy bands, returned if 4s bit of mode is set
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   01 Jul 10 Bug fix, read mode 2, spin polarized case
Cu   12 Aug 09 Can read files in qp list format.  Also reads Ef.
Cu             New argument list
Cu   06 Jun 07 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ifi,nsp,nq,nbf,nc,lde,nb1,nb2
      logical ls1,ls2
      double precision ef,qp(3,nq),eb(lde,nsp,nq)
C ... Local parameters
      logical ltmp
      integer nq1,nql,mode1,mode2,mode4,iq1,i,nline,isp,fmode,ios
      integer nqvar,NULLI
      double precision ql(3),evl(3000),xx
      real(8),allocatable :: ewk(:,:)
      character strn*120
      parameter (NULLI=-99999)
      procedure(logical) :: rdstrn,a2bin,isanrg
      procedure(integer) :: rdm

C ... Setup
      if (mod(mode,10) == 0) return
      mode1 = mod(mod(mode,10),2)
      mode2 = mod(mod(mode,10)/2,2)
      mode4 = mod(mod(mode,10)/4,2)
      fmode = mod(mode/10,10)
      if (fmode == 0) fmode = 1

C --- Find in file: nq=# kpoints, nbf=maximum number of bands ---
      rewind ifi
C     In this mode, we can't determine nq from header info; do later
      if (fmode == 1) then
        read (ifi,*,iostat=ios) i,xx
        if (ios /= 0) call rx(
     .   'suqlsr (abort): failed to read number of bands (line format)')
        if (mode1 == 1) then
          nbf = i
          ef = xx
        endif
      elseif (fmode == 2) then
        if (mode1 == 1) then
C         nqvar = 1 => flag # of k-points depend on q
  825     if (.not. rdstrn(ifi,strn,len(strn),.false.)) goto 999
          if (strn(1:1) == '#') goto 825
          nqvar = 0
          if (strn(1:1) == '%' .and. index(strn,'nlb') > 0) nqvar=1
C         Find ef
          i = index(strn,'efermi=') + 6
          if (i == 6) i = index(strn,'ef=') + 2
          if (i > 2) then
            if (.not. a2bin(strn,ef,4,0,' ',i,-1))
     .        call rxs('SUQLSR: failed to parse Fermi level: ',strn(2:))
          else
            ef = NULLI
          endif
C         Check that nsp matches passed value.
C         If file value is missing, assume it matches passed nsp
          i = index(strn,'nsp=') + 3
          if (i > 3) then
            if (.not. a2bin(strn,isp,2,0,' ',i,-1))
     .        call rxs('SUQLSR: failed to parse nsp= ',strn(2:))
          else
            isp = nsp
          endif
          ltmp = isanrg(isp,nsp,nsp,'suqlsr:','file number of spins',.true.)

C         Find nq,nbf
          rewind ifi
          nql = 0
          nc = 0
          if (rdm(ifi,1000,0,' ',xx,nql,nc) /= 1) call rx(
     .   'suqlsr (abort): could not read bands from file (list format)')
          nq = nql / nsp
          if (mod(nql,2) /= 0 .and. nsp == 2) then
            call info0(20,0,0, ' suqlsr (warning): '//
     .        'odd number of bands encountered but nsp=2')
          endif
          nbf = nc-3 - nqvar
          if (nbf <= 0) call rx(
     .   'suqlsr (abort): could not read bands from file (list format)')
        endif
      else
        call rxi('suqlsr:  unrecognized file format mode',fmode)
      endif

      if (mode4 /= 0) then
C       Error if attempt to seek more bands than available
        ltmp = isanrg(nb2,1,nb1-1+min(lde,nbf),'suqlsr:','top band index',.true.)
C       Error if bottom band index out of range
        ltmp = isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
      endif

C --- Bands in line mode ---
      if (fmode == 1) then
C ... For each panel, do
      nql = 0
      nline = 0
   91 continue
      read(ifi,*) nq1
      if (nq1 <= 0) goto 90
      isp = 0
      do  iq1 = 1, nq1
        isp = mod(isp,2)+1
        nql = nql+1
        if (mode1 == 0 .and. nql > nq) call rxi(
     .    'suqlsr: file more q-points than allocated: nqmx=',nq)
        read(ifi,*,END=999,ERR=999) ql(1),ql(2),ql(3)
        if (mode2 /= 0) call dcopy(3,ql,1,qp(1,nql),1)
        read(ifi,*,END=999,ERR=999) (evl(i),i=1,nbf)
        if (mode4 /= 0) then
C         Copy only if appropriate spin
          if (isp == 1 .and. ls2) then
          elseif (isp == 2 .and. ls1) then
          else
            call dcopy(nb2-nb1+1,evl(nb1),1,eb(1,1,nql),1)
          endif
        endif
      enddo
      nline = nline+1
      goto 91
C     End of loop over lines
   90 continue

      if (mode1 == 0) then
        call info2(30,1,0,' suqlsr: found %i qp in %i lines from file',
     .    nql,nline)
      else
        nq = nql/nsp
        call info2(30,1,0,' suqlsr: read %i qp in %i lines from file',
     .    nql,nline)
      endif

      if (mod(nql,2) /= 0 .and. nsp == 2) then
        call info0(20,0,0, ' suqlsr (warning): '//
     .    'odd number of bands encountered but nsp=2')
      endif

C --- Bands in list mode ---
      else if (mode2 /= 0 .or. mode4 /= 0) then

        nql = nq*nsp
        allocate(ewk(nql,nc))

C       nqvar = 1 => flag # of k-points depend on q
        rewind ifi
  826   if (.not. rdstrn(ifi,strn,len(strn),.false.)) goto 999
        if (strn(1:1) == '#') goto 826
        nqvar = 0
        if (strn(1:1) == '%' .and. index(strn,'nlb') > 0) nqvar=1
        nc = nbf+3 + nqvar

C       Read bands into ewk; eb(1:nb2-nb1+1,1:nq) = ewk(1:nq,nb1:nb2)
        rewind ifi
        if (rdm(ifi,1000,nql*nc,' ',ewk,nql,nc) /= 1) call rxi(
     .    'suqlsr (abort): could not read bands from file, mode',fmode)
        do  i  = 1, nql
          if (mode2 /= 0) then
            qp(1:3,1+(i-1)/nsp) = ewk(i,1:3)
          endif
          if (mode4 /= 0) then
            call dcopy(nb2-nb1+1,ewk(i,3+nqvar+nb1),nql,eb(1,i,1),1)
          endif
        enddo
C       if (mode2 /= 0) call prmx('qp',qp,3,3,nq)
C       call prmx('eb',eb,nb2-nb1+1,nb2-nb1+1,nql)
        deallocate(ewk)

      endif
      return

  999 continue
      call rxi('suqlsr: failed to read bands file, nq=',nql)

      end

      subroutine suqlsz(ndimhx,ifbn,nfbn,ifblst,ndlst,evec,wt)
C- Return "spin texture" color weights for each eigenvector
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimhx:dimensions evec
Ci   ifbn  :index to color list (nfbn,ifblst)
Ci         :0  Spin texture (4 weights generated)
Ci         :>0 use ifblst(ifbn)
Ci   nfbn  :number of elements for color weights projection
Ci   ifblst:list of elements for  color weights projection
Ci   ndlst :leading dimension of ifblst
Ci   evec  :eigenvectors
Co Outputs
Co   wt    :wt(n,0:3) = weight for eigenvector n
Co         :weight 0   = Mulliken projection of sigma 0 (i.e. charge)
Co         :weight 1:3 = Mulliken projection of sigma x,y,z (i.e. spin)
Cl Local variables
Cr Remarks
Cr   Inefficient, but it works
Cr   Spin texture:
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma:
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr     M_z =  (rho11)-(rho22)
Cr   Second (symmetrized) form is used because for better stability
Cr
Cr   Debugging: break into spin components
Cr   Show normalized
Cr   mch -qr z -coll 4 -split z 1,1+289,1+578 1,2 -pop -qr z -i -t -coll 4 -split zd 1,1+289,1+578 1,2 -pop zd11 z11 -xe -rsum zd21 z21 -xe -rsum -+
Cr   11 - 22,  12, and 21
Cr   mch -qr z -coll 4 -split z 1,1+289,1+578 1,2 -pop -qr z -i -t -coll 4 -split zd 1,1+289,1+578 1,2 -pop zd11 z11 -xe -rsum zd21 z21 -xe -rsum --
Cr   mch -qr z -coll 4 -split z 1,1+289,1+578 1,2 -pop -qr z -i -t -coll 4 -split zd 1,1+289,1+578 1,2 -pop zd21 z11 -xe -rsum
Cr   mch -qr z -coll 4 -split z 1,1+289,1+578 1,2 -pop -qr z -i -t -coll 4 -split zd 1,1+289,1+578 1,2 -pop zd11 z21 -xe -rsum
Cu Updates
Cu   09 Jun 17 Add spin texture
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifbn,nfbn(2),ndlst,ifblst(ndlst,*)
      integer ndimhx
      double precision wt(ndimhx,0:3)
      double complex evec(ndimhx,ndimhx)
C ... Dynamically allocated arrays
      complex(8),allocatable:: evecc(:,:),work(:,:)
C ... Local parameters
      integer n,j,k,ifbni,ipiv(ndimhx),ndimh
      real(8) :: wt0
      complex(8) :: rho(2,2)
!     procedure(real(8)) :: dlength

C     call zprm('z',2,evec,ndimhx,ndimhx,ndimhx)
      allocate(evecc(ndimhx,ndimhx),work(ndimhx,ndimhx))
      call zcopy(ndimhx**2,evec,1,evecc,1)
      call zgetrf(ndimhx,ndimhx,evecc,ndimhx,ipiv,j)
      if (j /= 0) call rx('SUQLST: failed to generate overlap')
      call zgetri(ndimhx,evecc,ndimhx,ipiv,work,ndimhx**2,j)
C     call zprm('zinv',2,evecc,ndimhx,ndimhx,ndimhx)

      do  n = 1, ndimhx
        wt(n,0) = 0
C       Normal color weights
        if (ifbn > 0) then
          ifbni = ifbn
C       Spin texture color weights
        else
          ifbni = 1
          ndimh = ndimhx/2
!         wt0 = 0; call dpzero(wt,size(wt))
        endif

        rho = 0
        do  j = 1, nfbn(ifbni)
          k = ifblst(j,ifbni)
          if (k <= 0 .or. k > ndimhx .or. (ifbn == 0 .and. k > ndimh)) then
            call fexit2(-1,111,' Exit -1 : suqlst: component %i outside range (1:%i)',k,ndimhx)
          endif
          if (ifbn > 0) then
            wt(n,0) = wt(n,0) + evecc(n,k)*evec(k,n)
          else
C           wt0     = wt0     + evecc(n,k)*evec(k,n) + evecc(n,k+ndimh)*evec(k+ndimh,n)
C           wt(n,0) = wt(n,0) + 2*dble(evecc(n,k)*evec(k+ndimh,n))
C           wt(n,1) = wt(n,1) + 2*dimag(evecc(n,k)*evec(k+ndimh,n))
C           wt(n,2) = wt(n,2) + evecc(n,k)*evec(k,n) - evecc(n,k+ndimh)*evec(k+ndimh,n)

C           M_x =  2 Re(rho(2,1)) = Re (rho(1,2)+rho(2,1))
C           M_y =  2 Im(rho(2,1)) = Im (rho(2,1)-rho(1,2))
C           M_z =    rho(1,1)-rho(2,2)
            rho(1,1) = rho(1,1) + evecc(n,k)*evec(k,n)
            rho(1,2) = rho(1,2) + evecc(n,k)*evec(k+ndimh,n)
            rho(2,1) = rho(2,1) + evecc(n,k+ndimh)*evec(k,n)
            rho(2,2) = rho(2,2) + evecc(n,k+ndimh)*evec(k+ndimh,n)
          endif
        enddo
        if (ifbn <= 0) then
          wt0      =  rho(1,1)+rho(2,2) ! Mulliken charge
          wt(n,0)  =  wt0
          wt(n,1)  =  dble (rho(1,2)+rho(2,1))
          wt(n,2)  =  dimag(rho(2,1)-rho(1,2))
          wt(n,3)  =  (rho(1,1)-rho(2,2))
        endif

C        print 333, n,wt(n,0:2),dlength(3,wt(n,0),ndimhx)
C  333   format(i4,3f15.10,f18.15)
      enddo

      deallocate(evecc,work)

      end
      subroutine suqlsc(ipr,qcnst,nqin,qpin,nqout,qpout)
C- Pare out qp satisfying constraint
C ----------------------------------------------------------------------
Ci Inputs
Ci   qcnst :constraint (character expression containing iq,q,qx,qy,qz)
Ci   nqin  :number of initial q-points
Ci   qpin  :initial q-points
Co Outputs
Co   nqout :number of final q-points
Co   qpout :final q-points
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   16 Feb 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character qcnst*80
      integer ipr,nqin,nqout
      double precision qpin(3,nqin),qpout(3,*)
C ... Local parameters
      logical a2bin,lsw
      integer nq,iq,ival,iv0,ip,stdo,nglob
      double precision ddot,xx

      call numsyv(iv0)
      nq = 0

      do  iq = 1, nqin

        call lodsyv('iq',1,dble(iq),ival)
        call lodsyv('qx', 1,qpin(1,iq),ival)
        call lodsyv('qy', 1,qpin(2,iq),ival)
        call lodsyv('qz', 1,qpin(3,iq),ival)
        xx = sqrt(ddot(3,qpin(1,iq),1,qpin(1,iq),1))
        call lodsyv('q', 1,xx,ival)
C       call shosyv(0,0,0,6)

        ip = 0
        call rxx(.not. a2bin(qcnst,lsw,0,0,' ',ip,len(qcnst)),
     .    'suqlsc:  cannot parse exprn')
        if (lsw) then
          nq = nq+1
          call dcopy(3,qpin(1,iq),1,qpout(1,nq),1)
        endif

      enddo
      call clrsyv(iv0)

      call info(ipr,0,0,' suqlsc: qp list shortened to %i from %i '//
     .  'by constraint:  '//trim(qcnst),nq,nqin)
      if (ipr > 40) then
        stdo = nglob('stdo')
C       write (stdo,10)
        do  iq = 1, nq
          write (stdo,20) iq, (qpout(ip,iq),ip=1,3)
        enddo
C  10   format(19x,'qp')
   20   format(i4,4f10.6)
      endif
      nqout = nq

      end

      subroutine suqlsf(mode,ifiq,nq0,nq,iv,qp)
C- Read qp from file
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 just return number of qp in file, without reading points
Ci         :  In this mode nq is output
Ci         :1 read qp, given
Ci   ifiq  :file logical unit for reading
Ci   nq0   :offset: read qp into qp(:,nq0+1..nq)
Ci Inputs/Outputs
Cio  nq    :number of qp:
Cio        :Output (mode 0)
Cio        :Input  (mode 1): 2nd dimension of qp; must be at least
Cio        :                 as large as nq0 + # points read
Co Outputs
Co   iv    :iv(1:3) = nkabc (if given)
Co         :iv(4:6) = qp offsets (if given)
Co         :iv(7)   = number of tetrahedra (if given)
Co   qp    :q-points, stored in qp(:,nq0+1..nq)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Feb 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ifiq,nq,nq0
      double precision qp(3,nq)
C ... Local parameters
      logical rdstrn
      character strn*120
      integer iv(7),nq1,rdm,j
      double precision xx

C ... Determine file format
      rewind ifiq
   15 continue
C     Check for standard qp file format (defined in getqp.f)
      if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 999
      if (strn(1:1) == '#') goto 15 ! Skip past comments

C ... Case standard qp mode
      if (strn(1:1) /= '%' .and. index(strn,'nkp=') > 0) then
        call getqp(0,ifiq,nq1,iv,iv(4),iv(7),xx,xx,xx)
        if (mode == 0) then
          nq = nq1
          return
        endif
        if (nq0+nq1 > nq) call rx2('suqlsf: qp array too'//
     .    ' small: reading %i qp but dimensioned for %i',nq0+nq1,nq)
        call getqp(1,ifiq,nq1,iv,iv(4),iv(7),qp(1,nq0+1),xx,xx)
C ... Case generic qp file
      else
        j = 3
        nq1 = 0
        rewind ifiq
        if (rdm(ifiq,10000,0,' ',xx,j,nq1) /= 1) goto 999
        if (mode == 0) then
          nq = nq1
          return
        endif
        if (nq0+nq1 > nq) call rx2('suqlsf: qp array too'//
     .    ' small: read %i qp but dimensioned for %i',nq0+nq1,nq)
        rewind ifiq
        if (rdm(ifiq,10000,3*nq1,' ',qp(1,nq0+1),j,nq1) /= 1) goto 999

C        call awrit1('%N suqlst: read %i qp from file '//fn//
C     .    '%a',' ',80,stdo,nq)
C        if (j /= 1) call rx('suqlst: failed to read qp')
        call ivset(iv,1,7,0)
      endif
      return

C ... Error exit : return with error (nq<0)
  999 continue
      nq = -1

      end

      logical function lmapqp(smapq,lmapq,q,qmap)
C- Map q according to supplied strng expressions
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmapq :lmapq(i) = true => map ith cartesian component of q
Ci         :lmapq(i) = false => ith cartesian component of q untouched
Ci   smapq :ascii strings with algebraic expressions for mapping given q
Ci         :into another q.  Expression can contain q,qx,qy,qz.
Ci         :smapq holds a sequences of algebraic expressions containing
Ci         :one or more of:
Ci         :   kx=expr  ky=expr  kz=expr
Ci         :Example: smapq = q1=qx+qy kx=qx+.1, ky=q1
Ci         :Then qmap(1) = q(1)+.1, qmap(2) = q(1)+q(2),  qmap(3)=q(3)
Ci         :Expression kx=.. (ky=.. or kz=..) is required if lmapq(1)
Ci         : (lmap2(2) or lmapq(3)) is true.
Ci   q     :input q
Co Outputs
Co   qmap  :q(1:3), or transformed component i for any lmapq(i) = .true.
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   02 Jul 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character smapq*(*)
      logical lmapq(3)
      double precision q(3),qmap(3)
C ... Local parameters
      integer i,ival,ip,a2vec,n
      double precision ql(3),xx,ddot
      character vn*2

C     Store initial values into ql, so as not to overwrite q
C     in case qmap and q point to the same addres space
      lmapqp = lmapq(1) .or. lmapq(2) .or. lmapq(3)
      if (.not. lmapqp) return

      call numsyv(n)
      ql = q
      call lodsyv('qx', 1,q(1),ival)
      call lodsyv('qy', 1,q(2),ival)
      call lodsyv('qz', 1,q(3),ival)
      xx = sqrt(ddot(3,q(1),1,q,1))
      call lodsyv('q', 1,xx,ival)
C     call shosyv(0,0,0,6)
      ip = 0
      ip = a2vec(smapq,len(smapq),ip,4,', ',2,-3,-1,ival,xx)
      if (ip < 0)
     .  call rx(' suqlst: failed to parse expression: '//trim(smapq))
C     call shosyv(0,0,0,6)

      do  i  = 1, 3
        ql(i) = q(i)
        if (.not. lmapq(i)) cycle
        vn = 'k'//char(ichar('x')+i-1)
        call getsyv(vn,ql(i),ip)
        if (ip == 0)
     .    call rx(' suqlst: failed to parse expression: '//trim(smapq))
      enddo
      do  i = 1, 3
        qmap(i) = ql(i)
      enddo

      call clrsyv(n)
C     call shosyv(0,0,0,6)

      end
      subroutine rdqlst(mode,s_lat,strn,onesp,nkp,qp,npanel)
C- Read vector of qp from file for generating bands
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  suqlst
Ci Inputs
Ci   mode  :0 return nkp
Ci         :1 return qp
Ci         :3 return qp and npanel (npanel returned -1 unless symm line mode)
Ci         :Add 4 to 1 or 3 to require number of qp match nkp
Ci   strn  :string of band modifying arguments passed to suqlst
Ci   onesp
Cio Inputs/Outputs
Cio  nkp   :(output if mode=0) number of k-points
Cio        :(input if mode=1) dimensions qp
Co Outputs
Co   qp    :k-point
Co   npanel:In symmetry line mode:
Co         :npanel(0) = number of panels;
Co         :npanel(i) = cumulative number of qp up to panel i
Co         :In contour mode:
Co         :npanel(0) = -3
Co         :npanel(1) = number of divisions along x
Co         :npanel(2) = number of divisions along y
Co         :In the remaining modes, qp are grouped as one block.
Co         :npanel(0) = -1
Co         :The rest of the npanel is not used
Co         :Note: it is the caller's responsibility to dimension npanel
Co         :large enough to handle the intended use
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   16 Aug 17 Extended to qp-list and 2D contour modes
Cu   16 Nov 16 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!       include 'structures.h'
C ... For structures
      type(str_lat)::   s_lat
C ... Passed parameters
      integer mode,nkp,onesp,npanel(0:*)
      character strn*(*)
      double precision qp(3,nkp),qpi(3)
C ... Local parameters
      logical ltmp
      character(len=1) :: dc
      integer i,nfbn(4),ix(1),ifblst(3),iopq,iq,i1,nqi,stdo,ipan,modeql
      integer procid,master,numprocs
      double precision xx(1)
      procedure(integer) :: iprint,nglob,mpipid
      procedure(logical) :: isanrg

      procid = mpipid(1)
      master = 0
      numprocs = mpipid(0)
      stdo = nglob('stdo')
      dc = strn(1:1)

      if (mode == 0) then
        nkp = 0 ; i = 1 ; nfbn = 0
        if (procid == master) then
          iopq = 6
          call pshpr(1)
          xx = 0  ! So the compiler doesn't complain
          call suqlst(s_lat,strn,iopq,0,0d0,i,xx,nfbn,ifblst,nkp,qpi,onesp)
          call poppr
        endif
        call mpibc1(nkp,1,2,.false.,'rdqlst','nkp')
        if (nkp <= 0) call rx0('rdqlst: no k-points specified')
        return
      endif

C ... Generate new q-points for spectral function calculation
C ... Use nqi in place of nkp to preserve nkp
      if (procid == master) then
        call suqlsm(modeql)  ! q-points generation mode
C nqi -> i2 qpi->qp
C       iq = running index to big qp list,
C       nqi = number of points in current line
C       i1 = index to current line
C        call info5(30,1,0,' rdqlst: %i points in %?#(n==1)#%i panels#%j#'//
C     .    '%?#(n==2)#qp list mode##%?#(n==3)#2D contour mode##',
C     .    iq,modeql,npanel,modeql,modeql)

        call info2(30,1,0,' rdqlst: %i points in %?#(n==1)#symmetry line##'//
     .    '%-1j%?#(n==2)#qp list##%-1j%?#(n==3)#2D contour## mode',
     .    nkp,modeql)

        call pshpr(min(11,iprint()))
        i = 1 ; iq = 0 ; nqi = 0 ; iopq = 5; ipan=0
        do
          ipan = ipan+1
          xx = 0  ! So the compiler doesn't complain
          call suqlst(s_lat,strn,iopq,0,0d0,i,xx,nfbn,ix,nqi,qpi,onesp)
          if (mod(mode,4) == 3 .and. modeql == 1) then
            npanel(ipan) = 0
            npanel(0) = ipan-1
          endif
          if (nqi <= 0) exit
          if (mod(mode,4) == 3) npanel(ipan) = npanel(ipan-1) + nqi
          do  i1 = 1, nqi
            iq = iq+1
            if (iq > nkp) call rx('insufficient memory allocated for qp')
            call suqlst(s_lat,strn,iopq,0,0d0,i,xx,nfbn,ix,nqi,qpi,onesp)
            call dpscop(qpi,qp,3,1,3*iq-2,1d0)
          enddo
          if (mod(mode,4) == 3 .and. modeql /= 1) npanel(0) = -1
          if (mod(mode,4) == 3 .and. modeql == 3) then
            npanel(0) = -3
            call suqlsxy(npanel(1),npanel(2)) ! 2D contour mode
          endif
          if (modeql /= 1) exit
        enddo
        call poppr
        if (iprint() > 50) then
          do  i1 = 1, iq
            write(stdo,902) i1,qp(1:3,i1)
  902       format(1x,i5,3f12.6)
          enddo
        endif
      endif
      call mpibc1(iq,1,2,.false.,'','')
      call mpibc1(qp,3*iq,4,.false.,'rdqlst','qp')
      if (mode > 4) ltmp = isanrg(nqi-1,nkp,nkp,'rdqlst','nkp',.true.)
      end
