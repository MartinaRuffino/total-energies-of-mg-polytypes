      subroutine subzi(s_bz,lmet,ltet,lwt,nwt,ndham,nsp,nspc,nkp,zval,
     .  nevmx,lwtkb,efermi,lswtk,ef0)
C- Brillouin-integration setup
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: *
Co     Stored:    numq
Co     Allocated: wtkb swtk
Cio    Elts Passed:n wtkb w def
Cio    Passed to: *
Ci Inputs
Ci   lmet  :See Remarks
Ci         :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file, if they are needed
Ci         :3 always make two band passes; weights never needed a priori
Ci         :4 BZ integration with 3-point scheme
Ci         :5 same as lmet=3
Ci   ltet  :T
Ci   ltet  :0 => sampling integration
Ci         :1 => tetrahedron integration
Ci         :11 => tetrahedron integration for bands, sampling for weights
Ci         :ltet > 0 => allocate space for tetrahedron weights
Ci   lwt   :F weights are not needed until all bands are obtained
Ci         :T weights are needed a priori (eg output density generated)
Ci   nwt   :if wtkb is allocated, scale size of allocation by nwt.
Ci         :Normally nwt=1, but used in cases (e.g. L.S. pert theory)
Ci         :weights for multiple Fermi levels are required.
Ci   ndham :leading dimension of s_bz%wtkb, s_bz%swtk
Ci         :Hamiltonian should not exceed this dimension
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   zval  :total valence charge
Ci   ef0   :(for printout only) trial fermi level (sampling)
Ci   def   :(for printout only) Fermi level window
Cio Inputs/Outputs
Cio  nevmx :On input, maximum number of eigenvectors to find
Cio         input nevmx<0 => do not generate eigenvectors
Cio         input nevmx=0 => choose a default value
Co Outputs
Co   lwtkb :0 weights are neither required nor available a priori
Co         :1 weights are required a priori, and are read from disk
Co         :-1 weights are required a priori, but were not read
Co         :2 weights were generated with constrained global moment
Co         (not set by this routine, but see bzwtsf)
Co   efermi:(lwtkb=1) Fermi level corresponding weights; otherwise
Co         :efermi is set to -99.
Co   lswtk :Flags whether to make 'spin weights' swtk
Co         :-2 do not make spin weights
Co         : 1 spin weights array allocated; make them
Cl Local variables
Cl       n : if n>0, polynomial order for M-P sampling integration (bzwts.f)
Cl         : n<0 => Fermi distribution
Cr Remarks
Cr   To integrate the output density in the Brillouin Zone, integration
Cr   weights are needed for each band and qp, but they are not known
Cr   until all the bands are obtained.  The problem is solved in one of
Cr   the following ways:
Cr
Cr     lmet=0 system assumed to be an insulator; weights known a priori
Cr
Cr     lmet=1 eigenvectors are written to disk, in which case the
Cr            integration for the charge density can be deferred until
Cr            all the bands are obtained
Cr
Cr     lmet=2 integration weights are assumed from a prior band pass
Cr
Cr     lmet=3 two band passes are made; the first generates only evals
Cr
Cr     lmet=4 information is retained for three distinct Fermi levels.
Cr            After the Fermi level is determined, the density is
Cr            obtained by interpolation of the three points.  (This
Cr            scheme is suitable for sampling only, since in that case
Cr            just the Fermi level is needed to set integration weights.
Cr            When this scheme is used in conjunction with the
Cr            tetrahedron method, the charge density is calculated with
Cr            sampling.
Cu Updates
Cu   09 Aug 13 ltet now an integer, with allowed values 0,1,11
Cu   13 Aug 12 Removed oevl from argument list
Cu   08 Feb 10 (D. Pashov) lmet=5 added
Cu   09 Jun 07 Setup for spin weights (noncollinear case)
Cu   25 Apr 04 subzi returns efermi=-99, or file value if lwtkb=1.
Cu             Altered argument list
Cu   11 Oct 02 (ATP) MPI
Cu   21 Mar 01 Added printout; argument list changed
Cu   23 Jan 01 set numq=3 for lmet=4 always
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lwt
      integer ltet,lmet,ndham,nsp,nspc,nkp,lwtkb,lswtk,nevmx,nwt
      double precision zval,ef0,efermi
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
C ... Local parameters
      integer ifi,lerr,stdo,isw,n,i
      integer fopna,iobzwt,iprint,lgunit,idalloc,allocvb
      character*11 strni(2)
      integer procid,master,mpipid
      double precision xx
      data strni /'sampling','tetrahedron'/

      procid = mpipid(1)
      master = 0

      n = isign(1,s_bz%n) * mod(iabs(s_bz%n),100) ! Polynomial order
      efermi = -99
      lwtkb = 0
      s_bz%numq = 1
      if (lmet == 4) s_bz%numq = 3

      if (lmet > 0 .and. (lmet /= 4 .or. ltet == 1)) then
        i = max(mpipid(0),1)
        i = idalloc('wtkb',allocvb()+2,ndham*nsp,nkp*nwt*i)
        call ptr_bz(s_bz,8+1,'wtkb',ndham*nsp*nkp*nwt,0,xx)
      else
        call ptr_bz(s_bz,8+1,'wtkb',1,0,xx)
      endif
      if (nevmx >= 0) then
        if (lmet == 2 .or. lmet == 3 .or. lmet == 5) then
          lwtkb = -1
C     ... Attempt to use existing weights
          if (lmet == 2 .and. lwt) then
            if (procid == master) then
              ifi = fopna('wkp',-1,4)
              if (nspc == 1)
     .          lerr = iobzwt(0,ndham,nkp,nsp,efermi,s_bz%wtkb,ifi)
              if (nspc == 2)
     .          lerr = iobzwt(0,ndham*2,nkp,1,efermi,s_bz%wtkb,ifi)
              call fclose(ifi)
              if (lerr == 0) lwtkb = 1
              if (lerr /= 0) efermi = -99
            endif
C           Broadcast lwtkb,efermi,wtkb
            call mpibc1(lwtkb,1,2,.false.,'subzi','lwtkb')
            if (lwtkb == 1) then
              call mpibc1(efermi,1,4,.false.,'subzi','efermi')
              call mpibc1(s_bz%wtkb,ndham*nsp*nkp,4,.false.,'subzi','wtkb')
            endif
          endif
        endif
      endif
      lswtk = -2
      if (nspc == 2) then
        if (lwtkb == 1 .or. lmet == 5) lswtk = 1
        i = max(mpipid(0),1)
        i = idalloc('swtk',allocvb()+2,ndham*nsp,nkp*i)
        call ptr_bz(s_bz,8+1,'swtk',ndham*nsp*nkp,0,xx)
      else
        call ptr_bz(s_bz,8+1,'swtk',1,0,xx)
      endif

      if (nevmx == 0) then
        nevmx = (int(zval) + 1)/2
        nevmx = int(zval)/2+1
        if (lmet /= 0) nevmx = max(nevmx+nevmx/2,9)
        nevmx = min(nevmx,ndham)
        if (nspc == 2) nevmx = 2*nevmx
      endif

C ... Printout
      if (nevmx >= 0 .and. iprint() > 30) then
        stdo = lgunit(1)
        if (lmet > 0) then
          call awrit0('%N subzi: '//strni(isw(mod(ltet,2))+1)//
     .      '%a integration of bands; '//
     .      strni(isw(lmet /= 4.and.ltet == 1)+1)//
     .      '%a integration of density',' ',80,stdo)
          if (lmet == 4 .or. ltet /= 1) then
            call info8(30,0,0,
     .        '%7p sampling integration uses:  '//
     .        '%?#(n<0)#Fermi distribution (T=%;4d)#%-1jN=%i  W=%;4d#'//
     .        '  ef0=%;6;6d%?#(n==4)#  def=%;6d##',
     .        n,s_bz%w,ef0,lmet,s_bz%def,0,0,0)
          endif
        else
          call info(20,0,0,' subzi : nonmetal',0,0)
        endif
        write(stdo,'(1x)')
      endif

      end
      integer function iobzwt(mode,nevx,nq,nsp,efermi,wtkb,ifi)
C- File I/O of dos weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 File I/O header and band weights
Ci         :  Header consists of parms nevx,nq,nsp,efermi
Ci         :  For file read, nevx,nq,nsp must match passed values
Ci         :1 File I/O header information only
Ci         :  For file read, nevx,nq,nsp must match passed values
Ci   nevx  :leading dimension of wtkb
Ci   nq    :number of irreducible k-points (bzmesh.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  efermi:ifi>0: Fermi level written to disk
Cio        :ifi<0: Fermi level read from disk
Cio  wtkb  :ifi>0: weights for k-point integration written to disk
Cio        :ifi<0: weights for k-point integration read from disk
Co Outputs
Co   iobzwt:0  File I/O was successful
Co   iobzwt:-1 File I/O was not successful
Cr Remarks
Cu Updates
Cu   16 May 01 when writing header info only, write nq=0
Cu   18 Feb 01 Added Fermi level to file (and argument list)
Cu    5 May 00 Adapted from nfp pvwts
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nevx,nsp,nq,ifi
      double precision efermi,wtkb(nevx,nsp,nq)
C ... Local parameters
      logical ltmp
      integer jfi,nevx0,nq0,nsp0,iprint,lgunit
      procedure(logical) isanrg
C#ifdefC DEBUG
C      integer isp,iv,j
C      double precision dasum
C#endif

      ltmp = isanrg(mode,0,1,' iobzwt','mode',.false.)

      iobzwt = 0
C ... File read
      if (ifi > 0) then
        rewind ifi
        read(ifi,end=80,err=80) nevx0,nq0,nsp0,efermi
        if (mode == 1) then
          call info(20,1,0,
     .      ' Read efermi from weights file : ef = %,6;6d',efermi,0)
C#ifdefC DEBUG
C          do  isp = 1, nsp
C            if (isp == 2) print *, 'spin 2'
C            do  iv = 1, nevx
C            if (dasum(nq,wtkb(iv,isp,1),2*nevx) > 1d-6) then
C              write(*,'(i4,10f12.6/(4x,10f12.6))')
C     .          iv,(wtkb(iv,isp,j), j=1,nq)
C              endif
C            enddo
C          enddo
C          stop
C#endif
          return

        endif
        if (nevx /= nevx0 .or. nq /= nq0 .or. nsp /= nsp0) goto 80
        read(ifi) wtkb
C       call info(20,1,0,' Read qp weights ...  ef=%;6d',efermi,0)
        call info(20,1,0,' Read qp weights ...  ef=%;6d,  %;6d occ states',
     .    efermi,sum(wtkb))
C#ifdefC DEBUG
C        do  isp = 1, nsp
C          if (isp == 2) print *, 'spin 2'
C          do  iv = 1, nevx
C            if (dasum(nq,wtkb(iv,isp,1),2*nevx) > 1d-6) then
C              write(*,'(i4,10f12.6/(4x,10f12.6))')
C     .          iv,(wtkb(iv,isp,j), j=1,nq)
C            endif
C          enddo
C        enddo
C        stop
C#endif
        return
   80   continue
        iobzwt = -1
        call info(20,1,0,' Incompatible or missing qp weights file ...',
     .    0,0)
        return

C ... File write
      else
        jfi = -ifi
        rewind jfi
        if (mode /= 1) write(jfi) nevx,nq,nsp,efermi
        if (mode == 1) write(jfi) nevx,0,nsp,efermi
        if (mode == 1) then
          if (iprint() >= 20) then
            call awrit1('%N Saved Fermi level to weights file ... '//
     .        'ef = %,6;6d',' ',80,lgunit(1),efermi)
          endif
          return
        endif
        write(jfi) wtkb
        call info(20,1,0,' Saved qp weights ...',0,0)
      endif
      end

