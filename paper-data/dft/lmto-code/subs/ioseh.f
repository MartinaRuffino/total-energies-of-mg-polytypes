      subroutine ioseh(optio,mode,nsp,nspse,nband,nqp,nw,ommin,ommax,ef,nblst,iblst,nqlst,ifi)
C- Read/write header information in sigma file
C ----------------------------------------------------------------------
Ci Inputs
Ci  optio  :1s digit (read only)
Ci         :0 read header data into
Ci            nsp,nband,nqp,nw
Ci         :1 File data must passed matched value for
Ci            nsp,nband,nqp,nw
Ci         :2 file mode must passed matched value
Ci         :4 Read also ommin,ommax,ef from header
Ci         :These digits may be taken in any combination
Ci         :10s digit
Ci         :2 file is I/O in ascii format
Ci         :otherwise, file in binary format
Ci         :100s digit
Ci         :0 units in Ry
Ci         :1 units in eV
Ci         :2 (read only) accept whatever units file is written in
Cio        :  Returns in 100s digit optio either 0 or 1
Ci         :1000s digit (read only)
Ci         :1 Read nblst
Ci         :2 Read nqlst(0)
Ci         :3 combination 1+2
Ci         :Add 4 to read iblst (if nblst is read) and nqlst if nqlst(0) is read
Ci         :It is the user's responsbility to allocate iblst and/or nqlst
Ci   ifi   :File logical unit: ifi>0 for read; ifi<0 for write
Cio Inputs/Outputs
Cio  mode  :switches with information about the contents of sigma
Cio        :1s digit : how sigm(w) is recorded
Cio        :0 sigm(1:nw,1:nband,1:nq,1:nsp) written as one record
Cio        :1 sigm(1:nw,1:nband,iq,isp) in nq*nsp records
Cio        :2 sigm written user-specified format (usually sigm(1:nw,1:nband))
Cio        :10s digit: what frequency-independent self-energies are stored
Cio        :Add 1 to include sgq (QSGW Vxc)
Cio        :Add 2 to include exchange self-energy
Cio        :Add 4 to include LDA XC potential
Cio        :100s digit: what frequency-dependent object is stored
Cio        :0 diagonal parts of self-energy (complex)
Cio        :1 spectral function (real)
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio  nspse :2 for spin-polarized case and collinear bands & sigma; otherwise 1
Cio  nband :number of bands stored
Cio  nqp   :number of irr. k-points
Cio  nw    :number of frequency points.  nw<0 => Matsubara frequencies
Cio  ommin :first frequency.  nw<0 => inverse temperature
Cio  ommax :last frequency.   nw<0 => not used
Cio  ef    :Chemical potential or Fermi level
Cio  nblst :If nonzero iblst is not empty
Cio  iblst :if 0, bands ordered 1..nband
Cio        :otherwise, list of bands included in file
Cl Local variables
Cr Remarks
Cr   Reads or writes header for sigma(omega), with sanity checks.
Cr   sigma file contents: (1s digit mode = 0 and 1)
Cr    1.  header, 3 records
Cr        record 1: 0, version-number, mode, units
Cr                : 0 is a placeholder for future use
Cr                : mode is a switch indicating file contents
Cr                : units is a string indicating energy units of contents
Cr        record 2: nsp, nband, nqp, nw, ommin, ommax, mu (nw>0)
Cr                : nsp, nband, nqp, nw, beta, mu (nw<0)
Cr        record 3: '# iblst' list.  list is a standad integer list (mkilsd)
Cr        record 4: One of the following:
Cr                : '# qlst' n n1 n2 n3 (for band mode along symmetry lines)
Cr                :          n = number of panels, n1,n2,.. cumulative number of points in panel
Cr                :          n = 0 => no separate panels
Cr                :          In binary mode, n, n1,n2,n3 are written as one vector
Cr                : '# qlst con' n1 n2 (for contour mode on a raster of points)
Cr                :          Raster has n1 divisions along x, n2 along y.
Cr                :          In binary mode, a vector [-3 n1 n2] is written
Cr    2.  qp, written in one record
Cr    3.  QP levels, Eqp, written in one record
Cr    4.  static QP sigm(nband,nqp,nsp), in one record (if 10s digit mode contains 1)
Cr    5.  Fock exchange sigmx(nband,nqp,nsp), in one record (if 10s digit mode contains 2)
Cr    6.  LDA exchange-correlation vxc(nband,nqp,nsp), in one record (if 10s digit mode contains 4)
Cr    7.  sigm, written in one or multiple records, dependeng on 1s digit mode
Cu Updates
Cu   30 Apr 17 (MvS) Data can be read/written for frequencies on Matsubara axis
Cu   14 Apr 17 (v4) chemical potential is read/written
Cu   09 Apr 17 Slight change in mode : sigma0 is now optionally I/O; iose2 has extra argument
Cu   13 Nov 16 New modes formatB formatC; new 1000s digit mode, nblst,iblst,nqlst
Cu   15 Jul 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nsp,nspse,nband,nw,nqp,ifi,optio,mode,nblst,iblst(*),nqlst(0:*)
      double precision ommin,ommax,ef
C ... Local parameters
      integer jfi,i1,i2,i3,i4,i5,modef,optio0,optio3,lascii,lunits,iblst0,ix(100)
      integer k1,k2,k3,k4,ifvsn
      character units*2,unitstrn(0:1)*2
      double precision x1,x2,x3,scale
      character strn*1024
      common /iosehd/ scale
      integer, parameter :: NULLI=-99999, ivsn=5
      real(8), parameter :: rytoev=13.60569193d0
      procedure(integer) :: a2vec
      data unitstrn /'Ry','eV'/

      optio0 = mod(optio,10)
      lascii = mod(optio/10,10)
      lunits = mod(optio/100,10)
      optio3 = mod(optio/1000,10)

C --- File WRITE ---
      if (ifi < 0) then
        call sanrg(.true.,lunits,0,1,'ioseh:','lunits')
C       if (mode == 0) ivsn = 0
        jfi = -ifi
        rewind jfi
C       Write 0 => file contains version number (versions 1 and later)
        if (lascii == 2) then
          write(jfi,"('# header version mode  units')")
          write(jfi,"(i6,2i7,4x,a2)") 0,ivsn,mode,unitstrn(lunits)
          if (nw > 0) then
            write(jfi,"('# nsp nspse nband  nqp   nw       ommin',16x,'ommax',16x,'ef')")
            write(jfi,"(i4,4i6,2x,1p,3e21.13)") nsp,nspse,nband,nqp,nw,ommin,ommax,ef
          else
            write(jfi,"('# nsp nspse nband  nqp   nw       beta',16x,'ef')")
            write(jfi,"(i4,4i6,2x,1p,3e21.13)") nsp,nspse,nband,nqp,nw,ommin,ef
          endif
          if (nblst > 0) then
            call ilst2a(iblst,nblst,strn)
            write(jfi,"('# iblst ',a)") trim(strn)
          else
            write(jfi,"('# iblst 0')")
          endif
          if (nqlst(0) > 0) then
            call awrit2('# qlst %n:1i',' ',128,jfi,nqlst(0)+1,nqlst)
          elseif (nqlst(0) == -3) then
            call awrit2('# qlst con %2:1i',' ',128,jfi,nqlst(1),2)
          else
            write(jfi,"('# qlst 0')")
          endif
          write(jfi,"('# header end')")
        else
          if (ivsn /= 0) write(jfi) 0, ivsn, mode, unitstrn(lunits)
          if (nw > 0) then
            write(jfi) nsp,nspse,nband,nqp,nw,ommin,ommax,ef
          else
            write(jfi) nsp,nspse,nband,nqp,nw,ommin,ef
          endif
          if (nblst /= 0) then
            write(jfi) nblst, iblst(1:nblst)
          else
            write(jfi) 0
          endif
          if (nqlst(0) > 0) then
            write(jfi) nqlst(0:nqlst(0))
          elseif (nqlst(0) == -3) then
            write(jfi) nqlst(0:2)
          else
            write(jfi) 0
          endif
        endif

C --- File READ ---
      else
        call sanrg(.true.,lunits,0,2,'ioseh:','lunits')
        rewind ifi
C       Read version, checking whether header is pre- version 1
        if (lascii == 2) then
          read(ifi,*)
          read(ifi,*,err=99,end=99) i1,ifvsn,modef,units
        else
          read(ifi,err=99,end=99) i1,ifvsn,modef,units
        endif
        if (ifvsn < ivsn-1) call rx1('ioseh: require file version %i, sorry',ivsn)
        i2 = -1
        if (units == unitstrn(0)) i2=0
        if (units == unitstrn(1)) i2=1
        call sanrg(.true.,i2,0,1,'ioseh:','units')
        scale = 1
        if (i2 > lunits) scale = 1/rytoev
        if (i2 < lunits .and. lunits /= 2) scale = rytoev

C       Assign mode or check with file correspondence
        if (mod(optio0/2,2) /= 0) then
          call sanrg(.true.,modef,mode,mode,'ioseh:','file''s mode')
        else
          mode = modef
        endif
        if (lunits == 2) optio = optio + 100*i2 - 200  ! overwrite 100s digit optio with given units

C   ... Read header parameters or check with file correspondence
        if (lascii == 2) then
          read(ifi,*,err=99,end=99)
          read(ifi,*,err=99,end=99) i1,i2,i3,i4
          backspace ifi
          i5 = 1              ! nspse versions prior to v5
          if (i4 > 0) then ! Real axis: read nsp, nspse, nband, nqp, nw, ommin, ommax, mu
            if (ifvsn > 4) then
              read(ifi,*,err=99,end=99) i1,i5,i2,i3,i4,x1,x2,x3
            else
              read(ifi,*,err=99,end=99) i1,i2,i3,i4,x1,x2,x3
            endif
          else                ! Imaginary axis : read nsp, nspse, nband, nqp, nw, beta, mu
            if (ifvsn > 4) then
              read(ifi,*,err=99,end=99) i1,i5,i2,i3,i4,x1,x3
            else
              read(ifi,*,err=99,end=99) i1,i2,i3,i4,x1,x3
            endif
            x2 = NULLI
          endif
C         Read iblst
          read(ifi,"(a)",err=99,end=99) strn
          if (mod(optio3,2) == 1) then
            call word(strn,3,k1,k2)
            if (k1 > k2) call rx('ioseh: no ib list')
            if (strn(k1:k1) == '0') then
              nblst = 0
            elseif (optio3 > 4) then
              call mkilst(strn(k1:),nblst,iblst)
            else
              call mkils0(strn(k1:),nblst,k2)
            endif
          endif
C         Read nblst
          read(ifi,"(a)",err=99,end=99) strn
          if (mod(optio3/2,2) == 1) then
            call word(strn,3,k1,k2)

            if (strn(k1:k2) == 'con') then
              nqlst(0) = -3
              if (optio3 > 4) then
                k3 = a2vec(strn,len(strn),k2,2,', ',2,-3,2,ix,nqlst(1))
                if (k3 /= 2) call rx('ioseh: failed to read q list')
              endif
            else
              k4 = k1-1
              k3 = a2vec(strn,len(strn),k4,2,' ',1,1,1,ix,k2)
              if (k3 /= 1) call rx('ioseh: failed to read q list')
              nqlst(0) = k2
              if (optio3 > 4 .and. nqlst(0) > 0) then
                k3 = a2vec(strn,len(strn),k4,2,', ',2,-3,k2,ix,nqlst(1))
                if (k3 /= k2) call rx('ioseh: failed to read q list')
              endif
            endif
          endif
          read(ifi,*,err=99,end=99)
        else
          i5 = 1              ! nspse versions prior to v5
          read(ifi,err=99,end=99) i1,ifvsn,i3,i4
          if (i4 > 0) then ! Real axis: read nsp, nspse, nband, nqp, nw, ommin, ommax, mu
            if (ifvsn > 4) then
              read(ifi,err=99,end=99) i1,i5,i2,i3,i4,x1,x2,x3
            else
              read(ifi,err=99,end=99) i1,i2,i3,i4,x1,x2,x3
            endif
          else ! Imaginary axis : read nsp, nband, nqp, nw, beta, mu
            if (ifvsn > 4) then
              read(ifi,err=99,end=99) i1,i5,i2,i3,i4,x1,x3
            else
              read(ifi,err=99,end=99) i1,i2,i3,i4,x1,x3
            endif
            x2 = NULLI
          endif
          read(ifi,*,err=99,end=99) iblst0 ! Doesn't read iblst yet
          if (mod(optio3,2) == 1) nblst = iblst0
          read(ifi,*,err=99,end=99) iblst0 ! Doesn't read nqlst yet
        endif
        if (mod(optio0,2) == 0) then  ! No matching requirement
          nsp = i1
          nspse = i5
          nband = i2
          nqp = i3
          nw = i4
        endif
        if (mod(optio0/4,2) /= 0) then ! Read also ommin,ommax,ef
          ommin = x1*scale
          ommax = x2; if (x2 /= NULLI) ommax = x2*scale
          ef    = x3; if (x3 /= NULLI) ef = x3*scale
        endif
        call sanrg(.true.,i1,nsp,nsp,'ioseh:','file''s nsp')
        call sanrg(.true.,i5,nspse,nspse,'ioseh:','file''s nspse')
        call sanrg(.true.,i2,nband,nband,'ioseh:','file''s nband')
        call sanrg(.true.,i3,nqp,nqp,'ioseh:','file''s nqp')
        call sanrg(.true.,i4,nw,nw,'ioseh:','file''s nw')
      endif
      return

C --- Error exit ---
   99 continue
      call rx('ioseh: failed to read header from file')
      end

      subroutine iose2(optio,mode,rmode,nspse,nband,nqp,nw,qp,eig,sgq,sigwq,sex,vxcl,ifi)
C- I/O body of SE file
C ----------------------------------------------------------------------
Ci Inputs
Ci  optio  :1s digit not used
Ci         :10s digit
Ci         :2 file is I/O in ascii format
Ci         :otherwise, file in binary format
Ci         :100s digit not used
Ci         :1000s digit not used
Ci   mode  :switches with information about the contents of sigma
Ci         :1s digit : defines how sigm(w) is recorded
Ci         :0 sigwq(1:nw,1:nband,1:nq,1:nsp) written as one record
Ci         :1 sigwq(1:nw,1:nband,iq,isp) in nq*nsp records
Ci         :2 sigwq not I/O by iose2
Ci         :10s digit: what frequency-independent self-energies are stored
Ci         :Add 1 to include sgq (QSGW sigma)
Ci         :Add 2 to include exchange self-energy
Ci         :Add 4 to include LDA XC potential
Ci         :100s digit: what frequency-dependent object is stored
Ci         :0 diagonal parts of self-energy or Green's function in eigenfunction basis (complex)
Ci         :1 spectral function (real)
Ci  rmode  :(read only; only 10s digit is used)  Analogous to 10s digit mode.
Ci         :specifies what of QSGW sigma, exchange self-energy, and lda VXc is to be read
Ci  nspse  :2 for spin-polarized case and collinear bands & sigma; otherwise 1
Ci  nband  :number of bands to read
Ci  nqp    :number of irr. k-points
Ci  nw     :number of frequency points.  Not used with format C
Ci  ifi    :File logical unit: ifi>0 for read; ifi<0 for write
Cio Inputs/Outputs
Cio  qp    :k-point
Cio  eig   :QP levels, relative to the Fermi level
Cio  sigwq :freq-dependent correlation self-energy for each each band, q, spin
Cio  sgq   :value of sigwq at the QP peak
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   The omega-dependent self-energy has sigma(eqp) subtracted so that
Cr   sigwq(omega) vanishes at the QP peak eig and G has a pole there:
Cr     G = [omega - (eig + sigwq(omega) + i*eps)]^-1
Cu Updates
Cu   13 Nov 16 New mode formatB
Cu   15 Jul 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer optio,mode,rmode,nspse,nband,nqp,nw,ifi
      double precision qp(3,nqp),eig(nband,nqp,nspse),sgq(nband,nqp,nspse),
     .  sex(nband,nqp,nspse),vxcl(nband,nqp,nspse)
      complex(8):: sigwq(nw,nband,nqp,nspse)
C ... Local parameters
      integer lascii,isp,iq,mode0,mode1,mode2,n,rmode1
      integer, parameter :: formatA=0, formatB=1, formatC=2
      double precision scale
      common /iosehd/ scale

      lascii = mod(optio/10,10)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      rmode1 = mod(rmode/10,10)
      mode2 = mod(mode/100,10)
      n = 2 ; if (mode2 == 1) n = 1  ! size of element (complex or real)

      if (mode0 == formatA .or. mode0 == formatB .or. mode0 == formatC) then
        if (lascii == 2) then
          call iose3(' --- qp ---',ifi)
          call dfdump(qp,3*nqp,ifi)
          call iose3(' --- eig ---',ifi)
          call dfdump(eig,nband*nqp*nspse,ifi)
          if (mod(mode1,2) >= 1) then
            call iose3(' --- sgq ---',ifi)
            call iose4(sgq,nband*nqp*nspse,.false.,mod(mode1,2)==mod(rmode1,2),ifi)
          endif
          if (mod(mode1/2,2) >= 1) then
            call iose3(' --- sex ---',ifi)
            call iose4(sex,nband*nqp*nspse,.false.,mod(mode1/2,2)==mod(rmode1/2,2),ifi)
          endif
          if (mod(mode1/4,2) >= 1) then
            call iose3(' --- vxc ---',ifi)
            call iose4(vxcl,nband*nqp*nspse,.false.,mod(mode1/4,2)==mod(rmode1/4,2),ifi)
          endif
          if (mode2 == 0) then
            call iose3(' --- sig ---',ifi)
          else
            call iose3(' --- A ---',ifi)
          endif
          if (mode0 == formatA) then
            call dfdump(sigwq,n*nw*nband*nqp*nspse,ifi)
          elseif (mode0 == formatB) then
            do  isp = 1, nspse
            do  iq = 1, nqp
              call dfdump(sigwq(1,1,iq,isp),n*nw*nband,ifi)
            enddo
            enddo
          endif
          if (scale /= 1 .and. ifi > 0) then
            call dscal(nband*nqp*nspse,scale,eig,1)
            call dscal(nband*nqp*nspse,scale,sgq,1)
            if (mode0 == formatA .or. mode0 == formatB) then
              call dscal(n*nw*nband*nqp*nspse,scale,sigwq,1)
            endif
          endif
C          print *, sgq(nband,nqp,nspse)
C          print *, sigwq(nw,nband,nqp,nspse)
C          print *, sigwq(nw,nband,nqp,nspse)-sgq(nband,nqp,nspse)
        else
          call dpdump(qp,3*nqp,ifi)
          call dpdump(eig,nband*nqp*nspse,ifi)
          if (mod(mode1,2) >= 1)
     .      call iose4(sgq,nband*nqp*nspse,.true.,mod(mode1,2)==mod(rmode1,2),ifi)
          if (mod(mode1/2,2) >= 1)
     .        call iose4(sex,nband*nqp*nspse,.true.,mod(mode1/2,2)==mod(rmode1/2,2),ifi)
          if (mod(mode1/4,2) >= 1)
     .      call iose4(vxcl,nband*nqp*nspse,.true.,mod(mode1/4,2)==mod(rmode1/4,2),ifi)
          if (mode0 == formatA) then
            call dpdump(sigwq,n*nw*nband*nqp*nspse,ifi)
          elseif (mode0 == formatB) then
            do  isp = 1, nspse
            do  iq = 1, nqp
              call dpdump(sigwq(1,1,iq,isp),n*nw*nband,ifi)
            enddo
            enddo
          endif
        endif
      else
        call rx('iose2: unrecognized mode')
      endif

      end
      subroutine iose3(strn,ifi)
      implicit none
      integer ifi
      character *(*) strn
      character*80 strn2

      if (ifi < 0) write(-ifi,'(a)') strn
      if (ifi > 0) read(ifi,'(a)') strn2
      if (trim(strn) /= trim(strn)) call rx('oops')
      end

      subroutine iose4(arr,n,lbin,lread,ifi)
C- Read from file, or read it and discard
      integer n
      logical lbin,lread
      real(8) :: arr(n)
      real(8), pointer :: wk(:)

      if (.not. lread) then
         allocate(wk(n))
         if (.not. lbin) then
           call dfdump(wk,n,ifi)
         else
           call dpdump(wk,n,ifi)
         endif
       else
         if (.not. lbin) then
           call dfdump(arr,n,ifi)
         else
           call dpdump(arr,n,ifi)
         endif
       endif

      end
