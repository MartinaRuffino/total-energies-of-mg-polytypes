      subroutine mksyml(instr,plat)
C- Build a symmetry lines from a string
C ----------------------------------------------------------------------
Ci Inputs
Ci   instr :string defining symmetry lines; see Remarks
Ci   plat  :primitive lattice vectors, in units of alat
Co Outputs
Co   symmetry lines file is written to disk.
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  This routine writes a set of symmetry lines to file 'syml.ext', using
Cr   instr to specify the symmetry lines.  The syntax is described below.
Cr  *instr must contain a list of at least two q-points q0 and q1, with the
Cr   syntax described below.  There may be as many q-points as desired.
Cr  *Symmetry lines are constructed for lines joining (q0,q1), (q1,q2), ...
Cr   until all the q-points are exhausted.
Cr  *The number of divisions on a line is determined by ndiv*length(qi+1 - qi).
Cr   ndiv may be optionally set; it defaults to 31.
Cr  *Labels may optionally be given.
Cr   They are appended to the end of each symmetry line.
Cr  *You can specify whether the q you supply are in Cartesian coordinates (default)
Cr   or as multiples of reciprocal lattice vectors.
Cr  *You can specify whether the q written to syml.ext are in Cartesian coordinates
Cr   (default) or as multiples of reciprocal lattice vectors.
Cr
Cr   Syntax:
Cr  *The minimum syntax for instr is q=#0x,#0y,#0z,#1x,#1y,#1z
Cr   which specifies the Cartesian components of q0 and q1.
Cr   The number of qp can be enlarged as desired by apppending ,#2x,#2y,#2z ...
Cr  *Options modifying the file are also specified through instr.
Cr   The various options are separate by a delimiter.  If the first character
Cr   of instr is not a letter, it becomes the delimiter.  Otherwise '~' is used.
Cr
Cr  *The full syntax is
Cr      [~n=ndiv][~lbl=KGLMGA][~mq][~wmq]~q=#0x,#0y,#0z,...
Cr   Option n=... specifies ndiv
Cr   Option lbl specifies labels for each of the q
Cr   Option mq tells this routine to rotate input q to Cartesian coordinates
Cr   Option wmq tells this routine to rotate output q to multiples of G
Cr  *Alternate syntax:
Cr      [~n=ndiv][~lbl=KGLMGA][~mq][~wmq]~lblq:S=#,#,#,S=#,#,#,...
Cr   In this case, q= is not used, but symmetry lines taken from lbl=.
Cr   You must specify q-points corresponding to each character in lbl with tag lblq.
Cr
Cr   Note that the joint use of ~mq~wmq does not affect the symmetry lines,
Cr   but it does affect the length used to compute the number of divisions.
Cr
Cr   Example:  symmetry lines for a hexagonal lattice
Cr      ~n=41~wmq~mq~lbl=KGLMGA~q=1/3,1/3,0,0,0,0,1/2,0,1/2,1/2,0,0,0,0,0,0,0,1/2
Cr   Equivalent using Alternate format:
Cr      ~n=41~mq~wmq~lblq:G=0,0,0,A=0,0,1/2,L=1/2,0,1/2,K=1/3,1/3,0,M=1/2,0,0~lbl=KGLMGA
Cu Updates
Cu   20 Oct 18  New alternate format (lblq)
Cu   11 Apr 18  Documentation and switch ~wmq added.
Cu   25 Feb 18  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) instr
      double precision plat(3,3)
C ... Dynamically allocated local arrays
!       character (len=:), allocatable :: outs ! when it works reliably (broken in gcc 9.1 on macos 10.14)
      character(len=2047) :: outs
      real(8), allocatable :: qp(:,:)
C ... Local parameters
      integer, parameter :: nqmx = 100
      character*1 dc,dc2,msg*64,fn*128,lbl(nqmx),lblq*(nqmx),sout*128
      integer i,ndiv,k,j1,ivec(nqmx),nq,nlbl,ifi,n
      double precision xv(3*nqmx),qlat(3,3)
      procedure(integer) :: wordsw,a2vec,fopna
      procedure(real(8)) :: dlength

C     Defaults
      ndiv = 31
      call mkqlat(plat,qlat,xv)
      if (instr == ' ') then
        msg = 'no given'
        goto 99
      endif
      lbl = ' '; lblq = ' '
      nq = 0

C     Get delimiter: use '~' if first character is a letter
      k = len(instr)+1
!       allocate(character(k)::outs)
      dc = instr(1:1)
      if (dc >= 'a' .and. dc <= 'z') then
        dc = '~'
        outs(1:1) = dc
        outs(2:) = instr
      else
        outs = instr
        k = len(outs) ! should not be necessary; avoid compiler bug
        outs(k:k) = ' '
      endif

C ... Get label-q-point pairs
      k = wordsw(outs,dc,'lblq','',j1) + 3
      nlbl = 0
      if (k > 2) then
        msg = 'failed to parse'
        dc2 = outs(j1:j1+1)
        if (dc2 == ' ') goto 99
        k = j1+2
        do  while (outs(k:k) == '=')
          nlbl = nlbl+1
          lblq(nlbl:nlbl) = outs(k-1:k-1)
          i = a2vec(outs,len_trim(outs),k,4,', '//dc//dc2,4,2,3,ivec,xv(3*nlbl-2))
          if (i /= 3) goto 99
          k = k+2
        enddo
      endif

C ... Get list of qp
      msg = ' '
      k = wordsw(outs,dc,'q=','',j1) + 2
      if (k > 2 .and. nlbl == 0) then
        i = a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,99,ivec,xv)
        if (mod(i,3) /= 0) msg = 'q vector must have multiple of 3 elements,'
        if (i < 6) msg = 'q list must contain at least 2 vectors,'
        if (i < 0) msg = 'failed to parse'
        nq = i/3
        allocate(qp(3,nq))
        call dcopy(3*nq,xv,1,qp,1)
      elseif (k > 2) then
        msg = 'conflict : select only one of "q=" and "lblq"'
      elseif (nlbl == 0) then
        msg = 'string must contain "q="'
      endif
      if (msg /= ' ') goto 99

C ... Get labels associated with each q
      k = wordsw(outs,dc,'lbl=','',j1) + 2
      if (k > 2) then
        if (nlbl > 0) then ! Count number of lines
          allocate(qp(3,nqmx))
          do  i = j1, len(outs)
            if (verify(outs(i:i),dc//' ') == 0) exit
            if (outs(i:i) == ',') cycle
            nq = nq+1
            k = scan(lblq,outs(i:i)); if (k == 0) exit
            lbl(nq) = outs(i:i)
            call dcopy(3,xv(3*k-2),1,qp(1,nq),1)
          enddo
          if (nq < 2) msg = 'q list must contain at least 2 vectors,'
          if (k == 0) msg = 'point '//outs(i:i)//' missing from lblq'
        else
          do  i = 1, nq
            do while (outs(j1:j1) == ','); j1 = j1+1; enddo
            lbl(i) = outs(j1:j1)
            j1 = j1+1
          enddo
        endif
      elseif (nlbl > 0) then
        msg = '"lblq" requires "lbl=" be present also'
      endif
      if (msg /= ' ') goto 99

C ... Get number of k divisions
      ndiv = 51
      k = wordsw(outs,dc,'n=','',j1) + 2
      if (k > 2) then
        if (a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,1,ivec,ndiv) < 1)
     .    msg = 'failed to parse'
      endif
      if (msg /= ' ') goto 99

C ... Convert qp to Cartesian coordinates
      if (wordsw(outs,dc,'mq','',j1) > 0) then
        call info0(1,0,0,'%9frotating qp to Cartesian coordinates ...')
        do  i = 1, nq
          call dcopy(3,qp(1,i),1,xv,1)
          call dgemm('N','T',1,3,3,1d0,xv,1,qlat,3,0d0,qp(1,i),1)
        enddo
      endif

C ... Convert qp to multiples of lattice vectors
      if (wordsw(outs,dc,'wmq','',j1) > 0) then
        call info0(1,0,0,'%9frotating qp to Multiples of lattice vectors ...')
C       Do rotation later to preserve lengths
C       do  i = 1, nq
C         call dcopy(3,qp(1,i),1,xv,1)
C         call dgemm('N','N',1,3,3,1d0,xv,1,plat,3,0d0,qp(1,i),1)
C       enddo
      endif

C ... Write syml file
      fn = 'syml'
      ifi = fopna(trim(fn),-1,0)
      call awrit0('# generated from text '//trim(instr),' ',len(instr)+25,ifi)
      do  i = 1, nq-1
        n = max(nint(ndiv*dlength(3,qp(:,i)-qp(:,i+1),1)),2)
        xv(1:3) = qp(:,i)
        xv(4:6) = qp(:,i+1)
        if (wordsw(outs,dc,'wmq','',j1) > 0) then
          call dgemm('N','N',1,3,3,1d0,qp(1,i),1,plat,3,0d0,xv,1)
          call dgemm('N','N',1,3,3,1d0,qp(1,i+1),1,plat,3,0d0,xv(4),1)
        endif
        call awrit3('%x%,4i %3;12,7D  %3;12,7D',sout,len(sout),0,n,xv,xv(4))
        if (lbl(i) /= ' ') call awrit0('%a'//'  '//lbl(i)//'  to  '//lbl(i+1),sout,len(sout),0)
        write(ifi,'(a)') trim(sout)
      enddo

      return
   99 continue
      call rx('mksyml: '//trim(msg)//' in text '//trim(instr))
      end

