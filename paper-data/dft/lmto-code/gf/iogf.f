      integer function iogfrs(lio,igf,errh,fnam,ifi,ic1,ic2,nbas,ntab,
     .  zp,qp,plat,pos,clp,clssl,iaxg,offg,liaxb,kcplx,gf)
C- Real-space / cluster Green function read/write
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Ci   lio  :controls file handling and header I/O
Ci         1s digit: file handling
Ci           0 for read
Ci           1 for write
Ci           2 add to use given ifi instead of obtaining it from fnam
Ci           4 add to suppress rewinding of file before I/O
Ci         10s-1000s digit: bits set conditions for file read:
Ci           1 (=MNC) file ic1,ic2 must match passed values
Ci             If bit is not set, ic1,ic2 are assigned to file values
Ci           2 (=MNB) file nbas must match passed value
Ci             If bit is not set, assigned to file value
Ci           4 (=MNT) file ntab must match passed value
Ci             If bit is not set, assigned to file value
Ci           8 (=MLB) file liaxb must match passed value
Ci             If bit is not set, assigned to file value
Ci          16 (=MPL) file plat,pos match passed values
Ci          32 (=MCD) file dimensions of g match those in clp
Ci          64 (=MZP) file zp,qp, and gf type match passed value
Ci                    gf type is specified by 100000s digit lio, below
Ci         128 (=MSC) spin coupled case.
Ci                    If set, true clp(3..5) correspond to orbital degrees
Ci                    of freedom only, and true dimensions are 2*clp(3..5)
Ci        10000s to 1000000s digit: bits flag what header info is I/O
Ci           1 (=ISH)  skip file opening and header I/O
Ci           2 (=ICLP) I/O clp(ic1:ic2)
Ci           4 (=ICLSSL) I/O  clssl
Ci           8 (=IIAX) I/O iaxg(ic1:ic2)
Ci          16 (=IPOS) I/O plat,pos
Ci          32 (=IDGF) sets how second dimension ndima of gf is fixed.
Ci                     bit=0 means ndima=clp(4,ic1)
Ci                     bit=1 means ndima=(clp(6,icl+1)-clp(6,icl))/ndimb
Ci
Ci   igf   :flags what is contained in array gf
Ci           0 no GF is I/O
Ci           1 gf contains g
Ci           2 gf contains s*g
Ci           3 gf contains S*g*S+
Ci
Ci   errh :controls error handling when data is missing or a required
Ci         match (described in lio, above) is not satisfied.
Ci         0 flags iogfrs to abort
Ci        >0 iogfrs returns <0 and prints a warning if verbosity > errh
Ci
Ci   fnam :file name holding gf.  See description of argument ifi.
Ci
Cio  ifi  :File logical unit.
Cio       *If 1s digit of lio contains 2's bit,
Cio        OR if 1000s digit of lio contains 1's bit,
Cio        iogfrs uses ifi as the file handle and does not open the file
Cio       *Otherwise, if 1s digit of lio does not have 2's bit set,
Cio        the file named 'fnam' is opened and ifi is assigned
Cio        the logical unit connected with that file.
Cio
Cio    ... On file read, the following are set to disk values as read
Cio        from the first record, or else the passed values must match
Cio        the disk values, depending on 10s+100s digit of lio.
Cio
Cio  ic1  :I/O starting with cluster ic1
Cio  ic2  :I/O ending with cluster ic2. ic2=0 => no cluster, just atom
Cio  nbas :dimensions clssl,pos.
Cio  ntab :dimensions iaxg.
Cio  zp   :complex energy
Cio  qp   :k-point
Cio
Cio    ... The following arrays are I/O on disk when appropriate
Cio        bits in lio are set.  On file read, iogfrs aborts or returns
Cio        a negative value when array is sought but missing from file.
Cio        errh controls how iogfrs handles read errors.
Cio  clp  :index and dimensioning information for R.S. cluster (ioordn.f)
Cio       :form of lmgf and also for the layer GF code lmpg
Cio       :In the latter case pgplp for a given PL can substitute for clp
Cio       :and only 1:5 are used.  (pass ic2=0 in that case)
Cio       : 1 cumulative number of basis atoms in CL=1 ... this CL
Cio           lmpg: number of basis atoms in all PL prior to this one
Cio       : 2 last pair in iaxg table for this cluster
Cio       :   lmpg: index labelling kind of potential of this PL
Cio       : 3 column (source) matrix dimension for this CL
Cio       :   In the spin-coupled case, the true dimension is 2*clp(3)
Cio       : 4 row (field) matrix dimension for this CL
Cio       :   In the spin-coupled case, the true dimension is 2*clp(4)
Cio       : 5 matrix dimension for direct zone of CL
Cio       :   lmpg: matrix dimension for this PL including i-waves
Cio       :   In the spin-coupled case, the true dimension is 2*clp(5)
Cio       : 6 offset to g for start of this cluster
Cio       : 7 last pair in iaxg table for direct zone, inner block
Cio       : 8 last pair in iaxg table for direct zone, outer block
Cio       : 9 last pair in iaxg table for which g is calculated
Cio  clssl:cluster list: site ib belongs to cluster clssl(ib)
Cio        clssl need not be passed if not I/O
Cio  iaxg :neighbor table for a cluster (paircl.f)
Cio        clssl need not be passed if not I/O
Ci
Ci   offg :offset to first entry in gf
Cio  liaxb:information about the cluster table
Cio        1s digit 0 => iax(1,*) holds site index
Cio        1s digit 1 => iax(1,*) holds cluster index
Cio  kcplx:mode of complex storage for g
Cio        0: g stored with imaginary following real
Cio        1: g in complex*16 format
Cio        2: g have real, imaginary separated by columns
Cio
Cio  gf   :Real-space or cluster Green's function, or related matrix
Cio        is I/O to disk
Cio        gf is dimensioned gf(ndima,ndimb) with ndimb = clp(3,icl)
Cio        and ndima = clp(4,icl) or (clp(6,icl+1)-clp(6,icl))/ndimb
Cio        depending how lio is set; see IDGF, above.
Co iogfrs :(file read) 0 if read is successful
Cr Remarks
Cr  *File consists of 1st record (header) containing :
Cr     *Principal Layer mode
Cr     lio   file's contents; see lio above
Cr     liaxb information about the cluster table; see liaxb above
Cr     *cluster mode
Cr     ic1 = 1st entry in cluster table to be I/O
Cr     ic2 = last entry in cluster table to be I/O. ic2=0 => no cluster, just atom
Cr     lio   file's contents; see lio above
Cr     liaxb information about the cluster table; see liaxb above
Cr     nbas  file's nbas
Cr     ntab  file's number of pairs in pair table
Cr  *Whether the next records are present or not depends on lio:
Cr     clp
Cr     clssl
Cr     iaxg
Cr     plat
Cr     pos
Cr  *Next follows a list of gf, the contents of each consisting of:
Cr     ndima, ndimb, kcplx, 0, igf, zp, qp  <- dimensioning info & energy
Cr                                            l6 = 100000s of lio
Cr     gf                                  <- Green's function array
Cu Updates
Cu   01 Feb 17 New switch IPLMOD to distinguish PL mode from cluster
Cu   01 May 16 New switch MSC for noncollinear case
Cu   26 Feb 02 iogfrs can read kcplx=0 and convert to kcplx=1
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lio,errh,niax,nclp,ic1,ic2,nbas,ntab,offg,kcplx,liaxb
      parameter (niax=10,nclp=9)
      integer iaxg(niax,1),clp(nclp,ic2+1),clssl(nbas)
      character*(*) fnam
      double precision zp(2),qp(3),gf(*),plat(3,3),pos(3,nbas)
C ... Local parameters
      logical ltmp
C     Bits for matching conditions (10-100s digit of lio)
      integer, parameter :: MNC=1,MNB=2,MNT=4,MLB=8,MPL=16,MCD=32,MZP=64,MSC=128
C     Bits for selecting which data to I/O (10000-100000s digit of lio)
      integer ISH,ICLP,ICLSSL,IIAX,IPOS,IDGF,IPLMOD
      parameter (ISH=1,ICLP=2,ICLSSL=4,IIAX=8,IPOS=16,IDGF=32,IPLMOD=64)
      integer icl,ifi,ipr,lio0,lio123,lio456,igf,ncp1,ndimb,ndima,ndimbl,ndimal,
     .  ndims,offnow,kio,kio123,kio456,ir(6),ntabl,kcplxl,cast,ii,lerr,igfl
      equivalence (ir(3),kio),(ir(5),ntabl)
      double precision zpl(2),qpl(3),tol
      parameter (tol=2d-7)
      character*8 locnam, readin*6
      procedure (logical) :: bittst
      procedure(integer) fopnx,lgunit,rdarry,isw

C ... Setup
      iogfrs = -1
      call getpr(ipr)
      call sanrg(.true.,kcplx,0,2,'iogfrs:','kcplx')
      locnam = fnam
      lio0   = mod(lio,10)
      lio123 = mod(lio/10,1000)
      lio456 = mod(lio/10000,1000)

      if (.not. bittst(lio456,ISH)) then
        if (.not. bittst(lio0,2)) ifi = fopnx(fnam,100,16+8+4+0,-1)
        if (lio0 < 4) rewind ifi
      endif
      ncp1 = ic2-ic1+2  ! Not used in IPL mode

C --- File write ---
      if (mod(lio0,2) == 1) then

C   ... Write header
        if (.not. bittst(lio456,ISH)) then
          if (bittst(lio456,IPLMOD)) then
            write(ifi) lio-mod(lio,10),nbas
            if (bittst(lio456,ICLP)) call ipdump(clp(1,ic1),6,-ifi)
          else
            ntabl = clp(2,ic2+1)
            write(ifi) ic1,ic2,lio-mod(lio,10),liaxb,nbas,ntabl
            if (bittst(lio456,ICLP)) call ipdump(clp(1,ic1),nclp*ncp1,-ifi)
            if (bittst(lio456,ICLSSL)) call ipdump(clssl,nbas,-ifi)
            if (bittst(lio456,IIAX))call ipdump(iaxg(1,ic1),niax*ntabl,-ifi)
          endif
          if (bittst(lio456,IPOS)) call dpdump(plat,9,-ifi)
          if (bittst(lio456,IPOS)) call dpdump(pos,3*nbas,-ifi)
        endif

C   ... Write gf
        ndims = 0
        if (igf /= 0) then
          offnow = offg
          do  icl = ic1, ic2
            ndimb = clp(3,icl)
            if (bittst(lio456,IDGF)) then
              ndima = (clp(6,icl+1) - clp(6,icl)) / ndimb
            else
              ndima = clp(4,icl)
            endif
C           Write cast instead of kcplx to comply with rdm.
            cast = kcplx + 2
            write(ifi)  ndima, ndimb, cast, 0, igf, zp, qp
            if (bittst(lio123,MSC)) then
              ndima = ndima*2
              ndimb = ndimb*2
              if (offg /= 0) call rx('iogfrs check offg nc case')
            endif
            call dpsdmp(gf,1+2*offnow,2*(offnow+ndima*ndimb),-ifi)
            offnow = offnow + ndima*ndimb
            ndims = ndims + ndima*ndimb
          enddo
        endif

        if (ipr >= 40 .and. igf > 0) call awrit3(
     .    ' IOGFRS: wrote to file '''//locnam//
     .    '%a'' %i elements from clusters %i...%i',' ',80,lgunit(1),
     .    ndims,ic1,ic2)

C --- File read ---
      else
        if (.not. bittst(lio456,ISH)) then
          readin = 'header'
          if (bittst(lio456,IPLMOD)) then
            ir = 0
            read(ifi,err=99,end=99) kio,ir(5)
          else
            read(ifi,err=99,end=99) ir
          endif
          kio123 = mod(kio/10,1000)
          kio456 = mod(kio/10000,1000)
          if (.not. bittst(lio123,MNC)) ic1 = ir(1)
          if (.not. bittst(lio123,MNC)) ic2 = ir(2)
          if (.not. bittst(lio123,MLB)) liaxb = ir(4)
          if (.not. bittst(lio123,MNB)) nbas = ir(5)
          if (.not. bittst(lio123,MNT)) ntab = ir(6)

          call sanrg(.true.,ir(5),nbas,nbas,'iogfrs:','file''s nbas')
          if (bittst(lio123,MSC) .neqv. bittst(kio123,MSC))
     .      call rx("iogfrs:  spin coupling treatment does not match file")

          lerr = errh
          if (errh >= ipr) lerr = -lerr
          if (bittst(lio456,IPLMOD)) then
            if (bittst(lio456,ICLP)) call ipdump(clp(1,ic1),6,ifi)
          else
            call sanrg(.true.,ir(1),ic1,ic1,'iogfrs:','file''s starting cluster index')
            call sanrg(.true.,ir(2),ic2,ic2,'iogfrs:','file''s ending cluster index')
            call sanrg(.true.,ir(4),liaxb,liaxb,'iogfrs:','file''s liaxb')
            call sanrg(.true.,ir(6),ntab,ntab,'iogfrs:','file''s table size')
            if (rdarry(isw(bittst(lio456,ICLP)),bittst(kio456,ICLP),
     .        'clp',clp(1,ic1),nclp*ncp1,.false.,lerr,ifi) < 0) goto 99
            if (rdarry(isw(bittst(lio456,ICLSSL)),bittst(kio456,ICLSSL),
     .        'clssl',clssl,nbas,.false.,lerr,ifi) < 0) goto 99
            if (rdarry(isw(bittst(lio456,IIAX)),bittst(kio456,IIAX),
     .        'iaxg',iaxg(1,ic1),niax*ntab,.false.,lerr,ifi) < 0) goto 99
          endif
          if (.not. bittst(lio123,MPL)) then
            if (rdarry(isw(bittst(lio456,IPOS)),bittst(kio456,IPOS),
     .        'plat',plat,9,.true.,lerr,ifi) < 0) goto 99
            if (rdarry(isw(bittst(lio456,IPOS)),bittst(kio456,IPOS),
     .        'pos',pos,3*nbas,.true.,lerr,ifi) < 0) goto 99
          else
            if (rdarry(2,bittst(kio456,IPOS),
     .        'plat',plat,9,.true.,lerr,ifi) < 0) goto 99
            if (rdarry(2,bittst(kio456,IPOS),
     .        'pos',pos,3*nbas,.true.,lerr,ifi) < 0) goto 99
          endif
        endif

C   ... Read gf
        if (igf /= 0) then
          iogfrs = -2
          readin = '  GF'
          offnow = offg
          do  icl = ic1, ic2
            ndimbl = clp(3,icl)
            if (bittst(lio456,IDGF)) then
              ndimal = (clp(6,icl+1) - clp(6,icl)) / ndimb
            else
              ndimal = clp(4,icl)
            endif
            read(ifi,err=99,end=99) ndima,ndimb,cast,ii,igfl,zpl,qpl
            if (ndimb /= clp(3,icl)) call rx('iogfrs: file mismatch')
            kcplxl = cast-2
            if (bittst(lio123,MCD)) then
              call sanrg(.true.,ndima,ndimal,ndimal,'iogfrs:','file''s dimension of g')
              call sanrg(.true.,ndimb,ndimbl,ndimbl,'iogfrs:','file''s dimension of g')
            endif
            if (bittst(lio123,MZP)) then
              ltmp = tol > abs(zp(1)-zpl(1))+abs(zp(2)-zpl(2))
              ltmp = ltmp .and. tol >
     .          abs(qp(1)-qpl(1))+abs(qp(2)-qpl(2))+abs(qp(3)-qpl(3))
              ltmp = ltmp .and. igf == igfl
              if (.not. ltmp) goto 99
            endif
            if (bittst(lio123,MSC)) then
              ndima = ndima*2
              ndimb = ndimb*2
              if (offg /= 0) call rx('iogfrs check offg nc case')
            endif
            ii = ndima
            ndima = ndima*ndimb
            call dpsdmp(gf,1+2*offnow,2*(offnow+ndima),ifi)
C           call yprm('g',kcplxl+2,gf(1+2*offnow),ii*ndimb,ii,ii,ndimb)
C           Turn g into kcplx sought by caller
            if (kcplx /= kcplxl) then
              if (kcplxl == 0 .or. kcplx == 0) then
                call ztoyy(gf(1+2*offnow),ii,ndimb,ii,ndimb,kcplxl,kcplx)
C                call yprm('g',kcplx+2,gf(1+2*offnow),ii*ndimb,ii,ii,
C     .            ndimb)
              else
C               To convert kcplx 2 to 1, set ztoy lbak=1, else 0
                call ztoy(gf(1+2*offnow),ndima,ndima,1,2-kcplx)
              endif
C             call yprm('g',kcplx,gf(1+2*offnow),ii*ndimb,ii,ii,ndimb)
            endif
            offnow = offnow + ndima
          enddo
          ndims = offnow - offg
          call info5(75,0,0,' IOGFRS: read from file '//
     .      trim(locnam)//' %i elements %?;(n>1); from clusters %i...%i;;',
     .      ndims,ic1*ic2,ic1,ic2,0)
        endif
      endif

      iogfrs = 0
      return

C ... Error handling
   99 continue
      if (errh == 0) call rxs('IOGFRS: file mismatch, file ',fnam)
      call info0(errh,0,0,'IOGFRS (warning) failed to read '//readin//' from file')
      return

      end
