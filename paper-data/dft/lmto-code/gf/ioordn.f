      subroutine suiord(lio,ifi,ncl,s_ctrl,s_ordn)
C- Allocate memory for, and read O-N structural data
Cu   08 May 13 Eliminate s_array
      use structures
      implicit none
      integer lio,ifi,ncl
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_ordn)::  s_ordn
C ... Local parameters
      integer niax,nclp,psuor2,lio01,nttabb,nttabg,ival2
      parameter (niax=10,nclp=9)
      logical bittst
      double precision xx

      if (ifi < 0) call rxi('suiord encountered negative ifi:',ifi)
      lio01  = mod(lio,100)

      if (bittst(lio01,1)) then
        call ptr_ctrl(s_ctrl,8+1,'clp',(ncl+1)*nclp,0,0,xx)
        rewind ifi
        call ioordn(1,ifi,ncl,s_ctrl%clp,xx,xx)
      endif
      nttabb = psuor2(1,ncl,10,s_ctrl%clp)
      nttabg = ival2(s_ctrl%clp,nclp,2,ncl+1)

C     if (bittst(lio01,2))  call defi(oiaxg, niax*nttabg)
C     if (bittst(lio01,4))  call defi(omapgv,3*nttabb)

      if (bittst(lio01,2)) call rx('create ptr_ordn')
C      if (bittst(lio01,2))
C     .  call ptr_ordn(s_ordn,1,'iaxg',niax*nttabg,0,w)
      if (bittst(lio01,4)) call rx('create ptr_ordn')
C      if (bittst(lio01,4))
C     .  call ptr_ordn(s_ordn,1,'mapgv',3*nttabb,,0,w)

      call ioordn(lio,ifi,ncl,s_ctrl%clp,s_ordn%iaxg,s_ordn%mapgv)

      end

      subroutine ioordn(lio,ifi,ncl,clp,iaxg,mapgv)
C- Order-N parameters I/O
C ----------------------------------------------------------------------
Ci   lio  :1s+10s digit: bits set conditions for what is read:
Ci           1 read clp
Ci           2 read iaxg
Ci           4 read mapgv
Ci   ifi  :logical unit for I/O.  ifi>0 for read, ifi<0 for write
Cio File read/write
Cio  clp   :index and dimensioning information for cluster (suordn)
Cio         1: cumulative number of basis atoms in CL=1 ... this CL
Cio         2: last pair in iaxg table for this cluster
Cio         3: column (source) matrix dimension for this CL
Cio         4: row (field) matrix dimension for this CL
Cio         5: matrix dimension for direct zone of CL
Cio         6: offset to g for start of this cluster
Cio         7: last pair in iaxg table for direct zone, inner block
Cio         8: last pair in iaxg table for direct zone, outer block
Cio         9: last pair in iaxg table for which g is calculated
Cio  iaxg  :cluster neighbor table (paircl.f).  for each entry:
Cio        iaxg(1) = cluster index
Cio        iaxg(2) = basis id of site connected to this cluster
Cio        iaxg(3..5)multiples of plat to be added to site position
Cio        iaxg(6) = ? cluster used to extract g for this pair.
Cio        iaxg(7) = uniquely orders elements connected to this site
Cio        iaxg(8) = site ID, used in merging different lattices
Cio        iaxg(9) = hamiltonian dimension of field point
Cio        iaxg(10)= effective site index of this site. ? Not used?
Cio  clssl :site ib belongs to cluster clssl(ib) (clsset.f)
Cio  mapgv :mapgv(*,i) indicates where GF for pair i in site-equivalant
Cio        to iaxg table are extracted from
Cio        mapg(1,i) = offset to iaxg table for field point
Cio        mapg(2,i) = offset to iaxg table for source point
Cio        mapg(3,i) = indicator of how this element is generated.
Cio                 -1 from (1 + dP gref)^-1 gref
Cio                -10 free-electron GF
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lio,niax,nclp,ncl,ifi
      parameter (niax=10,nclp=9)
      integer clp(nclp,*),iaxg(niax,1),mapgv(3,*)
C ... Local parameters
      logical bittst
      integer lio01,psuor2,nttabb,jfi,nttabg

      lio01  = mod(lio,100)

C     Size of site-based iax table, corresponding to iaxg
      nttabb = psuor2(1,ncl,10,clp)
      nttabg = clp(2,ncl+1)

C --- File write ---
      if (ifi < 0) then
        jfi = -ifi
        rewind jfi
        call ywrm(1,' ',0,jfi,' ',clp,0,nclp,nclp,ncl+1)
        call ywrm(1,' ',0,jfi,' ',iaxg,0,niax,niax,nttabg)
        call ywrm(1,' ',0,jfi,' ',mapgv,0,3,3,nttabb)

C --- File read ---
      else
C       For now, assume everything is on disk
        call crdbf(bittst(lio01,1),.true.,'clp',clp,nclp,ncl+1,0,ifi)
        call crdbf(bittst(lio01,2),.true.,'iaxg',iaxg,niax,nttabg,0,ifi)
        call crdbf(bittst(lio01,4),.true.,'mapgv',mapgv,3,nttabb,0,ifi)
      endif

      end
