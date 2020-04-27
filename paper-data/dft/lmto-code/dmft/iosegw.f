      subroutine iosegw(mode,s_gwse,isp,nsp,nblst,iblst,ifi)
C- Read spectral functions in the se file format
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 read dimensioning parameters
Ci            s_gwse%nband,s_gwse%nqp,s_gwse%nomg,s_gwse%ommin,s_gwse%ommax,s_gwse%nblst
Ci         :1 Check that file dimensioning parameters agree with s_gwse
Ci         :2 Allocate and read s_gwse%iblst,s_gwse%qp,s_gwse%eigiq,s_gwse%sigiq
Ci         :  Note: iblst will be merged with passed iblst if passed nblst>0
Ci         :  Check that file dimensioning parameters agree with s_gwse
Ci         :3 Read data as in 2, but arrays are not allocated.
Ci         :  Also qp are checked against passed s_gwse%qp
Ci         :4 same as 0+2 in combination
Ci         :5 same as 1+3 in combination
Ci         :10s digit
Ci         :1 allocate and read sigwq
Ci   isp   : read data for spin channel isp
Ci   nsp   : number of spin channels in the system; needed for allocation
Ci   nblst : number of bands to read.  nblst=0 => read all bands
Ci   iblst : list of bands to read, if nblst>0
Ci   ifi   : file logical unit, but >0 for read, <0 for write
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  The se file format is documented here: /docs/input/data_format/#the-se-file
Cu Updates
Cu   20 Nov 16
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,isp,nsp,nblst,ifi,iblst(nblst)
C ... For structures
!       include 'structures.h'
      type(str_gwse)::    s_gwse
C ... Dynamically allocated arrays
      integer, allocatable :: nqlst(:)
      real(8),allocatable:: eigiq(:,:,:),sex(:,:,:),vxcl(:,:,:),sigiq(:,:,:),qp(:,:)
      complex(8),allocatable:: seq(:,:)
C ... Local parameters
      real(8), parameter:: tolq=1d-7
      logical ltmp
      integer fmode,i,lascii,mod0,mod1,i1,i2,i3,i4,iq,nomg,nfblst,nband,nqp
      real(8) :: ef
      procedure(logical) :: isanrg
      procedure(integer) :: iinear

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)

C --- Read header ---
      if (mod0 <= 5) then
        lascii = 20
        i = 3000 + 200 + lascii + 4+2
        fmode = 22
        rewind ifi
        call ioseh(i,fmode,i1,i2,nband,nqp,nomg,s_gwse%ommin,ef,s_gwse%ommax,nfblst,[0],i4,ifi)
        if (mod0 == 0 .or. mod0 == 4) then
          s_gwse%nband = nband
          s_gwse%nkp = nqp
          s_gwse%nomg = nomg
          s_gwse%nblst = nfblst
        elseif (mod0 == 1 .or. mod0 == 2 .or. mod0 == 5) then
          ltmp = isanrg(nband,s_gwse%nband,s_gwse%nband,'iosegw:','file''s nband',.true.)
          ltmp = isanrg(nqp,s_gwse%nkp,s_gwse%nkp,'iosegw:','file''s nqp',.true.)
          ltmp = isanrg(nomg,s_gwse%nomg,s_gwse%nomg,'iosegw:','file''s nomg',.true.)
          ltmp = isanrg(nfblst,s_gwse%nblst,s_gwse%nblst,'iosegw:','file''s nblst',.true.)
        endif

        if (mod0 < 2) return

C       Get list of bands
        if (mod0 == 2 .or. mod0 == 4) then
          allocate(s_gwse%iblst(2,nfblst))
          call iinit(s_gwse%iblst,2*nfblst)
        endif

        if (mod0 == 3) call rx('not ready mode 3')
        allocate(nqlst(0:i4))
        i = i + 4000
        rewind ifi
        call ioseh(i,fmode,i1,i2,nband,nqp,nomg,s_gwse%ommin,s_gwse%ommax,nfblst,s_gwse%iblst(2,1:nfblst),nqlst,ifi)
        deallocate(nqlst)
C       Collect list of bands to read in s_gwse%iblst(1,:)
        if (nblst > 0) then
          i3 = 0
          do  i1 = 1, nfblst
            i2 = s_gwse%iblst(2,i1)
            i = iinear(nblst,i2,iblst,1)
            if (iblst(i) /= i2) i2 = 0
            if (i2 == 0) cycle
            i3 = i3+1
            s_gwse%iblst(1,i3) = i2
          enddo
          s_gwse%nblst = i3
        else
          s_gwse%iblst(1,:) = s_gwse%iblst(2,:)
        endif
        call info2(30,0,0,' IOSEGW: read data for %i bands out of %i from se file',s_gwse%nblst,nfblst)
        if (s_gwse%nblst == 0) call rx('no bands to read')

C       Read qp, sigiq, eigiq into temporary arrays
        if (mod0 == 2 .or. mod0 == 4) then
          allocate(s_gwse%qp(3,nqp),s_gwse%sgq(s_gwse%nblst,nqp,nsp),
     .      s_gwse%eig(s_gwse%nblst,nqp,nsp),s_gwse%sxq(s_gwse%nblst,nqp,nsp))
        endif
        allocate(qp(3,nqp),sigiq(nfblst,nqp,nsp),eigiq(nfblst,nqp,nsp),sex(nfblst,nqp,nsp),vxcl(nfblst,nqp,nsp))
        call iose2(lascii,fmode,1,nfblst,nqp,nomg,qp,eigiq,sigiq,[0d0],sex,vxcl,ifi)
        call dcopy(3*nqp,qp,1,s_gwse%qp,1)

        i3 = 0
        do  i1 = 1, nfblst
          i2 = s_gwse%iblst(2,i1)
          i = iinear(s_gwse%nblst,i2,s_gwse%iblst,2)
          if (s_gwse%iblst(1,i) /= i2) i2 = 0
          if (i2 == 0) cycle
          i3 = i3+1
          s_gwse%eig(i3,:,isp) = eigiq(i1,:,isp)
          s_gwse%sgq(i3,:,isp) = sigiq(i1,:,isp)
          s_gwse%sxq(i3,:,isp) = sex(i1,:,isp)
        enddo
        deallocate(qp,sigiq,eigiq)
      endif                     ! end of header input

C --- Read frequency dependent self-energy ---
      if (mod1 == 0) return
      nfblst = size(s_gwse%iblst)/2
      nomg = s_gwse%nomg
      nqp = s_gwse%nkp
      allocate(seq(nomg,nfblst))
      allocate(s_gwse%sigwq(nomg,s_gwse%nblst,nqp,nsp))
      do  iq = 1, nqp
        call dfdump(seq,2*nomg*nfblst,ifi)
        i3 = 0
        do  i1 = 1, nfblst
          i2 = s_gwse%iblst(2,i1)
          i = iinear(s_gwse%nblst,i2,s_gwse%iblst,2)
          if (s_gwse%iblst(1,i) /= i2) i2 = 0
          if (i2 == 0) cycle
          i3 = i3+1
          call dcopy(2*nomg,seq(1,i1),1,s_gwse%sigwq(1,i3,iq,isp),1)
        enddo
      enddo
      deallocate(seq)

      end
