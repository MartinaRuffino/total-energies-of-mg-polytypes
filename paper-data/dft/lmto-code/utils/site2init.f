      subroutine fmain
      implicit none
C ... Dynamically allocated arrays
C      integer, parameter :: mxrec=5000,recl0=128,nspmx=200,nbmx=2000,ngmx=48,reclprs=500
C      integer, parameter :: nrmx=5001, nxi0=10, n0=10
C      integer,allocatable:: istab(:),iprm(:),lockc(:),ics(:),nrclas(:)
C      real(8),allocatable:: wk(:),zc(:),rmtc(:)
C      character,allocatable:: instr(:)*(recl0)
C ... Local parameters
      integer, parameter :: nspmx=200,nbmx=10000
      integer :: lcart=1
      integer ifin,j1
      double precision xx
      character fnam*128,header*128
      character(len=8),target :: slabl(nspmx)
      integer, target :: ips(nbmx)
      real(8),allocatable :: pos(:,:),posl(:,:)
      integer nspec,ib,nbas,i,j
      integer,parameter :: NULLI=-99999
C      procedure(logical) cmdopt,a2bin
      procedure(integer) iosite,fopng,nargf,fopnx,wordsw
      double precision alat,plat(3,3),qlat(3,3)

C ... Setup
      i = 1
      header = '-'
      do
        call getarf(i,fnam)
        if (fnam(1:1) /= '-' .or. fnam == ' ') exit

        j = wordsw(fnam,'','-header','= ',j1)
        if (j /= 0) then
          if (fnam(j1:j1) == ' ') then
            i = i+1
            call getarf(i,fnam)
            j1 = 0
          endif
          i = i+1
          header = fnam(j1+1:)
        endif

        j = wordsw(fnam,'','--wsitex','',j1) + wordsw(fnam,'','-wsitex','',j1)
        if (j /= 0) then
          lcart = 0
          i = i+1
        endif
      enddo

C     Get file name if one is specified
      if (nargf() == i) then
        fnam = 'sitein'
      else
        call getarf(i,fnam)
      endif

      if (fopnx(fnam,72,-1,-1) == 0)
     .  call rx('site2init : could not find site file '//trim(fnam))

C ... Open site file for reading.  init file goes to stdout
      ifin = fopng(trim(fnam),-1,1)

C     Load species table
      nspec = 0
      j1 = iosite(135002,3d0,0,trim(fnam),ifin,slabl,alat,plat,nbas,
     .  nspec,xx,xx,xx,xx,xx,xx,xx)
      if (j1 < 0) call rx('site2init failed to read site file')
      allocate(pos(3,nbas),posl(3,nbas))

      rewind ifin
C     Read pos, ips
      j1 = iosite(8002,3d0,0,trim(fnam),ifin,slabl,alat,plat,nbas,
     .  nspec,posl,xx,xx,xx,ips,xx,xx)

      if (lcart == 0) then
        call dinv33(plat,1,qlat,xx)
        call dgemm('T','N',3,nbas,3,1d0,qlat,3,posl,3,0d0,pos,3)
        print *, posl(:,2)
        print *, pos(:,2)
      else
        pos = posl
      endif

      if (header == '-') then
        call awrit2('init file constructed from file '//trim(fnam)//
     .    ' ... %i species and %i atoms',header,len(header),0,nspec,nbas)
      endif

      call info0(1,0,0,'HEADER '//trim(header))
      call info0(1,0,0,'LATTICE')
      call info2(1,0,0,'%% const a=%;7d',alat,2)
      call info0(1,0,0,'        ALAT={a}')
      call info5(1,0,0,'        PLAT=%3;12,7D%N%13f%3;12,7D%N%13f%3;12,7D',plat(1,1),plat(1,2),plat(1,3),4,5)
      call info0(1,0,0,'SITE')
      do  ib = 1, nbas
        call info2(1,0,0,'     ATOM='//trim(slabl(ips(ib)))//
     .    '%19p%?;n;POS=;X=; %3;12,7D',lcart,pos(1,ib))
      enddo

      end
