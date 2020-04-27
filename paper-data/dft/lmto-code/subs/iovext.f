      integer function iovxta(mode,nc,nsp,nl,lmx,z,ics,pp,lrqd,ifi)
C-
C ----------------------------------------------------------------
Ci Inputs
Ci   mode  :1 read 2nd gen C, delta
Ci   nc    :number of classes
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   z     :nuclear charge
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci         :Set to zero to suppress information for core hole
Ci   ifi   :file index
Ci         :ifi>0 for READ
Ci   lrqd  :if T, exit with error when
Co Outputs
Co   qc:    core electronic charge
Co   qt:    total charge (nuclear + electronic) within sphere
Co   dq:    difference between spin up and spin down charge
Cl Local variables
Cl    icnt :read from disk: indicates what parameters are written
Cl         :to disk.
Cl         :1 shifts in 2nd gen ASA C, delta
Cl    ifmt :read from disk: format of file contents
Cl         :0 PPAR written as:  class label l isp  shft-C  shft-Delta
Cr Remarks
Cr   This is the input routine for adjustments to ASA parmeters.
Cu Updates
Cu   27 Aug 09 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lrqd
      integer mode,ifi,nl,nsp,nc,lmx(nc),ics(nc)
      double precision pp(6,nl,nsp,*),z(nc)
C ... Local parameters
      logical scat
      integer i,ifmt,nxtlin,partok,read_token,ic,l,isp,
     .  icnt,n
      double precision c,delta,xx,xx2
      character s*120, class*8

C --- File read ---
      if (ifi > 0) then

C     Check conditions for file to be a valid input file
      iovxta = -1
      rewind ifi
      if (nxtlin(ifi,s) < 0) goto 10
      iovxta = -2
      if (s(1:1) /= '%') goto 10
      if (read_token(s,'contents=',' ',1,2,1,icnt) <= 0) goto 10
      if (read_token(s,'format=',' ',1,2,1,ifmt) <= 0) goto 10

      call info5(20,1,0,
     .  ' IOVEXT: reading external pot. parameters file, '//
     .  'format %i.'//
     .  '%N File Contains: *'//
     .  '%?#(n%2)#%b PPAR,##%-1j'//
     .  '%?#(flor(n/2)%2)#%b XYZ,##%-1j'//
     .  '%a%b%j'//
     .  '%N Reading: *'//
     .  '%?#(n%2)#%b PPAR,##%-1j'//
     .  '%?#(flor(n/2)%2)#%b XYZ,##%-1j'//
     .  '%a%b',
     .  ifmt,icnt,mode,0,0)


C --- Read 2nd generation ppar shifts ---
      iovxta = -5
      if (mod(mode,2) /= 0) then
      if (.not. scat(ifi,'PPAR:',':',.true.)) goto 10

      n = 0
      do while (nxtlin(ifi,s) == 0)

        if (ifmt == 0) then
          iovxta = -6
          read(s,*,iostat=i) ic,class,l,isp,c,delta
          if (i /= 0) goto 10
          pp(2,l+1,isp,ic) = pp(2,l+1,isp,ic) + c
          xx = sign(1d0,pp(3,l+1,isp,ic))
          xx2 = dsqrt(pp(3,l+1,isp,ic)**2 + delta)
          pp(3,l+1,isp,ic) = xx*xx2
        else
          iovxta = -4; goto 10
        endif

      enddo

      endif

C --- Read contents for 2's bit mode ---
      if (mod(mode/2,2) /= 0) then
      endif

      return
      endif  ! File read

      goto 10

      return
   10 continue
      if (.not. lrqd) return
      select case( iovxta )
        case(-1) ; s = 'it has no records'
        case(-2) ; s =
     .      '1st line requires header format "%% contents=# format=#"'
        case(-4) ; s = 'format style not recognized'
        case(-5) ; s = 'no category PPAR found'
        case(-6) ;
        case default ; s = '%a%b '
      end select
      call rx('iovxta: failed to read file:  '//trim(s))

      end

      integer function read_token(instr,token,sep0,i1,cast,n,res)
C- Parse for algebraic expressions following a token in a string
C ----------------------------------------------------------------------
Ci Inputs
Ci   instr : string to parse
Ci   token : token to flag start parse region
Ci    sep0 :the character immediately to the left of the token
Ci         :must be one of the characters in sep0
Ci    i1   :start search at i1 (i1=1 for 1st char in string)
Ci    cast :0=logical, 2=int, 3=real, 4=double
Ci       n :number of elements following token
Co Outputs
Co     res :results from parsing string
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Aug 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) token,sep0,instr
      integer i1,n,cast
      integer res(*)
C ... Local parameters
      integer lstr,i2,rdtk2,ii,a2vec,iarr(n)

      read_token = -999
      lstr = len(instr)
      i2 = rdtk2(instr,lstr,token,sep0,.false.,i1)
      if (i2 < 0) return

      ii = 0
      read_token =
     .  a2vec(instr(i2:),lstr-i2+1,ii,cast,', ',2,-3,n,iarr,res)

      end

