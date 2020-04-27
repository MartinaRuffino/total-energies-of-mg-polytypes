      integer function lmfitmr1(mode,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,
     .  z,iprmb,ldim,ncof,nvar,ivcof,pp,ftcof)
C- Setup for Levenberg-Marquardt fit of ASA params to given bands
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec class clabel
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: iq1 iq2
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   mode  :1s digit controls what to return
Ci         :0 return ncof,nvar
Ci         :1 return ivcof; copy ASA parameters to ftcof,
Ci         :2 Copy ftcof to ASA parameters; printout change this iter
Ci         :3 Write difference between ftcof and ASA parameters to disk
Ci         :4 Read adjustment to ASA parameters from file; add to pp's
Ci         :5 same as 0, but constraints taken from vext file
Ci         :6 same as 1, but constraints taken from vext file
Ci         :10s digit controls what to vary
Ci         :1 ASA C and Delta
Ci         :2 ASA C+enu and Delta
Ci   nl    :(global maximum l) + 1
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Cio Inputs/Outputs
Cio   ncof  :number of coefficients that can potentially vary; see Remarks
Cio         :1s digit mode 0 or 5: output; otherwise input
Cio         :Not used if 1s digit mode = 4
Cio   nvar  :number of (independent) parameters to actually vary; see Remarks
Cio         :1s digit mode 0 or 5: output; otherwise input
Cio         :Not used if 1s digit mode = 4
Cio   pp    :Potential parameters
Cio         :mode  Usage
Cio         :0     Not used
Cio         :1     pp is input, copied to ftcof
Cio         :2     pp is output, copied from ftcof
Cio         :3     difference between ftcof and pp is written to disk
Cio         :4     file adjustment to pp's read from disk, added to pps
Cio   ivcof :Contains information about coefficients that may be varied.
Cio         :Not used if 1s digit mode = 4
Cio         :See Remarks
Cio         : ivcof(ivar,1): vector of indices to coefficients are varied
Cio         : ivcof(1:ncof,2) : contains information about constraints on
Cio         :              : each coefficient.
Cio         :  ivcof(i,2)  Value signifies:
Cio         :    0         coefficient i is free to vary
Cio         :   <0         coefficient i is free to vary, but other
Cio         :              coefficients are linked to it
Cio         :-NULLI        coefficient i is frozen; otherwise:
Cio         :   >0         change in coefficient i is linked to
Cio         :              coefficient ivcof(i,2)
Co Outputs
Co    ftcof :(mode 1 or 6) vector of coefficients that may be varied
Cl Local variables
Cl  nspf    :number of spins in vext file
Cr Remarks
Cr   There are ncof coefficients all told.  They are one of three types:
Cr   - coefficient is free to vary
Cr     nvar is the number of such coefficients.
Cr   - coefficient is frozen
Cr   - change in coefficient is linked to change in another coefficient.
Cu Updates
Cu   26 Nov 17 Added sp2cls tag
Cu   10 Nov 11 Begin migration to f90 structures
Cu   17 Nov 10 L-M fitting with C,enu in tandem
Cu   27 Aug 10 (H Donfack) Made linking constraints work
Cu   09 Aug 10 vext file has spin info, to enable spin splitting
Cu   01 Jul 10 Implement for spin-pol case
Cu   12 May 10 Allow format styles 1,2
Cu   25 Mar 10 New 1s digit modes 5,6
Cu   10 Aug 09 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nl,nbasp,nclasp,nsp,ldim,ncof,nvar
      integer iprmb(nl**2*nbasp),ivcof(ncof,2)
      double precision ftcof(ncof),pp(6,nl,nsp,nclasp),z(nclasp)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamical arrays
      integer, allocatable :: ibc(:)
      integer, allocatable :: clist(:), llink(:)
      real(8),allocatable:: delC(:,:,:),delD(:,:,:)
C ... Local parameters
      logical scat,lrqd
      integer n0,stdo,ipr
      parameter(n0=10,lrqd=.true.)
      integer i,k,ib,l,is,ic,jc,ifrzc(n0,2),ifrzd(n0,2),lmr,
     .  ncofl,nvarl,isp,ii(10),NULLI
      parameter (NULLI =-99999)
      integer mode0,mode1,ifi,fopna,nlist
      character clabl*8,clabl2*8,frzc(0:2),s*120
      double precision xx,xx2,delCl
      integer ifmt,nxtlin,read_token,icnt,lstyle
      integer mpipid,procid,master,nspf,sp2cls
      procedure(integer) :: nglob,iinear,lmfitmr4
      data frzc /' ','*','+'/

C ... Setup
      stdo = nglob('stdo')
      call getpr(ipr)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      lmfitmr1 = 0
      procid = mpipid(1)
      master = 0

      ncofl = 0
      nvarl = 0
      lmr = 0
      sp2cls = 0

C     Check conditions for file to be a valid input file
      if (mode0 == 4 .or. mode0 == 5 .or. mode0 == 6) then
        lmfitmr1 = -1
        ifi = fopna('vext',-1,0)
        rewind ifi
        if (nxtlin(ifi,s) < 0) goto 10
        lmfitmr1 = -2
        if (s(1:1) /= '%') goto 10
        if (read_token(s,'contents=',' ',1,2,1,icnt) <= 0) goto 10
        if (read_token(s,'format=',' ',1,2,1,ifmt) <= 0) goto 10
        if (read_token(s,'nsp=',' ',1,2,1,nspf) <= 0) goto 10
        if (read_token(s,'sp2cls=',' ',1,2,1,sp2cls) <= 0) sp2cls = 0
      endif

      if (ipr >= 30 .and.
     .  (mode0 <= 1 .or. mode0 == 5 .or. mode0 == 6)) then
        write(stdo,333)
  333   format(' LMFITMR: parameters to include in fit:'/
     .    '   class',8x,'l isp    block   vary(C) vary(D)  ivar',
     .    5x,'C',10x,'Delta')
      else if (ipr >= 30 .and. mode0 <= 3) then
        write(stdo,335)
  335   format(' LMFITMR: changes in parameters'/
     .    '   class',8x,'l     old C     old Delta',7x,
     .    'new C     new Delta',13x,'Change')
      else if (mode0 == 4) then

C       Check conditions for file to be a valid input file
        lmfitmr1 = -1
        ifi = fopna('vext',-1,0)
        rewind ifi
        if (nxtlin(ifi,s) < 0) goto 10
        lmfitmr1 = -2
        if (s(1:1) /= '%') goto 10
        if (read_token(s,'contents=',' ',1,2,1,icnt) <= 0) goto 10
        if (read_token(s,'format=',' ',1,2,1,ifmt) <= 0) goto 10
        if (read_token(s,'nsp=',' ',1,2,1,nspf) <= 0) goto 10

        call info5(20,1,0,
     .  ' LMFITMR: reading external pot. parameters file, format %i'//
     .  '%?#(n==1)# (spec->class)##.'//
     .  '%N File Contains: *'//
     .  '%?#(n==1|n==2)#%b PPAR,##%-1j'//
     .  '%a%b%j'//
     .  '%N Reading: *'//
     .  '%?#(n==1|n==2)#%b PPAR,##%-1j'//
     .  '%a%b',
     .  ifmt,sp2cls,icnt,mode1,5)

        if (nspf > nsp) call info0(30,0,0,
     .    ' (warning) 2nd spin in vext file ignored')
        if (nspf < nsp) call info0(30,0,0,
     .    ' (warning) vext file not spin pol .. splitting spins')
        if (icnt /= mode1)
     .    call rx('LMFITMR: contest mismatch file vext with input')

C   --- mode1<=2: 2nd generation ppar shifts ---
        lmfitmr1 = -5
        if (mode1 == 1 .or. mode1 == 2) then
          if (.not. scat(ifi,'PPAR:',':',.true.)) goto 10

          allocate(clist(nbasp),delC(nl,nsp,nclasp),delD(nl,nsp,nclasp))
          call dpzero(delC,nl*nsp*nclasp)
          call dpzero(delD,nl*nsp*nclasp)
          do while (nxtlin(ifi,s) == 0)
            if (lmfitmr4(1,ifmt,s,ic,clabl2,l,isp,xx,xx2,ii(2),ii(3)) < 0) goto 10
            if (l >= nl) then
              call info2(40,0,0,' LMFITMR (warning) skipping (l=%i, '//
     .          'ic=%i) in vext file',l,ic)
              cycle
            endif
            nlist = 1
            if (ifmt == 2) then
              lstyle = 2
              call word(s,1,ii(1),ii(2))
              call slist(lstyle,s(ii(1):ii(2)+1),clabl,z,nclasp,nlist,clist)
            elseif (sp2cls == 1) then
              nlist = 0
              do  jc = ic, nclasp
                if (jc == ic .or. s_ctrl%ics(jc) == ic) then
                  nlist = nlist+1
                  clist(nlist) = jc
                endif
              enddo
            endif
C           print *, nlist,clist(1:nlist)
            do  jc = 1, nlist
              if (ifmt == 2 .or. sp2cls == 1) ic = clist(jc)
              delC(l+1,isp,ic) = xx
              delD(l+1,isp,ic) = xx2
              if (nsp == 2 .and. nspf == 1) then
                delC(l+1,3-isp,ic) = xx
                delD(l+1,3-isp,ic) = xx2
              endif
            enddo
          enddo

          if (ipr >= 40) write(stdo,235)
  235     format(' LMFITMR: parameter shifts read from file:'/
     .    '  ic  l isp     del-C',7x,'del-D')

          do  ic = 1, nclasp
            do  l = 0, nl-1
            do  isp = 1, nsp

              if (ipr >= 40)
     .        write(stdo,331) ic,l,isp,delC(l+1,isp,ic),deld(l+1,isp,ic)
  331         format(i4,2i3,2f12.6)

              pp(2,l+1,isp,ic) = pp(2,l+1,isp,ic) + delC(l+1,isp,ic)
              if (mode1 == 2) then
              pp(1,l+1,isp,ic) = pp(1,l+1,isp,ic) + delC(l+1,isp,ic)
              endif
              delC(l+1,isp,ic) = sign(1d0,pp(3,l+1,isp,ic))
              delD(l+1,isp,ic) = dsqrt(pp(3,l+1,isp,ic)**2 +
     .                           delD(l+1,isp,ic))
              pp(3,l+1,isp,ic) = delC(l+1,isp,ic)*delD(l+1,isp,ic)
            enddo
            enddo
          enddo
          deallocate(clist,delC,delD)

C   --- Alternate variables : nothing yet ---
        else
          call rxi(' LMFITMR: mode not recognized, mode ',mode1)
        endif

        return
      endif

      if (mode0 == 3) then
        ifi = fopna('vext0',-1,0)
        rewind ifi
        write(ifi,360) mode1,nsp
  360   format('% contents=',i1,' format=0 nsp=',i1)
        write(ifi,361)
  361   format('PPAR: ASA2')
        write(ifi,358)
  358   format('# class',9x,'l isp   shft C    shft Delta  frz')
      endif

C --- Loop over sites (skip classes already touched inside loop) ---
      allocate(ibc(nbasp))
      call iinit(ibc,nbasp)
      allocate(llink(2*nl*nsp*nclasp))
      call iinit(llink,nl*nsp*nclasp*2)
      do  ib = 1, nbasp

        is = s_site(ib)%spec
        ic = s_site(ib)%class
        clabl = s_site(ib)%clabel

C   --- mode 1 : vary C and Delta ---
        if (mode1 == 1 .or. mode1 == 2) then

C   ... Skip classes already touched
        if (ibc(ic) /= 0) then
          lmr = lmr + nl**2
        else

C       Read constraints from file
        if (mode0 == 5 .or. mode0 == 6) then

          call ivset(ifrzc,1,n0*2,0)
          call ivset(ifrzd,1,n0*2,0)
          rewind ifi
          if (nxtlin(ifi,s) < 0) goto 10
          if (nxtlin(ifi,s) < 0) goto 10
          if (ifmt < 0 .or. ifmt > 2) then
            lmfitmr1 = -4; goto 10
          endif
          allocate(clist(nbasp))
          do while (nxtlin(ifi,s) == 0)
            if (lmfitmr4(2,ifmt,s,ii(1),clabl2,l,isp,xx,xx2,
     .        ii(2),ii(3)) < 0) goto 10
            nlist = 1
            if (ifmt == 2) then
              lstyle = 2
              call word(s,1,ii(4),ii(5))
              call slist(lstyle,s(ii(4):ii(5)+1),clabl,z,nclasp,nlist,clist)
C             print *, nlist,clist(1:nlist)
            elseif (sp2cls == 1) then
              nlist = 0
              do  jc = ic, nclasp
                if (jc == ic .or. s_ctrl%ics(jc) == ic) then
                  nlist = nlist+1
                  clist(nlist) = jc
                endif
              enddo
            endif
            if (ifmt == 2 .or. sp2cls == 1) then
              ii(1) = 0
              if (nlist > 0) then
                call hunti(clist,nlist,ic,0,jc)
                if (jc < nlist) ii(1) = clist(jc+1)
              endif
            endif

            if (ii(1) == ic) then
              ifrzc(l+1,isp) = ii(2)
              ifrzd(l+1,isp) = ii(3)
            endif

          enddo
          deallocate(clist)
        else
          ifrzc(1:n0,1) = s_spec(is)%iq1(1:n0)
          call icopy(n0,ifrzc,1,ifrzc(1,2),1)
          ifrzd(1:n0,1) = s_spec(is)%iq2
          call icopy(n0,ifrzd,1,ifrzd(1,2),1)
        endif

        do  l = 0, nl-1
          do  isp = 1, nsp

C           Indices to those parameters to be varied will be set later
            if (mode0 == 1 .or. mode0 == 6) then
              ivcof(ncofl+1,1) = 0
              ivcof(ncofl+2,1) = 0
              ftcof(ncofl+1) = 0
              ftcof(ncofl+2) = 0
            endif

C           Case Rl channel is in lower block
            if (iprmb(lmr+1) < ldim) then

              ii(1) = l
              ii(2) = isp
              ii(3) = ic

C             mode=1: Copy pp's to ftcof
              if (mode0 == 1 .or. mode0 == 6) then
                ftcof(ncofl+1) = pp(2,l+1,isp,ic)
C               pp(3) holds sqrt(Delta); copy Delta, preserving sign
                xx = sign(1d0,pp(3,l+1,isp,ic))
                ftcof(ncofl+2) = xx*pp(3,l+1,isp,ic)**2
C             mode=2: Copy ftcof to pp's; mode=3: write ftcof-pp's to disk
              elseif (mode0 == 2 .or. mode0 == 3) then
C               Recover ifrzc,ifrzd.  ii(5:6) = (-1,0,1)
                do  i = 1, 2
                  k = ivcof(ncofl+i,2)
                  ii(4+i) = min(max(k,0),1)
                  if (k == -NULLI) then
                    k = 1
                  elseif (k > 0) then
                    k = ivcof(k,2)
                    ii(4+i) = 2
                  endif
C                 print *, ic,l,k,i,ii(4+i)
                  if (i == 1) ifrzc(l+1,isp) = k
                  if (i == 2) ifrzd(l+1,isp) = k
                enddo
                xx = sign(1d0,pp(3,l+1,isp,ic))*pp(3,l+1,isp,ic)**2
                if (ipr >= 30) then
                  write(stdo,336) ic,clabl,l,
     .              pp(2,l+1,isp,ic),frzc(ii(5)),xx,frzc(ii(6)),
     .              ftcof(ncofl+1),ftcof(ncofl+2),
     .              ftcof(ncofl+1)-pp(2,l+1,isp,ic), ftcof(ncofl+2)-xx
  336             format(i4,2x,a,i3,f12.6,a1,f11.6,a1,1x,2f12.6,2x,
     .              2f12.6)
                endif
                if (mode0 == 2) then
                  delCl = ftcof(ncofl+1) - pp(2,l+1,isp,ic)
                  pp(2,l+1,isp,ic) = ftcof(ncofl+1)
                  if (mode1 == 2) then
                    pp(1,l+1,isp,ic) = pp(1,l+1,isp,ic) + delCl
                  endif
                  xx = sign(1d0,ftcof(ncofl+2))
                  pp(3,l+1,isp,ic) = xx*dsqrt(abs(ftcof(ncofl+2)))
                else
                  write(ifi,357) ic,clabl,l,isp,
     .              ftcof(ncofl+1)-pp(2,l+1,isp,ic),ftcof(ncofl+2)-xx,
     .              ifrzc(l+1,isp),ifrzd(l+1,isp)
  357             format(i4,2x,a,2i3,2f12.6,2i3)
                endif
              endif

C         ... Generate internal linked coefficients table
              llink(ncofl+1) = ifrzc(l+1,isp)
              llink(ncofl+2) = ifrzd(l+1,isp)
              do  i = 1, 2
                if (llink(ncofl+i) > 0) then
                  llink(ncofl+i) = -NULLI
                elseif (llink(ncofl+i) < 0) then
                  k = (iinear(ncofl/2+1,llink(ncofl+i),llink(i),2)-1)*2 + i
                  if (k /= ncofl+i) then
                    llink(ncofl+i) = k
                  endif
                endif
              enddo

C         ... Specific to modes 0,1,5,6
              if (mode0 <= 1 .or. mode0 == 5 .or. mode0 == 6) then

C           ... Store ivcof(:,2)
                if (mode0 == 1 .or. mode0 == 6) then
                  ivcof(ncofl+1,2) = llink(ncofl+1)
                  ivcof(ncofl+2,2) = llink(ncofl+2)
C                  print "(10i6)",ib,l,ncofl,ivcof(ncofl+1,2),
C     .              ivcof(ncofl+2,2),ifrzc(l+1,isp),ifrzd(l+1,isp)
                endif

                if (ipr >= 30) then
                if (llink(ncofl+1) <= 0 .or.
     .              llink(ncofl+2) <= 0) then
                  write(stdo,334) ic,clabl,l,isp,' low',ifrzc(l+1,isp),
     .              ifrzd(l+1,isp),nvarl+1,ftcof(ncofl+1),ftcof(ncofl+2)
                else
                  write(stdo,334) ic,clabl,l,isp,' low',ifrzc(l+1,isp),
     .            ifrzd(l+1,isp),0,ftcof(ncofl+1),ftcof(ncofl+2)
                endif
                endif

C               Vary C parameter
                if (llink(ncofl+1) <= 0) then
                  nvarl = nvarl + 1
                  if (mode0 == 1 .or. mode0 == 6) then
                    ivcof(nvarl,1) = ncofl+1
                  endif
                endif
C               Vary D parameter
                if (llink(ncofl+2) <= 0) then
                  nvarl = nvarl + 1
                  if (mode0 == 1 .or. mode0 == 6) then
                    ivcof(nvarl,1) = ncofl+2
                  endif
                endif
              endif
C           Case Rl channel is not in lower block
            else
              if (ipr >= 30 .and. (mode0 == 2 .or. mode0 == 3)) then
                write(stdo,332) ic,clabl,l,'high'
              elseif (ipr >= 30) then
                write(stdo,334) ic,clabl,l,isp,'high'
              endif
            endif
  332       format(i4,2x,a,i3,5x,a,L8,L8,i8,2f12.6)
  334       format(i4,2x,a,2i3,5x,a,i8,i8,i8,2f12.6)
  339       format(i4,2x,a,2i3,5x,a,a16,i8,2f12.6)
            ncofl = ncofl + 2

          enddo ! Loop over spin
          lmr = lmr + 2*l+1
        enddo ! Loop over l
        endif ! done with this class

C   --- Error if mode not recognized ---
        else
          call rxi(' LMFITMR: mode not recognized, mode ',mode1)
        endif

        ibc(ic) = ib ! Flags class has been treated

      enddo

C --- 1s digit mode zero: assign ncof,nvar ---
      if (mode0 == 0 .or. mode0 == 5) then
        ncof = ncofl
        nvar = nvarl
      endif
      if (mode0 == 3) then
        call fclose(ifi)
      endif

C      if (mode0 == 1 .or. mode0 == 6) then
C        call icopy(ncof,llink,1,ivcof(1,2),1)
C      endif

      deallocate(llink)
      deallocate(ibc)

C --- Error exit ---
      return
   10 continue
      if (.not. lrqd) return
      select case( lmfitmr1 )
        case(-1) ; s = 'it has no records'
        case(-2) ; s =
     .  '1st line requires header format "%% contents=# format=# nsp=#"'
        case(-4) ; s = 'format style not recognized'
        case(-5) ; s = 'no category PPAR found'
        case(-6) ;
        case default ; s = '%a%b '
      end select
      call rx('lmfitmr1: failed to read file:  '//trim(s))

      end
      subroutine lmfitmr2(mode,nlist,ii,nvar,nl,nsp,ivfit,pp)
C- Set or change current parameter in list of parameters to vary
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :-1 Information mode
Ci         :   return eps in pp(1,1,1,1)
Ci         :   nl,nsp should be > 0
Ci         :   ii,nvar,ivfit are not used.
Ci         :0 Initialization mode: set ivfit(1,2) and return
Ci         :1 Return in ii information set by init step for one parameter
Ci         :2 Increment current list of parameters by epsilon
Ci         :  (parameter details specified by ivfit(1:nlist,:))
Ci         :4 Decrement current list of parameters by epsilon
Ci         :6 Increment current list of parameters by epsilon
Ci         :8 Decrement current list of parameters by epsilon
Ci   nlist :List of current parameters, or index to current parameter
Ci         :mode 0,1: nlist points to the entry in ivfit that
Ci         :is to be modified or used
Ci         :mode >1  nlist is the number of parameters whose
Ci         :         valued is to be adjusted
Ci   ii    :If 1s digit mode = 0 or 1: vector of indices that specify
Ci         :detailed information about parameter type.
Ci         :Case ivfit(:,1) = 1: parameter is potential parameter
Ci         :  ii(1) = l quantum number
Ci         :  ii(2) = spin index isp
Ci         :  ii(3) = class index ic
Ci         :If 1s digit mode = 2: not used
Ci   nvar  :Leading dimension of ivfit
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cio Inputs/Outputs
Cio  ivfit :indices to coefficients that are varied
Cio        :ivfit(:,1) flags type of coefficient to vary (ALWAYS input)
Cio        :ivfit(:,2) flags which kind of coefficient within a type.
Cio        :           Value should not be altered by calling program.
Cio        :           Set internally if 1s digit mode = 0
Cio        :           Used as input when 1s digit mode is nonzero.
Cio        :Case ivfit(:,1) = 1..6: parameter is potential parameter
Cio        : ivfit(1)   ivfit(2)                Parameter
Cio        :    2     1l+10isp+100ic        C (for given l,isp,ic)
Cio        :    3     1l+10isp+100ic        Delta (for given l,isp,ic)
Co Outputs
Ci   pp    :potential parameters, which may be modified by this routine
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Aug 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nlist,nl,nsp,nvar,ii(4),ivfit(nvar,2)
      double precision pp(6,nl,nsp,*)
C ... Local parameters
      integer l,isp,ic,it,mode0,mode1,ilist
      double precision eps,shft,xx
      parameter (eps=1d-4)

C ... Information mode
      if (mode == -1) then
        pp(1,1,1,1) = eps
        return
      endif

      it = ivfit(nlist,1)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      if (mode0 == 0) then
        if (it <= 6) then
          l = ii(1)
          isp = ii(2)
          ic = ii(3)
          ivfit(nlist,2) = l+10*isp+100*ic
        else
          call rxi('no rule for ivfit(1)=',ivfit(nlist,1))
        endif
      elseif (mode0 == 1) then
        if (it <= 6) then
          l   = mod(ivfit(nlist,2),10)
          isp = mod(ivfit(nlist,2)/10,10)
          ic  = mod(ivfit(nlist,2)/100,1000)
          ii(1) = l
          ii(2) = isp
          ii(3) = ic
        else
          call rxi('no rule for ivfit(1)=',ivfit(nlist,1))
        endif
      elseif (mode0/2 /= 0) then
        it = ivfit(1,1)
        do ilist=1, nlist
        shft = eps
        if (mode0/2 >= 3) shft = 2*eps
        if (mod(mode0/2,2) == 0) shft = -shft
        if (it <= 6) then
          l   = mod(ivfit(ilist,2),10)
          isp = mod(ivfit(ilist,2)/10,10)
          ic  = mod(ivfit(ilist,2)/100,1000)
C         pp(3) holds sqrt(Delta)
          if (it == 3) then
            xx = sign(1d0,pp(it,l+1,isp,ic))
            pp(it,l+1,isp,ic) = xx*
     .        dsqrt(pp(it,l+1,isp,ic)**2 + shft)
          else
            pp(it,l+1,isp,ic) = pp(it,l+1,isp,ic) + shft
            if (it == 2 .and. mode1 == 2) then
              pp(1,l+1,isp,ic) = pp(1,l+1,isp,ic) + shft
            endif
          endif
        else
          call rxi('no rule for ivfit(1)=',ivfit(ilist,1))
        endif
        enddo
      endif

      end

      subroutine lmfitmr3(iv,ldim,dH,dS,eval,ndg,gradE)
C- Make gradient(eval) from dH, dS, eval
C ----------------------------------------------------------------------
Ci Inputs
Ci   iv    :Index to current parameter
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   dH    :Change in Hamiltonian for parameter iv
Ci   dS    :Change in Overlap for parameter iv
Ci   eval  :Eigenvalues
Ci   ndg   :Leading dimension of gradE
Co Outputs
Co   gradE :gradient eval for parameter iv
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   11 Aug 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iv,ldim,ndg
      double precision dH(ldim,ldim),dS(ldim,ldim),eval(ldim)
      double precision gradE(ndg,ldim)
C ... Local parameters
      double precision eps
      integer i

C     Get epsilon used for finite difference numerical differentiation
      call lmfitmr2(-1,1,[1],1,1,1,1,eps)

      do  i = 1, ldim
        gradE(iv,i) = (dH(i,i) - eval(i)*dS(i,i))/eps
      enddo
C      print 333, grade(iv,:)
C  333 format(6f12.6)
      end

      subroutine lmrqit(mode,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,z,iprmb,
     .  ldim,pp,npfit,ncof,nvfit,ivcof,dat,y,dy,sig,wk,cof,alsc,
     .  alp,alam,cov,chi,iter)
C- Next Levenberg-Marquardt iteration, interactive
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit controls what to return
Ci         :2 Copy ftcof to ASA parameters
Ci         :10s digit controls what to vary
Ci         :1 ASA C and Delta
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec class clabel
Ci     Stored: *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: iq1 iq2
Ci     Stored: *
Ci     Passed to: *
Ci   nl    :(global maximum l) + 1
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   pp    :potential parameters (atomsr.f)
Ci   npfit :number of points to fit
Ci   ncof  :total number of coefficients in fitting function;
Ci         :leading dimension of alp and cov
Ci   nvfit :number of coefficients to vary out of ncof parameters
Ci         :defining function.  The rest are frozen or constrained;
Ci         :see Remarks.
Ci         :ivcof below points to list of 1:nvfit parameters to be varied.
Ci         :Otherwise, mrqstp will check for consistency in array ivcof.
Ci   ivcof :Contains information about coefficients that may be varied.
Ci         :See Remarks
Ci         : ivcof(ivar,1): vector of indices to coefficients are varied
Ci         : ivcof(1:ncof,2) : contains information about constraints on
Ci         :              : each coefficient.
Ci         :  ivcof(i,2)  Value signifies:
Ci         :    0         coefficient i is free to vary
Ci         :   <0         coefficient i is free to vary, but other
Ci         :              coefficients are linked to it
Ci         :-NULLI        coefficient i is frozen; otherwise:
Ci         :   >0         change in coefficient i is linked to
Ci         :              coefficient ivcof(i,2)
Ci   dat   :data to be fit
Ci   y     :current values of fitting function
Ci   dy    :derivatives of y wrt each of 1..nvfit fitting coefficients
Ci   sig   :standard deviation sigma for each data point
Ci         :sig(i)=0 => point gets no weight (same as sigma = infinity)
Cio Inputs/Outputs
Cio  wk    :work array of dimension (ncof,5).  wk contains:
Cio        :wk(:,1) = beta vector; NR Eqns 15.5.8 and 15.5.6
Cio        :wk(:,2) = tryc: trial vector of cof
Cio        :wk(:,3) = dcof: change in cof
Cio        :wk(:,4) = pivot array in matrix inversion
Cio        :wk(:,5) = work array in matrix inversion
Cio        :wk should be preserved between successive calls to mrqstp
Cio  cof   :cof(1:ncof) contain coefficients to the fitting function
Cio        :Only a subset (nvfit<ncof) will be varied; the elements
Cio        :that are varied is set by ivcof above.
Cio        :On input, trial values for cofficients
Cio        :On output, updated trial values for coefficients
Cio        :cof must be preserved between successive calls to mrqstp
Cio  alsc  :alam scaling factor (Numerical Recipes uses alsc=10)
Cio        :If initially zero, alsc is set to 10
Cio  alp   :curvature (alpha) matrix, NR 15.5.8 and 15.5.11
Cio        :alp must be preserved between successive calls to mrqstp
Cio  alam  :L-M lambda parameter (see Remarks)
Cio        :On the first call: alam = initial value for lamdba
Cio        :If initially zero, alam starts with initial value of 0.001
Cio        :alam should be preserved between successive calls to mrqstp
Cio  iter  :number of iterations.
Cio        :On the first call, iter must be zero, which generates
Cio        :an initialization step.
Cio        :iter is incremented after each successive call.
Cio        :iter returned < 0 means an error has occured.
Cio        :iter = -1 => something wrong with ivcof array
Cio        :iter = -2 => mrqstp could not invert alpha' matrix
Co Outputs
Co   cov   :covariance matrix; see Remarks
Co   chi   :chi(1) chi-squared value for fit with current parameters
Co         :chi(2) prior chi-squared value
Cl Local variables
Cr Remarks
Cr   This routine uses the Levenberg-Marquardt algorithm for
Cr   least-squares fitting of a nonlinear function of coefficients cof
Cr   to given data.  It closely follows Numerical Recipes 2nd edition,
Cr   Chapter 15.5.
Cr   Relation between variables in NR and those here:
Cr      NR      NR fortran       here
Cr     lamdba    alamda          alam
Cr     C         covar           cov (intermediate iter: holds alpha')
Cr     a         a               cof
Cr     alpha     alpha           alp
Cr     beta      beta            wk(:,1)
Cr               atry            wk(:,2)
Cr               da              wk(:,3)
Cr     N         ndata           npfit
Cr               nca             ncof
Cr     M         mfit            nvfit
Cr               10              alsc
Cr               chisq           chi(1)
Cr               ochisq          chi(2)
Cr   Algorithm proceeds as follows.  For each iteration:
Cr     step      function
Cr      First iteration :
Cr      1        y,dy -> initial alpha, beta, chi
Cr      2        (no step 2)
Cr      3        alpha -> alpha' = alpha(1+lambda) (NR 15.5.13)
Cr               NB: alpha' stored in cov
Cr      4        alpha', beta -> trial cof (see mrqcof)
Cr               exit seeking y,dy of trial cof
Cr      Subsequent iterations have new y,dy, preserved :
Cr      1        y,dy -> trial alpha, beta, chi
Cr      2        If trial chi < old chi:
Cr                 replace cof,alpha,beta,chi with trial values
Cr                 decrease lambda
Cr               Otherwise:
Cr                 increase lambda
Cr               End of step 2
Cr      3        alpha -> alpha' = alpha(1+lambda) (NR 15.5.13)
Cr      4        alpha', beta -> trial cof (see routine mrqcof)
Cr   Calling this routine:
Cr     Iterate calls to this routine until convergence criteria is met.
Cr     (The caller must specify what that criterion is; see below)
Cr     A new set of trial cofs will be returned after each call.
Cr     Every iteration requires npfit,ncof,nvfit,ivcof,dat,y,dy,sig.
Cr     Set iter = 0 for the first iteration.
Cr     On each successive call, update y,dy; keep unchanged all other quantities
Cr
Cr     Example criteria:
Cr       1.  alam increases after e.g. 5 consecutive calls
Cr       2.  chi falls below a certain tolerance
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   30 Jul 09 First created: adapted from NR mrqmin
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nl,nbasp,nsp,ldim,nclasp,iprmb(nl**2*nbasp)
      double precision pp(6,nl,nsp,nclasp),z(nclasp)
      integer npfit,ncof,nvfit,ivcof(ncof,2),iter
      double precision alam,alsc,chi(2),alp(ncof,ncof),cov(ncof,ncof),
     .  cof(ncof),wk(ncof,5),sig(npfit),dat(npfit),y(npfit),
     .  dy(nvfit,npfit)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      real(8),allocatable :: cofl(:),wkl(:,:),alpl(:,:),ppl(:)
C ... Local parameters
      integer i,lmfitmr1,iterl,i2,ix2
      integer NULLI
      parameter (NULLI =-99999)
      integer,allocatable:: ivfit(:,:)
      double precision ppo(nl*nsp*nclasp)
      double precision xx,chil(2),alscl,alaml,alamnw

C     Hold on to everything in case user changes alam
      allocate(cofl(ncof),wkl(ncof,5),alpl(ncof,ncof),ivfit(ncof,2))
      allocate(ppl(6*nl*nsp*nclasp))
      call dcopy(ncof*5,wk,1,wkl,1)
      call dcopy(ncof,cof,1,cofl,1)
      alscl = alsc
      call dcopy(ncof*ncof,alp,1,alpl,1)
      alaml = alam
      alamnw = alam
      iterl = iter
      call dcopy(2,chi,1,chil,1)
      call dcopy(6*nl*nsp*nclasp,pp,1,ppl,1)

C     Re-entry to remake with altered lambda
   10 continue

C     Get new trial parameters, L-M algorithm
      call mrqstp(npfit,ncof,nvfit,ivcof,dat,y,dy,sig,wk,cof,alsc,alp,
     .  alam,cov,chi,iter)

C --- Copy the changes in cof (new trial parameters) to all linked-partners
C      do i1 = 1, nvfit! less than ncof
C         do i2 = 1, ncof
C           ix1 = ivcof(i1,1) !index to fitted var in ivcof
C           ix2 = ivcof(i2,2)
C           if (ix1 == ix2) then
C             print*, i2, ' is linked to ', ix2
C           endif
C         enddo
C       enddo
       do  i2  = 1, ncof
         ix2 = ivcof(i2,2)      !If >0, points to index of var actually fit
         if (ix2 > 0 .and. ix2 /= -NULLI) then
C          print*, i2, ' is linked to ', ix2
           cof(i2) = cof(i2) + (cof(ix2)-cofl(ix2))
         endif
       enddo

C     pp's in orthogonal representation
      call pptrns(1,nl,xx,nclasp,nsp,xx,nbasp,pp,ppo)

C     Update coefficients
      i = lmfitmr1(mode,s_ctrl,s_site,s_spec,
     .  nl,nbasp,nsp,nclasp,z,iprmb,
     .  ldim,ncof,nvfit,ivcof,pp,cof)

      call querym('New lambda',4,alamnw)
      if (alamnw /= alaml) then
        call dcopy(ncof*5,wkl,1,wk,1)
        call dcopy(ncof,cofl,1,cof,1)
        alsc = alscl
        call dcopy(ncof*ncof,alpl,1,alp,1)
        alaml = alamnw
        alam = alamnw
        iter = iterl
        call dcopy(2,chil,1,chi,1)
        call dcopy(6*nl*nsp*nclasp,ppl,1,pp,1)
        if (alam <= 0) goto 12
        goto 10
      endif

   12 continue
      deallocate(cofl,wkl,alpl,ppl,ivfit)

      end

      integer function lmfitmr4(job,ifmt,s,
     .  ic,clabl,l,isp,delC,delD,ifrzC,ifrzD)
C- Parses string for ppar channel info, constraints, shifts
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0  do nothing; just return
Ci         :1  read parameters ic,clabl,l,isp,delC,delD
Ci         :2  read parameters ic,clabl,l,isp,delC,delD,ifrzC,ifrzD
Ci         :   Note: if ifmt is 2, ic is not read.
Ci   ifmt  :ASA file format
Ci         :format 0 : s contains the following seven columns
Ci         :   ic clabl       l isp   shft-C    shft-D   frzC frzD
Ci         :   Example:
Ci         :     1  Zn        2  1   -0.120538   -0.000386  0  0
Ci         :format 1 : same as format 0, but all numbers may be
Ci             algebraic expressions
Ci         :format 2 : same as format1, but the first column is
Ci         :   an expression that identifies one or more classes
Ci         :   through routine slist.
Ci         :   This routine does not read ic when ifmt=2; instead
Ci         :   data is poked into any class for which the first
Ci         :   column evaluates to nonzero.  For example, all
Ci         :   classes with atomic number 25 will be read using
Ci         :       z==25
Ci         :   in the first column
Ci   s     :string containing channel data
Co Outputs
Co   ic    :class for current channel (ifmt < 2 only)
Co   clabl :class name
Co   l     :l index for current channel
Co   isp   :spin index for current channel
Co   delC  :Shift in C
Co   delD  :Shift in Delta
Co   ifrzC :flag to freeze delC (0 => vary, 1=> freeze)
Co   ifrzD :flag to freeze delD (0 => vary, 1=> freeze)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 May 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,ifmt,ic,l,isp,ifrzC,ifrzD
      double precision delC,delD
      character s*120
C ... Local parameters
      character clabl*8
      integer ix(10),a2vec,k,i,j
      character ss*120
      double precision xv(10)

      if (job == 0) then
        lmfitmr4 = 0
        return
      endif
      if (ifmt < 0 .or. ifmt > 2 .or.
     .     job < 0 .or.  job > 2) then
        lmfitmr4 = -4
        return
      endif
      lmfitmr4 = -6
      if (ifmt == 0) then
        if (job == 1) then
          read(s,*,iostat=i) ic,clabl,l,isp,delC,delD
        elseif (job == 2) then
          read(s,*,iostat=i) ic,clabl,l,isp,delC,delD,ifrzC,ifrzD
        endif
      else
        ss = s
        call word(ss,2,j,k)
        clabl = ss(j:k)
        ss(j:k) = ' '
        i = 0
        if (ifmt /= 2) then
          k = a2vec(ss,len(ss),i,4,' ',1,-2,-7,ix,xv)
        else
          call word(ss,1,j,k)
          ss(j:k) = ' '
          k = 1+a2vec(ss,len(ss),i,4,' ',1,-2,-7,ix(2),xv(2))
        endif
        if (job == 1 .and. k < 5) return
        if (job == 2 .and. k < 7) return
        if (ifmt /= 2) then
          ic = nint(xv(1))
        endif
        l    = nint(xv(2))
        isp  = nint(xv(3))
        delC = xv(4)
        delD = xv(5)
        if (job == 2) then
          ifrzC = nint(xv(6))
          ifrzD = nint(xv(7))
        endif
      endif
      lmfitmr4 = 0

      end

      integer function lmfitmr5(iv,nvar,isp,nspc,ivfit)
C- Return 1 if variable contributes to this hamiltonian, 0 otherwise
C ----------------------------------------------------------------------
Ci Inputs
Ci   iv    :Index to current parameter
Ci   nvar  :Leading dimension of ivfit
Ci   isp   :spin index of current hamiltonian
Ci   nspc  :2 for spin-coupled case, otherwise 1
Ci   ivfit :indices to coefficients that are varied
Ci         :All that is relevant here is 10's digit,
Ci         :which is spin index of parameter to be varied
Co Outputs
Co   lmfitmr5 = 1 if potential parameter has spin index contributing
Co                to hamiltonian, 0 otherwise
Cu Updates
Cu   30 Jun 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iv,isp,nspc,nvar,ivfit(nvar,2)
C ... Local parameters
      integer ispv

      ispv = mod(ivfit(iv,2)/10,10)
      lmfitmr5 = 0
      if (nspc == 2 .or. ispv == isp) lmfitmr5 = 1

      end

      subroutine lmfitmr6(mode,s_site,nl,nbasp,nsp,
     .  iprmb,ldim,ivar,ncof,ivcof,nlink,ivfit)
C- Return information about coefficients associated with variable ivar
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :10s digit controls what to vary
Ci         :1 ASA C and Delta
Ci         :2 ASA C+enu and Delta
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec class clabel
Ci     Stored: *
Ci     Passed to: *
Ci   nl    :(global maximum l) + 1
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ivar  :independent variable for which information is returned
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   ncof  :number of coefficients that can potentially vary
Ci   ivcof :Contains information about coefficients that may be varied.
Ci         :See Remarks
Ci         : ivcof(ivar,1): vector of indices to coefficients are varied
Ci         : ivcof(1:ncof,2): contains information about constraints on
Ci         :              each coefficient.
Ci         :  ivcof(i,2)  Value signifies:
Ci         :    0         coefficient i is free to vary
Ci         :   <0         coefficient i is free to vary, but other
Ci         :              coefficients are linked to it
Ci         :-NULLI        coefficient i is frozen; otherwise:
Ci         :   >0         change in coefficient i is linked to
Ci         :              coefficient ivcof(i,2)
Co Outputs
Co   nlink :1 + number of other coefficients linked to the
Co         :ivar'th coefficient to be varied.
Co         :nlink=1 => only the coefficient is to be varied (no linking)
Co   ivfit :vector of information about the ivar'th coefficient that is
Co         :varied, including all coefficients linked to the ivar'th one.
Co         :nlink elements are returned
Co         :ivfit(1:nlink,1) flags type of coefficient to vary
Co         :ivfit(1:nlink,2) flags which coefficient within a type
Co         :So far, only ASA variations C and Delta are implemented (mode=10)
Co         :In that case: pp depends on l,isp, and class
Co         :ivfit specifies which pp as follows:
Co         : ivfit(:,1)   ivfit(:,2)        type
Co         :    2         l+10*isp+100*ic   C parameter
Co         :    3         l+10*isp+100*ic   Delta
Cl Local variables
Cr Remarks
Cr         :Coefficients are one of three types:
Cr         : - coefficient is free to vary
Cr         : - coefficient is frozen
Cr         : - change coefficient is linked to change in another coefficient.
Cr         :Let nvar be the number of such coefficients.
Cr         :Knowledge of nvar is not required by this routine, but
Cr         :ivar must not exceed nvar.
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Aug 10 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nl,nbasp,nsp,ldim,ncof,ivar,nlink
      integer iprmb(nl**2*nbasp),ivcof(ncof,2),ivfit(ncof,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: ibc(:)
C ... Local parameters
      integer lmr,ib,l,ic,ncofl,isp,ii(10),nvarl,mode1,jcof
      double precision xx
C     integer NULLI
C     parameter (NULLI =-99999)

      if (ncof == 0) return

C --- Loop over sites (skip classes already touched inside loop) ---
      ncofl = 0
      nvarl = 0
      lmr = 0
      nlink = 0
      mode1 = mod(mode/10,10)
C     Doubles as a flag ... if nonzero, one cofficient has been found
      ivfit(1,1) = 0
      allocate(ibc(nbasp)); call iinit(ibc,nbasp)
      do  ib = 1, nbasp

        ic = s_site(ib)%class

C   --- mode 1 : vary C and Delta ---
        if (mode1 == 1 .or. mode1 == 2) then

C   ... Skip classes already touched
        if (ibc(ic) /= 0) then
          lmr = lmr + nl**2
        else

        do  l = 0, nl-1
          do  isp = 1, nsp

C           Case Rl channel is in lower block
            if (iprmb(lmr+1) < ldim) then

              ii(1) = l
              ii(2) = isp
              ii(3) = ic

C             C parameter
              if (ivcof(ncofl+1,2) <= 0 .and. ivfit(1,1) == 0) then
                nvarl = nvarl + 1
                if (nvarl == ivar) ivfit(1,1) = 2
                jcof = ivcof(nvarl,1)
              endif
C             D parameter
              if (ivcof(ncofl+2,2) <= 0 .and. ivfit(1,1) == 0) then
                nvarl = nvarl + 1
                if (nvarl == ivar) ivfit(1,1) = 3
                jcof = ivcof(nvarl,1)
              endif

C         ... Increment number of links; load ivfit
              if (ivfit(1,1) /= 0) then
                if (ivcof(ncofl+1,2) == jcof .or.
     .              ivcof(ncofl+2,2) == jcof .or. nlink == 0) then
                  nlink=nlink+1
                  ivfit(nlink,1) = ivfit(1,1)
                  call lmfitmr2(0,nlink,ii,ncof,nl,nsp,ivfit,xx)
                endif
              endif

            endif
            ncofl = ncofl + 2

          enddo ! Loop over spin
          lmr = lmr + 2*l+1
        enddo ! Loop over l
        endif ! done with this class

C   --- Error if mode not recognized ---
        else
          call rxi(' LMFITMR: mode not recognized, mode ',mode1)
        endif
        ibc(ic) = ib ! Flags class has been treated

      enddo
      deallocate(ibc)

C     print *, ivfit(1:nlink,2)

      end
