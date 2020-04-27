      subroutine lmfopb(sopts,s_lat,s_spec,ehf,ninit)
C- Optimize lmf basis
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  gmax
Co     Stored:     gmax
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name lmxb z rmt lmxa pz orbp idxdn p norp ntorb
Co     Stored:     orbp norp ntorb p pz idxdn lmxb
Co     Allocated:  *
Cio    Elts passed:pz name
Cio    Passed to:  uspecb lmfop2 ioorbp
Ci Inputs
Ci   sopts :Switches controlling optimization process.
Ci         :1st char is delimiter separating switches.  Switches:
Ci         :wbas  write current basis parameters to file basp2 and exit
Ci         :atm   (not implemented)
Ci         :init  density initialization (not checked for a long time)
Ci         :etol  only adjust parameter if delta E > etol
Ci         :sort  sort parameters according to increasing rs
Ci         :rs    specify RSM in optimization
Ci         :e     specify EH in optimization
Ci         :es    Ditto for empty spheres
Ci         :spec  Optimize for particular species only
Ci   ehf   :energy this iteration
Co Outputs
Co   ninit :number of initial iterations to converge density before
Co         :optimizing
Cl Local variables
Cl   ir    :0 optimization completed
Cl   ivar  :index to current variable for minimization
Cl   ncall :number of function calls for which energy aleady calculated
Cl         :for variable being minimized.
Cl   iter  :number of function calls so far in search for min.
Cl         :When iter reaches ncall, lmfopb returns with ir=-1
Cl   vhist :table of energy values for prior 1..ncall points
Cr Remarks
Cr
Cu Updates
Cu   31 Dec 16 wbas now writes only after lpass is nonzero (so P is taken from rst)
Cu   07 Feb 14 global --optbas can be --optbas:rs or --optbas:e or --optbas:rs:e
Cu   10 Nov 11 Begin migration to f90 structures
Cu   14 Nov 09 Some bug fixes
Cu   09 Jul 09 Some bug fixes, new etol
Cu   10 May 09 Some improvements for version 7
Cu   28 Aug 04 Interactive mode queries whether to restore
Cu   07 Sep 03 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
      double precision ehf
      integer ninit
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical lrsm,leh,ltmp,ioorbp
      integer mxcall,nvmax
      parameter (mxcall=20,nvmax=100)
      integer a2vec,i,is,j1,j2,j3,nspec,nglob,llst(6),lmxb,nvar,iter,it,ncall,
     .  ivar,ifi,fopna,ir,lsort,linces,lpass,ipr,l,lt,stdo,stdl,nkaph,lpz
      integer procid,master,mpipid
      integer optlst(3,nvmax),optls0(3,nvmax),iwk(nvmax,5),isw
      character dc*1,slabl*8,strnl*6,pnam(2)*3
      double precision val,rmt,ax,bx,cx,xmin,xmax,fa,fb,fc,fopb,z
      double precision vhist(2,mxcall),wk(nvmax),estart
      common /lmfop/ vhist,iter,ncall
      external fopb
C ... for lmf basis
      integer nkap0,n0
      parameter (nkap0=4,n0=10)
      integer lh(nkap0),nkapi,ik,ninit0
      double precision rsmh(n0,nkap0),eh(n0,nkap0),e,rsm,tol,gmax,etol
      save ivar,xmin,xmax,tol,lpass,optlst,estart,gmax,ninit0,etol
      data ivar /0/ lpass /0/ estart /0d0/ ninit0 /0/ etol /5d-4/


C     call shstru('spec',sspec,1,1)

      pnam(1) = 'rsm'
      pnam(2) = 'eh'
      nspec = nglob('nspec')
      nkaph = nglob('nkaph')
      lpz = 0 ; if (nglob('lpz') > 0) lpz = 1
      call getpr(ipr)
      stdo = nglob('stdo')
      stdl = nglob('stdl')
      if (ninit0 == 0) ninit = 0  ! Assume no initialization pases
      procid = mpipid(1)
      master = 0
      lrsm = .false.
      leh  = .false.
      lsort = 0
      linces = 0

C --- Parse arguments in sopts ---
      j1 = 1
      dc = sopts(j1:j1)
      j1 = j1+1
      nvar = 0

C ... Return here to resume parsing for arguments
   40 continue
      call nwordg(sopts,0,dc//' ',1,j1,j2)
      if (j2 >= j1) then
C       print *, sopts(j1:j2)
        if (sopts(j1:j2) == 'wbas')  then
          if (lpass /= 0) then
            if (procid == master) then
              ifi = fopna('basp2',-1,0)
              rewind ifi
              ltmp = ioorbp(110001,nkaph-lpz,1,nspec,s_spec,0,-ifi)
            endif
            call rx0('saved basis in file basp2')
          endif
        elseif (sopts(j1:j2) == 'atm')  then
          call rx('lmfopb has not yet implemented atm switch')
C       Handle initial cd passes
        elseif (sopts(j1:j1+4) == 'init=')  then
          if (ninit0 == 0) then
            i = j1+4
            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,llst,ninit) /= 1)
     .        goto 999
            ninit = ninit + 1
            ninit0 = ninit
            gmax = s_lat%gmax
          endif
          s_lat%gmax = gmax ! Preserve variable each iter
          ninit = ninit - 1
          if (ninit > 0) then
            call info2(10,1,0,' LMFOPB: start density initialization,'//
     .        ' iteration %i',ninit0-ninit,0)
            return
          elseif (ninit == 0) then
            call info2(10,1,0,' LMFOPB: completed %i iterations '//
     .        'density initialization',ninit0-1,0)
          endif
C       If optimization step improves E by < etol, keep par unchanged
        elseif (sopts(j1:j1+4) == 'etol=')  then
          i = j1+4
          if (a2vec(sopts,j2,i,4,' '//dc,2,1,1,llst,etol) /= 1)
     .      goto 999
        elseif (sopts(j1:j2) == 'sort')  then
          lsort = 1
        elseif (sopts(j1:j2) == 'rs')  then
          lrsm = .true.
        elseif (sopts(j1:j2) == 'e')  then
          leh = .true.
        elseif (sopts(j1:j2) == 'es')  then
          linces = 1
        elseif (sopts(j1:j1+4) == 'spec=')  then
          call nwordg(sopts,0,dc//', ',1,j1,j3)
          call locase(sopts(j1+5:j3))
          do  is = 1, nspec
            slabl = s_spec(is)%name
            call locase(slabl)
C           Found a species: add orbitals to list
            if (slabl == sopts(j1+5:j3)) then
              strnl = '012345'
              lrsm = .false.
              leh  = .false.

C         ... Return here to resume parsing for species arguments
  140         continue
              if (j3 >= j2) goto 150
              j1 = j3+2
              call nwordg(sopts,0,','//dc//' ',1,j1,j3)
              if (sopts(j1:j3) == 'rs')  then
                lrsm = .true.
                goto 140
              elseif (sopts(j1:j3) == 'e')  then
                leh = .true.
                goto 140
              elseif (sopts(j1:j1+1) == 'l=')  then
                strnl = sopts(j1+2:j3)
                goto 140
              else
                goto 999
              endif

C         ... Add specification to list
  150         continue
              if (.not. lrsm .and. .not. leh) lrsm = .true.
              call strip(strnl,j1,j3)
              read(strnl,'(6I1)') (llst(j1), j1=1,j3)
              lmxb = s_spec(is)%lmxb
C             Exclude states l>lmxb
              do  j1 = 1, j3
                if (llst(j1) > lmxb) llst(j1) = -1
              enddo
              call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
C             call lmfop2(0,is,sspec,s_spec,rsmh,eh)
C             Exclude states rsm=0
              ik = 1
              do  j1 = 1, j3
                if (nkapi == 0) llst(j1) = -1
                l = llst(j1)
                if (l >= 0) then
                  e = eh(l+1,ik)
                  rsm = rsmh(l+1,ik)
                  if (rsm <= 0) llst(j1) = -1
                endif
              enddo

C         ... Add to list, smoothing radii
              if (lrsm) then
                do  j1 = 1, j3
                  if (llst(j1) >= 0) then
                    nvar = nvar+1
                    optls0(1,nvar) = is
                    optls0(2,nvar) = llst(j1)
                    optls0(3,nvar) = 1
                  endif
                enddo
              endif

C         ... Add to list, smoothing energies
              if (leh) then
                do  j1 = 1, j3
                  if (llst(j1) >= 0) then
                    nvar = nvar+1
                    optls0(1,nvar) = is
                    optls0(2,nvar) = llst(j1)
                    optls0(3,nvar) = 2
                  endif
                enddo
              endif

              j1 = j2+2
              goto 40
            endif
          enddo
          call rxs('lmfop: specified nonexistent species: ',
     .      sopts(j1:j3))
        else
          goto 999
        endif
        j1 = j2+2
        goto 40
      endif

C ... If no variables explicitly specified, optimize all rs
      if (nvar == 0) then
      if (.not. lrsm .and. .not. leh) then
        lrsm = .true.  ! Default is to optimize wrt rs only
      endif
      if (lrsm) then
        do  is = 1, nspec
          lmxb = s_spec(is)%lmxb
          z = s_spec(is)%z
          call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
          if (nkapi == 0) cycle
          if (z == 0 .and. linces == 0) cycle
          ik = 1
          do  l = 0, lmxb
            if (rsmh(l+1,ik) <= 0) cycle
            nvar = nvar+1
            optls0(1,nvar) = is
            optls0(2,nvar) = l
            optls0(3,nvar) = 1
          enddo
        enddo
      endif
      if (leh) then
        do  is = 1, nspec
          lmxb = s_spec(is)%lmxb
          z = s_spec(is)%z
          call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
          if (nkapi == 0) cycle
          if (z == 0 .and. linces == 0) cycle
          ik = 1
          do  l = 0, lmxb
            if (rsmh(l+1,ik) <= 0) cycle
            nvar = nvar+1
            optls0(1,nvar) = is
            optls0(2,nvar) = l
            optls0(3,nvar) = 2
          enddo
        enddo
      endif
      endif

      if (nvar == 0)
     .  call rx0('lmfopb: no parameters to optimize ... quitting')

C --- Order variables to be minimized ---
      if (lsort == 1 .and. lpass == 0) then

C       Order rs variables before e variables
        call icopy(nvar,optls0(3,1),3,iwk,1)
        call ivheap(1,nvar,iwk,iwk(1,2),101)
        call ivprm(3,nvar,optls0,iwk(1,3),iwk(1,2),1)

C       Count the number of rs variables
        j1 = 0
        do  i = 1, nvar
          if (optls0(3,i) == 1) j1 = i
        enddo

C       Put rs variables in ascending order
        do  i = 1, j1
          is = optls0(1,i)
          l  = optls0(2,i)
          lt = optls0(3,i)
          rmt = s_spec(is)%rmt
          slabl = s_spec(is)%name
          call lmfop2(0,is,s_spec,rsmh,eh)
          if (rmt > 0) then
            wk(i) = min(rsmh(l+1,ik),rmt)
          else
            wk(i) = rsmh(l+1,ik)
          endif
        enddo
        call dvheap(1,j1,wk,iwk(1,2),1d-3,101)
        call ivprm(3,j1,optls0,iwk(1,3),iwk(1,2),1)

C       Put e variables in ascending order
        do  i = j1+1, nvar
          is = optls0(1,i)
          l  = optls0(2,i)
          lt = optls0(3,i)
          rmt = s_spec(is)%rmt
          slabl = s_spec(is)%name
          call lmfop2(0,is,s_spec,rsmh,eh)
          wk(i-j1) = eh(l+1,ik)
        enddo
        call dvheap(1,nvar-j1,wk,iwk(1,2),1d-3,101)
        call ivprm(3,nvar-j1,optls0(1,j1+1),iwk(1,3),iwk(1,2),1)
      endif

C --- First pass: copy sorted table to permanent array ---
      if (lpass == 0) then

C       Keep static record of these variables; reset each new pass
        gmax = s_lat%gmax

        call icopy(3*nvar,optls0,1,optlst,1)
        lpass = 1
      else
C       Restore these variables
        s_lat%gmax = gmax
      endif

C --- Minimization procedure for next variable ---
   60 continue
      if (nvar == 0) then
        call info0(10,1,0,
     .    ' LMFOPB (warning): no orbitals specified to optimize')
      endif

C ... Exit when optimization complete
      if (ncall == 0 .and. ivar == nvar) then
        call info2(20,1,0,' LMFOPB:  before optimization ehf=%,4;4d.'//
     .    '  After optimization ehf=%,4;4d',estart,vhist(2,iter))
        call rx0('Optimization complete.  See file basp2')
      endif

C ... Setup for new variable, or re-entrance after generation of ehf
      if (ncall == 0) then
        vhist(1,1) = 0
        vhist(1,2) = 0
        ivar = ivar+1
      else
        vhist(2,ncall) = ehf
        if (estart == 0) estart = ehf
      endif

C --- Initial Setup ---
      if (ncall == 0 .and. ivar == 1 .and. ipr >= 20) then

        call info5(10,1,0,
     .    ' LMFOPB:  optimizing energy wrt %i parameters'//
     .    '%?#n#, etol=%;3g#%j#:'//
     .    '%N   spec%7fl  type     start',
     .    nvar,isw(etol > 0),etol,0,0)
        ik = 1
        do  j1 = 1, nvar
          is = optlst(1,j1)
          l  = optlst(2,j1)
          lt = optlst(3,j1)

          rmt = s_spec(is)%rmt
          slabl = s_spec(is)%name
C         call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
          call lmfop2(0,is,s_spec,rsmh,eh)
C         Starting guesses; also set limit to range
          if (lt == 1) then
            if (rmt > 0) then
              val = min(rsmh(l+1,ik),rmt)
            else
              val = rsmh(l+1,ik)
            endif
            xmin = min(1d0,val-.1d0)
            xmax = rmt
          endif
          if (lt == 2) then
            val = eh(l+1,ik)
            xmin = min(-1d0,val-.5d0)
            xmax = -.05d0
          endif
          write(stdo,333) is,slabl,l,pnam(lt),val
  333     format(i4,':',a8,i2,3x,a3,2x,f8.3)
        enddo
      endif

C ... Get parameters for this variable
      is = optlst(1,ivar)
      l  = optlst(2,ivar)
      lt = optlst(3,ivar)
C     call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
      call lmfop2(0,is,s_spec,rsmh,eh)
      rmt = s_spec(is)%rmt

C ... Starting guesses; also set limit to range
      if (ncall == 0) then
        if (lt == 1) then
          if (rmt > 0) then
            val = min(rsmh(l+1,ik),rmt)
          else
            val = rsmh(l+1,ik)
          endif
          xmin = min(1d0,val-.1d0)
          xmax = rmt
          if (rmt == 0) xmax = 2.5d0
          tol = .05d0
        endif
        if (lt == 2) then
          val = eh(l+1,ik)
          xmin = min(-1d0,val-.5d0)
          xmax = -.05d0
          tol = .04d0
        endif
        ax = val
        if (val < (xmin+xmax)/2) then
          cx = min(val+.2d0,xmax)
          if (abs(cx-ax) < .1d0) ax = ax-.1d0
        else
          cx = max(val-.2d0,xmin)
          if (abs(cx-ax) < .1d0) ax = ax+.1d0
        endif
        vhist(1,1) = ax
        vhist(1,2) = cx

        slabl = s_spec(is)%name
        call info5(20,1,0,' LMFOPB: begin optimization of var #%i in'//
     .    ' species '//slabl//'%a: start=%;3d range=(%,3;3d,%,3;3d)',
     .    ivar,ax,xmin,xmax,0)
      else
        slabl = s_spec(is)%name
        call info2(20,1,0,' LMFOPB: continue optimization of var #%i,'//
     .    ' species '//slabl//'%a val=%d',ivar,vhist(1,ncall))
      endif

C --- Either bracket minimum or bump against constraint ---
      iter = 0
      ax = vhist(1,1)
      cx = vhist(1,2)
      call mnbrak(ax,bx,cx,fa,fb,fc,fopb,xmin,xmax,1,mxcall,it,ir)
      if (ir == -2) vhist(1,iter) = cx
      if (lt == 1) then
        rsmh(l+1,ik) = vhist(1,iter)
      endif
      if (lt == 2) then
        eh(l+1,ik) = vhist(1,iter)
      endif
      call lmfop2(1,is,s_spec,rsmh,eh)
      if (ir == -2) then
        call info8(20,1,0,' LMFOPB: var #%i: min outside range'//
     .    ' (%,3;3d,%,3;3d) ... use %,3;3d, ehf=%;6d,'//
     .    '  delta ehf=%;3g',ivar,xmin,xmax,
     .    vhist(1,iter),vhist(2,iter),vhist(2,iter)-vhist(2,1),0,0)
        if (ipr > 1) call awrit4('op var #%i limit to %;3d:  ehf=%;6d, dehf=%;3g',
     .    ' ',120,stdl,ivar,vhist(1,iter),vhist(2,iter),vhist(2,iter)-vhist(2,1))
        ncall = 0
        if (procid == master) then
          ifi = fopna('basp2',-1,0)
          rewind ifi
          ltmp = ioorbp(110001,nkaph-lpz,1,nspec,s_spec,0,-ifi)
          call fclr(' ',ifi)
        endif
        j1 = 1
        if (etol > 0 .and. -(vhist(2,iter)-vhist(2,1)) < etol) j1 = 0
        call query('enter <ret> or 0<ret> to '//
     .    'continue or restore starting value before continuing',2,j1)
        if (j1 == 0) then
          if (lt == 1) then
            rsmh(l+1,ik) = vhist(1,1)
          endif
          if (lt == 2) then
            eh(l+1,ik) = vhist(1,1)
          endif
          call info2(20,0,0,' LMFOPB: var #%i restored to %;3d ',ivar,
     .      vhist(1,1))
          call lmfop2(1,is,s_spec,rsmh,eh)
          if (procid == master) then
            ifi = fopna('basp2',-1,0)
            rewind ifi
            ltmp = ioorbp(110001,nkaph-lpz,1,nspec,s_spec,0,-ifi)
            call fclr(' ',ifi)
          endif
        endif
        goto 60
      endif
C ... mnbrak wants another function call
C     if (iter >= ncall) then
      if (ir == -1) then
        ncall = ncall+1
        return
      endif

C --- Minimize bracketed function ---
      call brent(ax,bx,cx,fa,fb,fc,fopb,1,tol,.01d0,0,ir)

      if (ir >= 0) vhist(1,iter) = bx
      if (ir >= 0) vhist(2,iter) = fb

      if (lt == 1) then
        rsmh(l+1,ik) = vhist(1,iter)
      endif
      if (lt == 2) then
        eh(l+1,ik) = vhist(1,iter)
      endif
      call lmfop2(1,is,s_spec,rsmh,eh)

      if (ir >= 0) then
        call info5(20,1,0,' LMFOPB: var #%i optimized to %;3d (%i'//
     .    ' iter):  ehf=%;6d, delta ehf=%;3g',
     .    ivar,bx,iter,fb,fb-vhist(2,1))
        if (ipr > 1) call awrit4('op var #%i opt to %;3d:  ehf=%;6d, dehf=%;3g',
     .    ' ',120,stdl,ivar,bx,fb,fb-vhist(2,1))
        if (procid == master) then
          ifi = fopna('basp2',-1,0)
          rewind ifi
          ltmp = ioorbp(110001,nkaph-lpz,1,nspec,s_spec,0,-ifi)
          call fclr(' ',ifi)
        endif
        j1 = 1
        if (etol > 0 .and. -(fb-vhist(2,1)) < etol) j1 = 0
        call query('enter <ret> or 0<ret> to '//
     .    'continue or restore starting value before continuing',2,j1)
        if (j1 == 0) then
          if (lt == 1) then
            rsmh(l+1,ik) = vhist(1,1)
          endif
          if (lt == 2) then
            eh(l+1,ik) = vhist(1,1)
          endif
          call info2(20,0,0,' LMFOPB: var #%i restored to %;3d ',ivar,
     .      vhist(1,1))
          call lmfop2(1,is,s_spec,rsmh,eh)
          if (procid == master) then
            ifi = fopna('basp2',-1,0)
            rewind ifi
            ltmp = ioorbp(110001,nkaph-lpz,1,nspec,s_spec,0,-ifi)
            call fclr(' ',ifi)
          endif
        endif
C      else
C        call info2(20,1,0,' LMFOPB:  found=F   x=%d',vhist(1,iter),0)
      endif


C ... brent wants another function call
      if (ir == -1) then
        ncall = ncall+1
        return
      endif

C     Done with this variable; start another
      ncall = 0
      goto 60

C     return
  999 call rxs('lmfopb: failed to parse switch: ',sopts(j1:j2))
      end
      subroutine lmfop2(lpack,is,s_spec,rsmh,eh)
C- Pack or unpack species parameters
      use structures
      implicit none
C ... Passed parameters
      integer lpack,is
      integer nkap0,n0
      parameter (nkap0=4,n0=10)
C     integer lh(nkap0),nkapi,npqni,ik
      double precision rsmh(n0,nkap0),eh(n0,nkap0)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ikap
      double precision orbp(n0,2,nkap0)

      call dpzero(orbp,n0*nkap0*2)
      orbp = s_spec(is)%orbp
C ... Unpack
      if (lpack == 0) then
C       call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
        ikap = 1
        call dcopy(n0,orbp(1,1,ikap),1,rsmh(1,ikap),1)
        call dcopy(n0,orbp(1,2,ikap),1,eh(1,ikap),1)
      else
        ikap = 1
        call dcopy(n0,rsmh(1,ikap),1,orbp(1,1,ikap),1)
        call dcopy(n0,eh(1,ikap),1,orbp(1,2,ikap),1)
        s_spec(is)%orbp = orbp
      endif
      end

      double precision function fopb(x,nfcn0,ir)
      double precision x(*)
      integer nfcn0,ir

      integer mxcall,ncall,iter
      parameter (mxcall=20)
      double precision vhist(2,mxcall)
      common /lmfop/ vhist,iter,ncall

      iter = iter+1
      if (iter <= ncall) then
        if (x(1) /= vhist(1,iter))
     .    call rxi('lmfopb: history corrupted iter %i .. aborting',iter)
        ir = iter
        fopb = vhist(2,iter)
      else
        ir = -1
        if (iter > mxcall)
     .    call rxi('lmfopb: failed to converge in %i calls',mxcall)
        vhist(1,iter) = x(1)
        fopb = 0
      endif
      end
