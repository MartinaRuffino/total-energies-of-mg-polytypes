C#define PRTNOCR
      subroutine popted(sopts,s_ctrl,s_lat,s_bz)
C- popt file editor
C ----------------------------------------------------------------------
Cio Structures
Ci   s_ctrl (not used)
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat qlat
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: nkabc nkp
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:wtkp qp ipq
Cio    Passed to: *
Ci Inputs/Outputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Ci   popted never returns.
Ci   rst file can be written.
Cr Remarks
Cl Local variables
Cl   havep :0 no polarizations read
Cl         :1 polarizations available for irreducible k
Cl         :2 polarizations available for some list of k
Cu Updates
Cu   08 Nov 09 Some extensions for dielectric response
Cu   19 Dec 08 New exch option
Cu   01 Jul 08 First created
C  ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer ifi,n0
      parameter (n0=10)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_bz)::    s_bz
C ... Local parameters
      integer fopna,fopng,stdo,nglob
      integer i,j,k,j1,j2,js1,js2,nr,nc,rdm,nw,a2vec,ios,iwk(10),iwk2(4)
      logical lnsave,lbin,lsopts,irr
      integer havep,nq
      integer n1,n2,n3,nkabc(3),nkp,npol,ipol,ifac(3),ix(10)
      integer i11,i12,i21,i22,i31,i32,iq,jq,iq1,i3,i2,i1,ival,
     .  jj1,jj2,jj3,iindx(3,2)
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      equivalence (i11,iindx(1,1)),(i21,iindx(2,1)),(i31,iindx(3,1))
      equivalence (i12,iindx(1,2)),(i22,iindx(2,2)),(i32,iindx(3,2))
      integer,allocatable :: iprm(:)
      real(8),allocatable :: pfile(:,:),poptk(:,:,:),poptk2(:,:,:),
     .  wtkp(:),y2(:)

      real(8),allocatable:: qp(:,:),emesh(:),wk(:)
C     real(8),allocatable:: qfbz(:,:)

      double precision xx,xx2,omega,dval,sumd,qpi(3)
C      double precision qk
      double precision qb(3,3),alat,plat(3,3),qlat(3,3),q(3)
      character dc*1, fn*120, outs*150, strn*150
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
C      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
C     .                    (jj2*ifac(2)-1)*qb(k,2) +
C     .                    (jj3*ifac(3)-1)*qb(k,3)
C ... data statements
C     data vec0 /0d0,0d0,0d0/

      stdo = nglob('stdo')
      havep = 0
      lnsave = .false.
      npol = 1
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      allocate(wtkp(nkp))
      do  i = 1, nkp
        wtkp(i) = dval(s_bz%wtkp,i)*0.5d0*n1*n2*n3
      enddo
C     call dinv33(plat,1,qlat,vol)
      qb(1:3,1) = qlat(1:3,1)/n1
      qb(1:3,2) = qlat(1:3,2)/n2
      qb(1:3,3) = qlat(1:3,3)/n3
      ifac = 1
      irr = .true.

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the popt file editor. ',
     .    'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the popt file editor.  ',
     .    'Enter ''?'' to see options.')
        lsopts = .false.
      endif

C ... Return here to resume parsing for arguments
   10 continue
      if (lsopts) then
        js2 = js2+1
        if (js2 > len(sopts)) then
          lsopts = .false.
          goto 10
        endif
        if (sopts(js2:js2) == dc) goto 10
        js1 = min(len(sopts),js2)
        call nwordg(sopts,0,dc,1,js1,js2)
        if (js2 < js1) lsopts = .false.
      endif

C 306 format(' Failed to parse string ',a,' ... try again.')
C 100 continue
C#ifdef PRTNOCR
      print '(/'' Option : '',$)'
C#elseC
C      print '(/'' Option : '')'
C#endif
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      call locase(outs)

C ... Parse and execute the next command
      if (.false.) then

      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',
     .    ' ''?'' to see menu')
        goto 10

      elseif (outs(1:8) == 'readtk ' .or. outs(1:9) == 'readtkl ') then

        call word(outs,2,j1,j2)
        if (j2 < j1) then
          irr = .true.
          nq = nkp
          call info2(0,0,0,'%10preading %i irr k-points, Kotani style',
     .      nq,0)
        else
          j = 0
          j = a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-3,1,ix,npol)
          call info2(2,1,0,' reading %i qp',j,0)
          if (j <= 0) goto 98
          irr = .false.
          nq = j
          call info2(0,0,0,'%10preading %i k-points, Kotani style',nq,0)
        endif
        if (allocated(qp)) deallocate(qp)
        allocate(qp(3,nq))


C   ... For each k-point, do
        do  iq = 1, nq
C         Make Kotani-style file name
          fn = 'EPS0000.nlfc.dat'
          if (outs(1:8) == 'readtk ') fn = 'EPS0000.dat'
          write(fn(4:7),'(i4)') iq
          do  j = 4, 7
            if (fn(j:j) == ' ') fn(j:j) = '0'
          enddo

C         Open file, read past header
          ifi = fopng(trim(fn),-1,1)
          rewind ifi
          read(ifi,*) strn

C         Get energy mesh ; put into wk
          if (allocated(wk)) deallocate(wk)
          allocate(wk(100000))
          j = 0
          do while (.true.)
            read(ifi,*,end=30,err=30) qpi,wk(j+1)
            j = j + 1
          enddo
   30     continue

C         First qp: allocate and assign
          if (iq == 1) then
            nr = j
            nc = nq+1
            if (allocated(emesh)) deallocate(emesh)
            allocate(emesh(nr))
            call dcopy(nr,wk,1,emesh,1)
            allocate(pfile(nr,nc))
            call info2(0,0,0,'%10pFile '//trim(fn)//
     .        ' has %i energy points',nr,0)
            if (allocated(iprm)) deallocate(iprm)
            allocate(iprm(nq)); call iinit(iprm,nq)
          endif
C         Subsequent qp: assignments and sanity checks
          if (nr /= j) then
            call info5(0,0,0,'%10pabort, qp %i:  expected %i '//
     .        'energy points from but read %i',iq,nr,j,0,0)
            goto 98
          endif
          call info2(2,0,0,'%10preading iq=%,4i q =%3;8,5D',iq,qpi)
          qp(:,iq) = qpi
          if (irr) then
            call findqbz(nkp,s_bz%qp,plat,1d-5,qpi,jq)
C            print *, iq,jq
            if (iprm(jq) == 1) then
              call info2(2,0,0,
     .          '%4poops! this qp already read: qp = %3;8,5D',qpi,0)
              goto 98
            endif
            if (jq == -1) then
              call info2(2,0,0,
     .          '%4poops! cannot match qp %i (= %3;8,5D) to any qp',
     .          iq,qpi)
              goto 98
            endif
            iprm(jq) = 1
          endif

C         Read eps
          rewind ifi
          read(ifi,*) strn
          do  j = 1, nr
            read(ifi,*,end=30,err=30) qpi,pfile(j,1),pfile(j,jq+1)
            if (emesh(j) /= pfile(j,1)) then
              call info2(2,0,0,
     .          '%4poops! energy mismatch, ie=%i iq=%i',j,iq)
              goto 98
            endif
          enddo


          wtkp(jq) = dble(n1)*dble(n2)*dble(n3)

C         Cleanup for this qp
          call fclr(' ',ifi)

        enddo
        havep = 1
        if (irr) havep = 2
        call info2(0,0,0,'%10p%i qp and %i energy points read',nq,nr)

        print *, 'wtkp set to constant.'
        print *, 'Future: replace wtkp in kmape with some general wt'


      elseif (outs(1:5) == 'read ' .or. outs(1:6) == 'readq '
     .                             .or. outs(1:6) == 'readb ') then

        lbin = outs(1:6) == 'readb '
        call word(outs,2,j1,j2)
        if (j2 < j1) then
          fn = 'popt'
          if (lbin) then
            fn = 'poptb'
            ifi = fopna(fn,-1,4)
          else
            ifi = fopna(fn,-1,1)
          endif
        else
          fn = outs(j1:j2)
          if (lbin) then
            ifi = fopna(fn,-1,4)
          else
            ifi = fopng(outs(j1:j2),-1,1)
          endif
        endif
        if (allocated(pfile)) deallocate(pfile)
        nr = 0; nc = 0
        rewind ifi
        j = 0; if (lbin) j = 2
        i = rdm(ifi,j,0,' ',xx,nr,nc)
        if (i == 1) then
          allocate(pfile(nr,nc))
          rewind ifi
          i = 0
          if (outs(1:6) == 'qread ') i = 1
          if (lbin) i = 2
          i = rdm(ifi,i,nr*nc,' ',pfile,nr,nc)
        endif
        if (i /= 1) then
          print *, 'failed to read file '//trim(fn)//'... nothing done'
          deallocate(pfile)
          goto 10
        endif
        havep = 2
        call info2(2,0,0,' read %i energies, %i spectra from file '//
     .    trim(fn),nr,nc-1)

        goto 10

C ... npol
      elseif (outs(1:5) == 'npol ') then

        call words(outs,nw)
        if (nw < 2) goto 98
        call word(outs,2,j1,j2)
        j = 0
        j = a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-3,1,ix,npol)

        call info2(2,1,0,' popt file has %i polarizations',npol,0)
        goto 10

C ... scale optical lines
      elseif (outs(1:7) == 'scalep ') then

        if (havep == 0) goto 97
        call words(outs,nw)
        if (nw < 2) goto 98
        call word(outs,2,j1,j2)
        j = 0
        j = a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-3,1,ix,xx)

        call info2(2,1,0,' scale popt by %d, %i spectra',xx,nc-1)
        do  i = 1, nr
          do  j = 2, nc
            pfile(i,j) = xx * pfile(i,j)
          enddo
        enddo
        goto 10

C ... kmape, kshowe
      elseif (outs(1:6) == 'kmape ' .or. outs(1:7) == 'kshowe ') then

        if (havep == 1) print *, 'popt not available at irr k'
        if (havep /= 2) goto 97
        call words(outs,nw)
        if (nw < 2) goto 98
        call word(outs,2,j1,j2)
        j = 0
        if (a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-3,1,ix,omega)
     . /= 1) goto 99

C       Read polarization
        ipol = 1
        if (nw >= 3) then
          call word(outs,3,j1,j2)
          j = 0
          if (a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-3,1,ix,ipol)
     . /= 1) goto 99
        endif

        if (nc-1 /= nkp*npol) then
          call info5(0,0,0,'    oops! number of file spectra (%i) '//
     .      ' doesn''t match nkp*npol (%i*%i)',nc-1,nkp,npol,0,0)
          goto 98
        endif

C   ... kmap branch
        if (outs(1:6) == 'kmape ') then

C       Optionally select a single index in one of 3 k vectors
        do  i = 1, 3
          iindx(i,1) = 1; iindx(i,2) = nkabc(i)
        enddo
        if (nw >= 4) then
          call word(outs,4,j1,j2)
          if (outs(j1:j1) /= 'i') goto 99
          if (outs(j1+2:j1+2) /= '=') goto 99
          read(outs(j1+1:j1+1),*,iostat=ios) i
          if (ios /= 0 .or. i < 1 .or. i > 3) goto 99
          j = 0
          if (a2vec(outs(j1+3:),len(outs(j1+3:)),j,2,', ',2,-3,1,ix,j)
     . /= 1) goto 99
          if (j < 1 .or. j > nkabc(i)) then
            call info2(0,0,0,'   oops! i%i%-1j > n%i (=%i)',i,nkabc(i))
            goto 98
          endif
          iindx(i,1) = j; iindx(i,2) = j

        endif

C       Allocate array
        if (allocated(poptk)) deallocate(poptk)
        allocate(poptk(iindx(1,1):iindx(1,2),
     .                 iindx(2,1):iindx(2,2),
     .                 iindx(3,1):iindx(3,2)))

        ix(1) = iindx(1,2)-iindx(1,1)+1
        ix(2) = iindx(2,2)-iindx(2,1)+1
        ix(3) = iindx(3,2)-iindx(3,1)+1
        call info8(2,0,0,' Map data omega=%,1d, pol=%i:  %i k-pts => '//
     .      '(%i %i %i) array',omega,ipol,nkp,ix(1),ix(2),ix(3),0,0)
        if (ipol > npol) then
          call info2(0,0,0,'   oops! You have ipol>npol(=%i)',npol,0)
          goto 98
        endif

        allocate(y2(nr))
        do  iq = 1, nkp
        iq1 = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          iq1 = iq1+1
          if (i3 < i31 .or. i3 > i32) cycle
          if (i2 < i21 .or. i2 > i22) cycle
          if (i1 < i11 .or. i1 > i12) cycle
C     ... Do when this qp is symmetry-equivalent to (i1,i2,i3)
          if (ival(s_bz%ipq,iq1) == iq) then
C           q(1) = qk(1,i1,i2,i3)
C           q(2) = qk(2,i1,i2,i3)
C           q(3) = qk(3,i1,i2,i3)
C           Second derivatives for spline
            j = 1 + (iq-1) + nkp*(ipol-1)
            call cspline(pfile,pfile(1,j+1),nr,1d99,1d99,y2)
C           Spline at omega
            call csplint(pfile,pfile(1,j+1),y2,nr,omega,xx,xx2)
            poptk(i1,i2,i3) = xx / wtkp(iq) * dble(n1)*dble(n2)*dble(n3)
          endif
        enddo
        enddo
        enddo
        enddo
        deallocate(y2)

C   ... kshow branch
        else

        if (allocated(poptk)) deallocate(poptk)
        allocate(poptk(4,nkp,1))
        allocate(y2(nr))
        do  iq = 1, nkp
          call dpscop(s_bz%qp,q,3,3*iq-2,1,1d0)
C         q(1) = qk(1,i1,i2,i3)
C         q(2) = qk(2,i1,i2,i3)
C         q(3) = qk(3,i1,i2,i3)
C         Second derivatives for spline
          j = 1 + (iq-1) + nkp*(ipol-1)
          call cspline(pfile,pfile(1,j+1),nr,1d99,1d99,y2)
C         Spline at omega
          call csplint(pfile,pfile(1,j+1),y2,nr,omega,xx,xx2)
          poptk(1,iq,1) = q(1)
          poptk(2,iq,1) = q(2)
          poptk(3,iq,1) = q(3)
          poptk(4,iq,1) = xx / wtkp(iq) * dble(n1)*dble(n2)*dble(n3)
C         poptk(4,iq,1) = xx / wtkp(iq)
        enddo
        deallocate(y2)

C       Default sort
        allocate(iprm(nr))
        call dvheap(4,nkp,poptk,iprm,1d-6,1)
        if (nw >= 4) then
          call word(outs,4,j1,j2)
          if (outs(j1:j2) /= 'sort') goto 99
          if (nw >= 5) then
            call word(outs,5,j1,j2)
C           print *, outs(j1:j2)

            jj1 = j1
            call ivset(iwk,1,4,99)
            do  k = 1, 4
              jj2 = j2
              call nwordg(outs,0,', ',1,jj1,jj2)
              call tokmat(outs(jj1:jj2),(/'q1','q2','q3','ds'/),
     .          4,2,' ',i,j,.false.)
              if (i < 0) goto 99
              iwk(k) = i+1
              jj1 = jj2+2
              if (jj1 >= j2) exit
            enddo
C            print *, iwk(1:4)
          endif
C         Fill in missing parts, sanity check
          deallocate(iprm); allocate(iprm(4))
          call ivheap(1,4,iwk,iprm,101)

          i = 0
          do  k = 2, 4
            if (iwk(iprm(k)) /= 99) cycle
            iwk2 = iwk(iprm(1:k))
            call ivheap(1,k,iwk2,iwk(5),0)
C           print *, iwk2(1:k)
            do  i = 1, k
              if (iwk2(i) /= i) then
                iwk(iprm(k)) = i
                exit
              endif
            enddo
          enddo

          deallocate(iprm); allocate(iprm(nr))
          allocate(poptk2(4,nkp,1))
          do  k = 1, 4
            poptk2(k,:,1) = poptk(iwk(k),:,1)
          enddo
          call dvheap(4,nkp,poptk2,iprm,1d-6,1)
          deallocate(poptk2)

        endif

        print 344
  344   format(4x,'iq',t12,'q1',t24,'q2',t36,'q3',t50,'D(k)',t64,
     .    'BZ contr')
        sumd = 0
        do  iq1 = 1, nkp
          iq = iprm(iq1)
          write(stdo,345) iq, (poptk(i,iq,1), i=1,3), poptk(4,iq,1),
     .      poptk(4,iq,1)*wtkp(iq)/(dble(n1)*dble(n2)*dble(n3))
  345     format(i6,3f12.6,2f14.5)
          sumd = sumd +
     .      poptk(4,iq,1)*wtkp(iq)/(dble(n1)*dble(n2)*dble(n3))
        enddo
        write(stdo,"(t51,' sum =',f14.5)") sumd
        deallocate(poptk)
        deallocate(iprm)
        endif

        goto 10

C ... show
      elseif (outs(1:5) == 'show ') then
        print *, ' '
        if (havep > 0) then
          call info8(2,0,0,' file:%16p'//
     .      'energy window = (%,1d;%,1d), %i pts, %i cols',
     .      pfile(1,1),pfile(nr,1),nr,nc-1,0,0,0,0)
        endif
        if (allocated(poptk)) then
          call info8(2,0,0,' k-mapped data:%16p'//
     .      '(%i:%i,%i:%i,%i:%i)',
     .      i11,i12,i21,i22,i31,i32,0,0)
        endif

c ... save popt data
      elseif (outs(1:5) == 'save ' .or. outs(1:6) == 'savea ') then

        lbin = outs(1:5) == 'save '
        call word(outs,2,j1,j2)
        if (j2 >= j1) fn = outs(j1:j2)
        if (lbin) then
          if (j2 < j1) fn = 'poptb'
          ifi = fopna(fn,-1,4)
          rewind ifi
          call ywrm(1,' ',1,ifi,' ',pfile,0,nr,nr,nc)
          call fclose(ifi)
          call info0(2,0,0,' wrote popt to binary file '//trim(fn))
        else
          if (j2 < j1) fn = 'popta'
          ifi = fopna(fn,-1,0)
          rewind ifi
C         call ywrm(0,' ',1,ifi,'(3f14.8)',pfile,0,nr,nr,nc)
          call ywrm(0,' ',1,ifi,
     .      '(f12.6,6f12.5/(12x,6f12.5))',pfile,0,nr,nr,nc)
          call fclose(ifi)
          call info0(2,0,0,' wrote popt to ascii file '//trim(fn))
        endif

        lnsave = .false.

c ... save k-mapped data
      elseif (outs(1:6) == 'savek ' .or. outs(1:7) == 'saveka ') then

        lbin = outs(1:6) == 'savek '
        if (.not. allocated(poptk)) then
          print *, 'Sorry, no poptk available to save.'
          goto 98
        endif
        if (i31 == i32) then
          i = i12-i11+1
          j = i22-i21+1
        else
          i = (i12-i11+1)*(i22-i21+1)
          j = i32-i31+1
        endif
        call word(outs,2,j1,j2)
        if (j2 >= j1) fn = outs(j1:j2)
        if (lbin) then
          if (j2 < j1) fn = 'pkb'
          ifi = fopna(fn,-1,4)
          rewind ifi
          call ywrm(1,' ',1,ifi,' ',poptk,0,i,i,j)
          call fclose(ifi)
          call info0(2,0,0,' wrote poptk to binary file '//trim(fn))
        else
          if (j2 < j1) fn = 'pka'
          ifi = fopna(fn,-1,0)
          rewind ifi
          call ywrm(0,' ',1,ifi,'((6f12.6))',poptk,0,i,i,j)
          call fclose(ifi)
          call info0(2,0,0,' wrote poptk to ascii file '//trim(fn))
        endif

C ... abort
      elseif (outs(1:2) == 'a ') then
        call rx0('aborting popt editor ... no file written')

C ... quit
      elseif (outs(1:2) == 'q '. or. outs(1:5) == 'quit ') then
        if (lnsave) then
          print '('' popt file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) == 'y' .or. outs(1:1) == 'q'))
     .      goto 10
        endif
        call rx0('exit popt editor')

C ... help
      elseif (outs == '?') then
        print 310
        print 311
C        print 312
C        print 313
        print 314
  310   format(
     .    ' Select one of these options:'/
     .    t4,'read  [fn]'/
     .    t4,'readq [fn]'/
     .    t4,'readb [fn]'/
     .    t4,'readtk [nk]'/
     .    t4,'readbtkl [nk]'/t15,
     .    'read spectra from file.'/
     .    t15,'Use "readq" for vanilla ASCII files ',
     .    '(no expressions or % directives)'/
     .    t15,'Use "readb" for binary files'/
     .    t15,'Optional "fn" specifies file name.  If not supplied:'/
     .    t15,'"fn" defaults to ',
     .    '"popt.ext" (ASCII) or "poptb.ext" (binary)'/
     .    t15,'Use "readtk" for Kotani style, EPSnnnn.dat'/
     .    t15,'Use "readtkl" for Kotani style, no LFC, EPSnnnn.nlfc.dat'
     .    /
     .    t15,'In the latter 2 cases, specify number of k to read.'/t15,
     .    'If nk is not specified, popted expects the irr BZ.'//
     .    t4,'npol n',t15,
     .    'stipulate number of polarizations in spectra file'//
     .    t4,'scalep fac',t15,
     .    'scale spectra by fac'
     .    )

  311   format(/
C    .  t4,'...The following maps popt to the full BZ.'/
     .  t4,'kmape omega [ipol] [i1=#|i2=#|i3=#]'/
     .    t15,'Map popt(omega) to full BZ'/
     .    t15,'Optional 2nd arg is polarization (default=1)'/
     .    t15,'Optional 3rd arg selects a 2D array with single ',
     .    'element'/
     .    t15,'from 1st, 2nd, or 3rd axis (index=#)'//
     .  t4,'kshowe omega [ipol] [sort=q1,q2,q3,ds]'/
     .    t15,'Show popt(omega) for each irr k'/
     .    t15,'Optional 2nd arg is polarization (default=1)'/
     .    t15,'Optional 3rd arg sorts.  Arrange order of'/
     .    t15,'q1,q2,q3,ds ',
     .    'to sort data according 1st, 2nd, 3rd, 4th priority')

  314   format(/
     .    t4,'save  [fn]',t15,'saves popt data in binary ',
     .    'file (name="poptb" unless fn supplied)'/
     .    t4,'savea [fn]',t15,'saves popt data in ascii  ',
     .    'file (name="popta" unless fn supplied)'/
     .    t4,'savek [fn]',t15,'saves k-mapped data in binary ',
     .    'file (name="pkb" unless fn supplied)'/
     .    t4,'saveka [fn]',t15,'saves k-mapped data in ascii  ',
     .    'file (name="pka" unless fn supplied)'/
     .    t4,'q',t15,'to quit the editor'/
     .    t4,'a',t15,'to abort')

      else
        call word(outs,1,j1,j2)
        call info0(0,0,0,' popted:  "'//trim(outs(j1:j2))//
     .    '" not recognized ... enter "?" for a list of options')
      endif
      goto 10

   97 call info0(0,0,0,'%10ppopt must be read before '
     .  //'invoking this command ')
   98 call info0(0,0,0,' popted:  improper usage of '//trim(outs)//
     .  ' ... nothing done')
      goto 10
   99 call info0(0,0,0,' popted:  failed to parse arg starting at: '//
     .  trim(outs(j1:))//' ... nothing done')
      goto 10

      end

C      subroutine findirr(nirr,qirr,qp,iq)
CC- Find iq in qirr that matches qp
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nirr
CCi   qirr
CCi   qp    :k-point
CCo Outputs
CCo   iq    :point iq that matches qirr.  No match -> iq = -1
CCl Local variables
CCl         :
CCr Remarks
CCr
CCb Bugs
CCb   should call latvec to compare
CCu Updates
CCu   08 Sep 05
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nirr,iq
C      double precision qirr(3,nirr),qp(3)
CC ... Local parameters
C      integer i
C      double precision tol
C      parameter (tol=1d-5)
C
C      do  i = 1, nirr
C        if (abs(qp(1)-qirr(1,i)) < tol) then
C        if (abs(qp(2)-qirr(2,i)) < tol) then
C        if (abs(qp(3)-qirr(3,i)) < tol) then
C          iq = i
C          return
C        endif
C        endif
C        endif
C      enddo
C      iq = -1
C
C      end
