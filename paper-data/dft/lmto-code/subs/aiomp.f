      logical function aiomp(alabl,pmpol,nl,lmxv,nsp,ifi)
C- File I/O for ASA multipole moments integrals
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   alabl : class label
Ci   nl    :(global maximum l) + 1
Ci   lmxv  :l-cutoff for potential
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  pmpol :intra-atomic density-density response matrix
Cr Remarks
Cr Bugs
Cr   Input not checked when file nsp mismatches passed nsp
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*8 alabl
      integer nl,lmxv,nsp,ifi
      double precision pmpol(nl,nl,0:lmxv,3,nsp)
C ... Local parameters
      double precision x2(10)
      integer i,n,ip,k,a2vec,ix2(10),ipr,nl2,nsp2,l1,l2,lf,idamax
      logical scat,rdstrn,pars1v
      character*100 s

      aiomp = .false.
      call getpr(ipr)
C --- File READ ---
      if (ifi > 0) then
C   ... return unless file has PMPOL: category
        if (.not. scat(ifi,'PMPOL:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (nl /= nl2 .and. ipr >= 10)
     .    print *, 'aiomp (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp /= nsp2 .and. ipr >= 10)
     .    print *, 'aiomp (warning) mismatch in nsp, class '//alabl
C   ... For now, abort
        if (nsp /= nsp2 .or. nl /= nl2) goto 18
        call dpzero(pmpol,nl*nl*(lmxv+1)*3*nsp)
C   ... read pmpol
        n = min(nsp,nsp2)
        do i = 1, n
        do k = 1, 3
          do l1 = 0, nl-1
          do l2 = 0, nl-1

   16     if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
          if (s == ' ') goto 16
          ip = 0
C    ...  Abort if fail to read lmxv+1 numbers
          ip = a2vec(s,len(s),ip,4,' ',1,-2,-(lmxv+3),ix2,x2)
          if (ip /= lmxv+3) call rxs('AIOMP: failed to parse ',s)
          call dcopy(lmxv+1,x2(3),1,pmpol(l1+1,l2+1,0,k,i),nl**2)
          enddo
          enddo
        enddo
        enddo
        aiomp = .true.
        return
   18   continue
        print *, 'aiomp: (input skipped) bad syntax, class '//alabl
        return

C --- File WRITE ---
      else
        write(-ifi,21) alabl, nl, nsp
        do i = 1, nsp
        do k = 1, 3
          do l1 = 0, nl-1
          do l2 = 0, nl-1
C           x2(1) = damax(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl)
            ip = idamax(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl)
            x2(1) = pmpol(l1+1,l2+1,ip-1,k,i)
            if (x2(1) > 9999d0 .or. x2(1) < -999d0) then
              call dcopy(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl,x2,1)
C             call awrit4('%,4i%,4i %n:1g',' ',100,6,l1,l2,lmxv+1,x2)
              call awrit4('%,4i%,4i %n:1;11F',' ',100,-ifi,l1,l2,lmxv+1,
     .          x2)
            else
              write(-ifi,333) l1,l2,(pmpol(l1+1,l2+1,lf,k,i),lf=0,lmxv)
            endif
              enddo
            enddo
  333     format(2i4,10f12.6)
          if (i /= nsp .or. k /= 3) write(-ifi,333)

          enddo
        enddo
        aiomp = .true.
      endif
   21 format('PMPOL: ',a4,'  nl=',i1,'  nsp=',i1)

      end