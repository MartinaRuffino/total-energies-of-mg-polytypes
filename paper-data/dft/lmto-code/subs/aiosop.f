      logical function aiosop(alabl,sop,nl,lmax,nsp,ifi)
C- File I/O for spin-orbit coupling parameters.
C ----------------------------------------------------------------------
Ci Inputs
Ci   alabl :class label
Ci   nl    :(global maximum l) + 1
Ci   lmax  :maximum l for a given site
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  sop   :spin-orbit coupling parameters (made in soprm.f)
Cio        :sop(l+1,is1,is2,i) are matrix elts between spins is1 and is2
Cio        :for quantum number l. Nine types of integrals are stored.
Cio        :(i=1) <phi so  phi> (i=2) <phi so  dot> (i=3) <dot  so dot>
Cio        :(i=4) <phi ||  phi> (i=5) <phi ||  dot> (i=6) <dot  || dot>
Cio        :(i=7) <phi Bxc phi> (i=8) <phi Bxc dot> (i=9) <dot Bxc dot>
Cio        :The first three are the terms need to make SO perturbation
Cio        :the last three are used when external field is applied.
Cr Remarks
Cr
Cu Updates
Cu    4 Apr 04 Extended  to include matrix elements for XC field
Cu   07 Feb 03 Extended  to include matrix elements for applied field
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmax,nsp,ifi
      double precision sop(0:nl-1,nsp,nsp,9)
C Local parameters
      integer l,is1,is2,ks1,ks2,j,k,lmx,nl2,nsp2,ipr,lmxf
      logical scat,pars1v
      character s*72,nam*8

      aiosop = .false.
      call getpr(ipr)
      if (ifi > 0) then
C   ... return unless file has SO category
        if (.not. scat(ifi,'SO:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        call word(s,2,is1,is2)
        nam = s(is1:is2)
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (min(nl,lmax+1) /= nl2 .and. ipr >= 10)
     .    print *, 'aiosop (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp /= nsp2 .and. ipr >= 10)
     .    print *, 'aiosop (warning) mismatch in nsp, class '//alabl
        read(ifi,'(a72)') s
        call dpzero(sop,nl*nsp*nsp*9)
        lmx = min(nl,nl2)-1
        lmxf = nl2-1
C   ... read SO parms
            do l = 1, lmxf
              do is2 = 1, min(nsp, nsp2)
                do is1 = 1, min(nsp, nsp2)
            if (l <= lmx) then
              read(ifi,*) k,ks1,ks2,(sop(l,is1,is2,j), j=1,3)
            else
              read(ifi,*) k,ks1,ks2
            endif
            if (ks1 /= is1 .or. ks2 /= is2)
     .        call rx('aiosop: spin mismatch')
                enddo
              enddo
            enddo
        read(ifi,'(a72)') s
            do l = 0, lmxf
              do is2 = 1, min(nsp, nsp2)
                do is1 = 1, min(nsp, nsp2)
            if (l <= lmx) then
              read(ifi,*) k,ks1,ks2,(sop(l,is1,is2,j), j=4,6)
            else
              read(ifi,*) k,ks1,ks2
            endif
            if (ks1 /= is1 .or. ks2 /= is2)
     .        call rxs('aiosop: spin mismatch, class ',nam)
                enddo
              enddo
            enddo
        read(ifi,'(a72)') s
            do l = 0, lmxf
              do is2 = 1, min(nsp, nsp2)
                do is1 = 1, min(nsp, nsp2)
            if (l <= lmx) then
              read(ifi,*) k,ks1,ks2,(sop(l,is1,is2,j), j=7,9)
            else
              read(ifi,*) k,ks1,ks2
            endif
            if (ks1 /= is1 .or. ks2 /= is2)
     .        call rxs('aiosop: spin mismatch, class ',nam)
                enddo
              enddo
            enddo

        if (nsp2 < nsp) then
              do j = 1, 9
                do is2 = 1, nsp
                  do is1 = 1, nsp
                    do l = 0, nl-1
                      sop(l,is1,is2,j) = sop(l,1,1,j)
                    enddo
                  enddo
                enddo
              enddo
        endif

        aiosop = .true.
        return
   18   continue
        print *, 'aiosop: (input skipped) bad syntax, class '//alabl
      else
        write (-ifi, 1) alabl, lmax+1, nsp
    1   format ('SO: ', a5, '  nl=', i1, '  nsp=', i1/
     .    '   l is js  < phi so phi>  < phi so dot>  < dot so dot>')
        do l = 1, lmax
          do is2 = 1, nsp
            do is1 = 1, nsp
              write (-ifi,4) l, is1, is2, (sop(l,is1,is2,j), j=1,3)
            enddo
          enddo
        enddo
C       write(-ifi,'(''   matrix elements of phi and phidot'')')
        write (-ifi,2)
    2   format(
     .    '   l is js  < phi || phi>  < phi || dot>  < dot || dot>')
        do l = 0, lmax
          do is2 = 1, nsp
            do is1 = 1, nsp
              write (-ifi,4) l, is1, is2, (sop(l,is1,is2,j), j=4,6)
            enddo
          enddo
        enddo
C       write(-ifi,'(''   matrix elements of phi B phidot'')')
        write (-ifi,3)
    3   format (
     .    '   l is js  <phi|Bxc|phi>  <phi|Bxc|dot>  <dot|Bxc|dot>')
        do l = 0, lmax
          do is2 = 1, nsp
            do is1 = 1, nsp
              write (-ifi,4) l, is1, is2, (sop(l,is1,is2,j), j=7,9)
            enddo
          enddo
        enddo

        aiosop = .true.
      endif
    4 format (i4,2i3,3f15.10)
      end
