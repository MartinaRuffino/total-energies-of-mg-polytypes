      subroutine suclst(nsitmx,nbas,nsp,s_site,s_spec,
     .  clsopt,isite,iclsl,iclsn,nsites)
C- Set up list of sites for core level spectra
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs: nsitmx (dimension of isite,iclsl,iclsn: set in bndfp)
Ci         nbas,clsopt (command line options for parsing)
Co Outputs:
Co         isite,iclsl,iclsn (list of site,core-l and core-n)
Co         nsites (total number of site triples returned)
Cr Remarks
Cr parses the string
Cr  --cls[.ib,l,n[.ib,l,n][.lst=list,l,n[.lst=list,l,n]][.fn=file]]
Cr       where '.' is any separator
Cr                 except a ':' or a ',' if lst= is used
Cr                 except a ',' if ib,l,n is used
Cr       ib,l,n is a triple of site,core-l,core-n
Cr       list is a list of sites having l and n in common
Cr       file is a file of nsites records of triples
Cr example
Cr       --cls/2,1,2/4,2,3/lst=3,5:7,0,1
Cr  produces
Cr         site  l  n
Cr           2   1  2
Cr           4   2  3
Cr           3   0  1
Cr           5   0  1
Cr           6   0  1
Cr           7   0  1
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nsitmx
      integer isite(nsitmx),iclsl(nsitmx),iclsn(nsitmx),nsites,nbas,nsp
      character*(*) clsopt
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer i,j,j1,j2,iprmin,nblst,ib,is,l,n,ifi,p,q
      integer iprint,lgunit,fopna
      character clabl*8,fn*16,dc*1,outstr*64,cc2*2
      logical a2bin

      iprmin = 10

      dc = clsopt(1:1)
      if (dc == ' ') call rx('suclst: found --cls but no site')
      nsites = 0
      j2 = 0
      fn =  ' '
C ... Return here to resume parsing for arguments
    1 continue
        j2 = j2 + 1
        if (clsopt(j2:j2) == dc) goto 1
        j1 = min(len(clsopt),j2)
        call nwordg(clsopt,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (clsopt(j1:j1+2) == 'fn=')  then
            if (j1+3 > j2) call rx('suclst: bad file name')
            fn = clsopt(j1+3:j2)
C       ... read in data from file
            ifi = fopna(fn,-1,1)
            i = 0
    2       continue
            read (ifi,*,end=3,err=3) isite(i+1),iclsl(i+1),iclsn(i+1)
            i = i + 1
            if (i > nsitmx) call rxi('bndfp needs nsitmx',nsites)
            goto 2
    3       continue
            call fclose(ifi)
            nsites = i
            goto 4
          elseif (clsopt(j1:j1+3) == 'lst=') then
            if (j1+4 > j2) call rx('suclst: bad list')
C   ... expect ,l,n at end of word; move pointer back to end of list
            j = j2 - 4
            call mkils0(clsopt(j1+4:j),nblst,i)
            call mkilst(clsopt(j1+4:j),nblst,isite(1+nsites))
            j = j + 1
            if (clsopt(j:j) /= ',') goto 5
            i = 0
            j = j + 1
            cc2 = clsopt(j:j)//' '
            if (.not. a2bin(cc2,l,2,0,' ',i,-1)) goto 5
            j = j + 1
            if (clsopt(j:j) /= ',') goto 5
            i = 0
            j = j + 1
            cc2 = clsopt(j:j)//' '
            if (.not. a2bin(cc2,n,2,0,' ',i,-1)) goto 5
            if (j /= j2) call rx('bug in suclst')
            nsites = nsites + nblst
            if (nsites > nsitmx) call rxi('bndfp needs nsitmx',nsites)
            call icopy(nblst,l,0,iclsl(nsites+1-nblst),1)
            call icopy(nblst,n,0,iclsn(nsites+1-nblst),1)
          else
            j = j1
            i = 0
            cc2 = clsopt(j:j)//' '
            if (.not. a2bin(cc2,ib,2,0,' ',i,-1)) goto 6
            j = j + 1
            if (clsopt(j:j) /= ',') goto 6
            j = j + 1
            i = 0
            cc2 = clsopt(j:j)//' '
            if (.not. a2bin(cc2,l,2,0,' ',i,-1)) goto 6
            j = j + 1
            if (clsopt(j:j) /= ',') goto 6
            j = j + 1
            i = 0
            cc2 = clsopt(j:j)//' '
            if (.not. a2bin(cc2,n,2,0,' ',i,-1)) goto 6
            if (j /= j2) call rx('bug in suclst')
            nsites = nsites + 1
            if (nsites > nsitmx) call rxi('bndfp needs nsitmx',nsites)
            isite(nsites) = ib
            iclsl(nsites) = l
            iclsn(nsites) = n
          endif
          goto 1
        endif
    4   continue
      if (nsites == 0) call rx('suclst: found --cls but no site')
C --- print table of sites, l and n
      if (nsites > 0 .and. iprint() >= iprmin) then
        if (nsp == 1) then
          call awrit0('%N suclst: set up table of site, l and n for CLS'
     .      //'%N   site   label     l  n   DOS channels',
     .      ' ',256,lgunit(1))
        else
          call awrit0('%N suclst: set up table of site, l and n for CLS'
     .      //'%N   site   label     l  n   DOS channels up(down)',
     .      ' ',256,lgunit(1))
        endif
        j = 1
        do  i = 1, nsites
          ib = isite(i)
          if (ib > nbas) call rx('suclst: ib > nbas')
          is = s_site(ib)%spec
          clabl = s_spec(is)%name
          if (nsp == 1) then
            call awrit6('%:3,3i     '//clabl//'%,2i%:1,2i%3f%i,%i,%i',
     .           ' ',128,lgunit(1),isite(i),iclsl(i),iclsn(i),j,j+1,j+2)
          else
            call awrit3('%:3,3i     '//clabl//'%,2i%:1,2i',
     .                  outstr,128,0,isite(i),iclsl(i),iclsn(i))
            call strip(outstr,p,q)
            call awrit6(outstr(1:q)//'%3f%i(%i),%i(%i),%i(%i)',' ',128,
     .                  lgunit(1),j,j+1,j+2,j+3,j+4,j+5)
          endif
          j = j + 3*nsp
        enddo
      endif
      return
    5 continue
C --- bad list ERROR
      call fexit(-1,119,
     .       'suclst: bad list, expecting --cls'//dc//'lst=list,l,n',0)
    6 continue
C --- bad triple ERROR
      call fexit(-1,119,
     .       'suclst: bad triple, expecting --cls'//dc//'site,l,n',0)
      end

