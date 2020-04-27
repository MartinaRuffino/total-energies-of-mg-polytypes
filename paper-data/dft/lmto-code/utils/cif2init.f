      subroutine fmain
      implicit none
      integer ifi,jfi,i,j,fopng,i1,i2,ix(9),a2vec,nargf,spcgrpno
      character id*256, header*128, fnam*256, outnam*256
      logical rdstrn
      double precision vol,plat(3,3),px(9)
      double precision abc(3),abg(3)
      integer nspec,mxspec,nbas,ib
      parameter (mxspec=256)
      integer,allocatable :: ips(:)
      character*(8),allocatable :: slabl(:)

C     Get file name if one is specified
      if (nargf() == 1) then
        fnam = 'cif2cell.out'
      else
        call getarf(1,fnam)
      endif

CJJ, optional second argument for output name instead of "init"
C needed if more than one copy is running in a given dirn at once
      if (nargf() .eq. 3) then
        call getarf(2,outnam)
      else
        outnam='init'
      endif 

C     Read header; check that first string in 1st line is CIF2CELL
      call word(fnam,1,i1,i2)
      ifi = fopng(fnam(i1:i2),-1,1)

      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      call word(id,1,i1,i2)
      if (id(i1:i2) /= 'CIF2CELL')
     .  call rx('cif2site: first line should start with CIF2CELL')

C     Header not needed but require that it be present for now
      call findline(ifi,'Output for',id)
      header = trim(id(12:))
      print *, trim(header)

C     Most data not needed but require that it be present for now
      call findline(ifi,'Lattice parameters',id)
      call findline(ifi,'a           b           c',id)
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      i = 0
      if (a2vec(id,len(id),i,4,', ',2,-3,3,ix,abc) /= 3)
     .  call rx('failed to read a,b,c: '//id)
      call findline(ifi,'alpha        beta       gamma',id)
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      i = 0
      if (a2vec(id,len(id),i,4,', ',2,-3,3,ix,abg) /= 3)
     .  call rx('failed to read alpha,beta,gamma: '//id)

      call info8(2,0,0,
     .  ' A = %d  B = %d  C = %d'//
     .  '   alpha = %d  beta = %d  gamma = %d',
     .  abc(1),abc(2),abc(3),abg(1),abg(2),abg(3),0,0)

CC     Move file pointer
C      call findline(ifi,'INPUT CELL INFORMATION',id)
C      call findline(ifi,'Space group number',id)
C      i = index(id,':')
C      print *, id(i:)
C      if (a2vec(id,len(id),i,2,' ',1,-1,1,ix,spcgrpno) /= 1) then
C        if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
C      endif

C     Read space group number
      call findline(ifi,'OUTPUT CELL INFORMATION',id)
      call findline(ifi,'Space group number',id)
      i = index(id,':')
      if (a2vec(id,len(id),i,2,' ',1,-1,1,ix,spcgrpno) /= 1) then
        if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      endif
      call findline(ifi,'Hermann-Mauguin',id)

      call info2(2,0,0,' Space group '//trim(id(26:))//
     .  '   number %i',spcgrpno,0)

      call findline(ifi,'Bravais lattice vectors',id)
      read(ifi,*) plat
      call info5(0,0,0,
     .  ' PLAT=%3:2;11,7D%N%6f%3:2;11,7D%N%6f%3:2;11,7D',
     .  plat(1,1),plat(1,2),plat(1,3),0,0)

      call findline(ifi,'lattice coordinates',id)
      call findline(ifi,'Atom',id)

C     Count number of sites and build species table
      call info0(2,1,-1,' ... read site data')
      nbas = 0; nspec = 0
      allocate(slabl(999999),ips(999999))
      do  while (rdstrn(ifi,id,len(id),.false.))
        if (len(trim(id)) == 0) cycle
        if (id(1:4) == 'Unit') exit
C       Species is 1st argument
        call word(id,1,i1,i2)
C       Vacancies: just take compound name
        if (id(i1:i2) == '(') then
          call wordg(id,100,'{}',2,i,j)
          if (j >= 0) then
            i1 = i
            i2 = j
          endif
        endif
        i = i2
        if (a2vec(id,len(id),i,4,', ',2,-3,3,ix,px) /= 3)
     .    call rx('failed to read basis vectors: '//id)
        nbas = nbas+1

        do  i = 1, nspec
          if (id(i1:i2) /= slabl(i)) cycle
          ips(nbas) = i; goto 10 ! Species already exists; use id
        enddo
        nspec = nspec+1
        slabl(nspec) = id(i1:i2)
        ips(nbas) = nspec
   10   continue ! Re-entry for already extant species

C        call info2(0,0,0,'%5pATOM='//slabl(i)//'%14p X= %3:2;11,7D',px,
C     .    0)

      enddo

      vol = dabs(
     .      plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .      plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .      plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))

      call info2(0,0,0,' ... %i atoms, %i species',nbas,nspec)

      call info5(0,0,0,'     Cell volume = %d A^3  = %d a0^3'//
     .  '  = %d a0^3/atom',
     .  vol*abc(1)**3,vol*(abc(1)/0.529177d0)**3,
     .  vol*(abc(1)/0.529177d0)**3/nbas,0,0)

      call info2(0,0,0,' ... writing init file',nbas,nspec)

C --- Create init file ---
      jfi = fopng(outnam,-1,0)

      rewind jfi
      call awrit1('HEADER '//trim(header),' ',120,jfi,spcgrpno)
      call awrit1('LATTICE%N#       SPCGRP=%i',' ',120,jfi,spcgrpno)
      call awrit6('#%7fA=%d  B=%d  C=%d   ALPHA=%d  BETA=%d  GAMMA=%d',
     .  ' ',120,jfi,abc(1),abc(2),abc(3),abg(1),abg(2),abg(3))
      call awrit1('%% const a=%d',' ',120,jfi,abc)
      call awrit0('        ALAT={a}  UNITS=A',' ',120,jfi)
      call awrit3(
     .  '        PLAT=%3:2;11,7D%N%13f%3:2;11,7D%N%13f%3:2;11,7D',
     .  ' ',256,jfi,plat(1,1),plat(1,2),plat(1,3))
      call info5(0,0,0,
     .  '        PLAT=%3:2;11,7D%N%13f%3:2;11,7D%N%13f%3:2;11,7D',
     .  plat(1,1),plat(1,2),plat(1,3),0,0)

C     Repeat loop over site data, writing to init file
      call awrit0('SITE',' ',120,jfi)
      rewind ifi
      call findline(ifi,'Bravais lattice vectors',id)
      call findline(ifi,'lattice coordinates',id)
      call findline(ifi,'Atom',id)
      ib = 0
      do  while (rdstrn(ifi,id,len(id),.false.))
        if (len(trim(id)) == 0) cycle
        if (id(1:4) == 'Unit') exit
        call word(id,1,i1,i2)
C       Vacancies: just take compound name
        if (id(i1:i2) == '(') then
          call wordg(id,100,'{}',2,i,j)
          if (j >= 0) then
            i1 = i
            i2 = j
          endif
        endif
        ib = ib+1
        if (id(i1:i2) /= slabl(ips(ib))) call rx('bug in cif2init')

        i = i2
        if (a2vec(id,len(id),i,4,', ',2,-3,3,ix,px) /= 3)
     .    call rx('failed to read basis vectors: '//id)

        call info2(0,0,0,'%5pATOM='//slabl(ips(ib))//
     .    ' X= %3:2;11,7D',px,0)

        call awrit1('%5pATOM='//slabl(ips(ib))//
     .    ' X= %3:2;11,7D',' ',120,jfi,px)


      enddo

      end
      subroutine errm(id)
      character id*(*)

      if (id /= ' ') print *, 'failed at line '//trim(id)
      call rx0('usage:  cif2init [filename]')

      end
      subroutine findline(ifi,strn,id)
      implicit none
      integer ifi
      character strn*(*),id*(*)
      logical rdstrn

      id = ' '
      do
      if (.not. rdstrn(ifi,id,len(id),.false.)) ! New line
     .    call rxs(' failed to find string containing : ',strn)
      if (index(id,strn) /= 0) return
      enddo

C      call word(id,1,i1,i2)
C      if (id(i1:i2) /= arg2) cycle
C      if (narg >= 2) then
C        call word(id,2,i1,i2)
C        if (id(i1:i2) /= arg2) cycle
C      endif
C      if (narg >= 3) then
C        call word(id,3,i1,i2)
C        if (id(i1:i2) /= arg3) cycle
C      endif
C

      end

