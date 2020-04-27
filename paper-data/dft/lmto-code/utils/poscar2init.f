      subroutine fmain
      implicit none
      integer ifi,i,j,i1,i2,ix(9)
      character id*256, header*128, sgnam*24, fnam*256
      logical lok
      double precision vol,vol0,alat,plat(3,3),px(9),posi(3),scale
      integer j1,j2,nspec,mxspec,nbas
      parameter (mxspec=256)
      integer,allocatable :: ibas(:)
      character*(8),allocatable :: slabl(:)
      real(8),allocatable :: pos(:,:)
      logical speclabels ! T => POSCAR file specifies species before sites
      procedure(logical) :: rdstrn,parstr
      procedure(integer) :: a2vec,nargf,fopng,wordsw

      i = 1
      scale = 1
      call getarf(i,fnam)
      j = wordsw(fnam,'','-scale','= ',j1);
      if (j /= 0) then
        if (a2vec(fnam,len_trim(fnam),j1,4,' ',2,2,1,ix,scale) /= 1) call rx('failed to parse '//fnam)
        i = i+1
      endif

C     Get file name if one is specified
      if (nargf() == i) then
        fnam = 'POSCAR'
      else
        call getarf(i,fnam)
      endif

C     Read header
      call word(fnam,1,i1,i2)
      ifi = fopng(fnam(i1:i2),-1,1)
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)

C --- Case header exists beginning with 'id'.  Assume uniform style ---
      sgnam = ' '
      if (id(1:2) == 'id') then

C       Get label
        i = 0
        lok = parstr(id,'formula=',len(id),8,'=',i,j)
        i1 = i+9+1
        call nwordg(id,101,']',1,i1,i2)
        header = id(i1:i2)

        call info0(0,0,0,'HEADER  '//id(i1:i2)//'%a')
        call info0(0,0,0,'CMD     --fixpos:tol=1e-3')

C       Get symmetry group name
        lok = parstr(id,'sg_name=',len(id),8,'=',i,j)
        i1 = i+9+1
        call nwordg(id,101,']',1,i1,i2)
        j = 0
        do  i = i1, i2
          if (id(i:i) == ' ') then
          else
            j = j+1
            sgnam(j:j) = id(i:i)
          endif
        enddo
        if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      endif

C     Check whether first string in line is a number
C     If not assume this line is a comment line; read the next one
      i = 0
      if (a2vec(id,len(id),i,4,' ',1,-1,1,ix,px) /= 1) then
        if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      endif

C     We have definitely moved past the header.
C     This line should contain the lattice constant or the volume
      i = 0
      if (a2vec(id,len(id),i,4,' ',1,-1,1,ix,alat) /= 1) call errm(id)
      alat = alat * scale

C     alat is lattice constant if > 0, but volume if < 0
      vol = -alat
C     Read plat
      read(ifi,*) plat
      plat = plat / scale
      vol0 =dabs(
     .      plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .      plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .      plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))
C     Get alat if volume was read
      if (vol > 0) alat = (vol/vol0)**(1d0/3d0)

      vol = alat**3*dabs(
     .      plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .      plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .      plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))

      call info0(0,0,0,'LATTICE%N#       SPCGRP='//sgnam//'%a')
      call info2(0,0,0,'        ALAT=%d  UNITS=A',alat,0)
      call info5(0,0,0,
     .  '        PLAT=%3:2;11,7D%N%13f%3:2;11,7D%N%13f%3:2;11,7D',
     .  plat(1,1),plat(1,2),plat(1,3),0,0)

C     Next line may be a sequence of numbers, or a sequence of species
C     In either case assume it specifies species
C     so number of words = number of species
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      call words(id,nspec)
      allocate(slabl(nspec),ibas(nspec))

C     Case first entry is a number => assume number of sites/species
      i = 0
      if (a2vec(id,len(id),i,4,' ',1,-1,1,ix,px) == 1) then
        backspace ifi
        read(ifi,*) ibas(1:nspec)
        speclabels = .false.

C     Otherwise, assume a sequence of species
      else
        j1 = 1
        do  i = 1, nspec
          call wordg(id,1,' ',i,j1,j2)
          slabl(i) = id(j1:j2)
        enddo
C       Next line contains number of sites per species
        read(ifi,*) ibas(1:nspec)
        speclabels = .true.
      endif

C     Read how site data is stored
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      call locase(id)
      if (id /= 'direct') call rx(" basis must be in 'direct' mode")

C ... Case labels not given a priori
      if (.not. speclabels) then
        call info0(0,0,0,'SITE')

C       Read until sites are exhausted
        do while (rdstrn(ifi,id,len(id),.false.))
C       Species is 4th argument
        call word(id,4,i1,i2)
C       Vacancies: just take compound name
        if (id(i1:i2) == '(') then
          call wordg(id,100,'{}',2,i,j)
          if (j >= 0) then
            i1 = i
            i2 = j
          endif
        endif
        i = 0
        if (a2vec(id,len(id),i,4,', ',2,-3,3,ix,px) /= 3)
     .    call rx('failed to read basis vectors: '//id)

        call info2(0,0,0,'%5pATOM='//id(i1:i2)//'%14p X= %3:2;11,7D',px,
     .    0)

      enddo

      else ! Species labels already given
        nbas = 0
        do  i = 1, nspec
          do  j = 1, ibas(i)
            nbas = nbas+1
          enddo
        enddo

        allocate(pos(3,nbas))
        call info2(0,0,0,'# %i atoms, %i species',nbas,nspec)
        call info0(0,0,0,'SITE')

C       Read through all species, and sites within species
        nbas = 0
        do  i = 1, nspec
          do  j = 1, ibas(i)
            nbas = nbas+1
            read(ifi,*) posi
C           pos(1:3,nbas) = posi(1:3)
            call info2(0,0,0,'%5pATOM='//trim(slabl(i))//
     .        '%14p X= %3:2;11,7D',posi,0)
          enddo
        enddo

      endif
      end
      subroutine errm(id)
      character id*(*)

      if (id /= ' ') print *, 'failed at line '//trim(id)
      call rx0('usage:  poscar2init [filename]')

      end
