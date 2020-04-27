      subroutine fmain
      implicit none

      integer ifi,i,j,fopng,i1,i2,ix(9),a2vec,nargf
      character id*256, fnam*256  ! header*128, sgnam*24
      logical rdstrn
      logical lcart  ! T for basis in Cartesian coordinates
      double precision vol,vol0,alat,plat(3,3),px(9),posi(3)
      integer j1,j2,nspec,nbas,isw
      integer,allocatable :: ibas(:)
      character*(8),allocatable :: slabl(:)
      real(8),allocatable :: pos(:,:)
      logical speclabels ! T => POSCAR file specifies species before sites
      integer awrite,ib,jfi

C     Get file name if one is specified
      if (nargf() == 1) then
        fnam = 'POSCAR'
      else
        call getarf(1,fnam)
      endif

C     Read header
      call word(fnam,1,i1,i2)
      ifi = fopng(fnam(i1:i2),-1,1)
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)

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

C     alat is lattice constant if > 0, but volume if < 0
      vol = -alat
C     Read plat
      read(ifi,*) plat
      vol0 =dabs(
     .      plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .      plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .      plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))
C     Get alat if volume was read
      if (vol > 0) alat = (vol/vol0)**(1d0/3d0)

C     Convert to a.u.
      alat = alat / 0.529177d0

      vol = alat**3*dabs(
     .      plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .      plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .      plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))

      call info2(2,0,0,' lattice constant %d A = %d a.u.',
     .  alat*0.529177d0,alat)
      call info5(2,0,0,
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
        i = 0
        if (a2vec(id,len(id),i,2,', ',2,-3,nspec,ix,ibas) /= nspec)
     .    call errm(id)
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

C     Read whether basis is stored in Cartesian coordinates or multiples of plat
      if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
      call locase(id)
      if (id == 'selective dynamics') then  ! Assume it has no effect; read again
        if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)
        call locase(id)
      endif
      lcart = id == 'cartesian'
      if (.not. lcart .and. id /= 'direct') call rx(" basis must be in 'cartesian' or 'direct' mode")

C     Count the number of basis sites
      nbas = 0
      do  i = 1, nspec
        do  j = 1, ibas(i)
          nbas = nbas+1
        enddo
      enddo

      allocate(pos(3,nbas))
      call info2(2,0,0,'# %i atoms, %i species',nbas,nspec)

C ... Case labels not given a priori
      if (.not. speclabels) then

        call info0(2,1,0,' Writing site data to file site')

        jfi = fopng('site',-1,0)
        rewind jfi
        i = awrite('%% site-data vn=3.0 xpos fast io=15 '//
     .    'nbas=%i alat=%d plat=%9:1;8,1d',' ',180,jfi,
     .    nbas,alat,plat,px,px,px,px,px)

C       Read through all sites
        ib = 0
        do  ib = 1, nbas

          if (.not. rdstrn(ifi,id,len(id),.false.)) call errm(id)

C         Species is 4th argument
          call word(id,4,i1,i2)
C         Vacancies: just take compound name
          if (id(i1:i2) == '(') then
            call wordg(id,100,'{}',2,i,j)
            if (j >= 0) then
              i1 = i
              i2 = j
            endif
          endif

          i = 0
          if (a2vec(id,len(id),i,4,', ',2,-3,3,ix,posi) /= 3)
     .      call rx('failed to read basis vectors: '//id)

          call info2(2,0,0,'%5pATOM='//id(i1:i2)//'%14p X= %3:2;11,7D',posi,0)

C         pos(1:3,ib) = posi(1:3)
          write(jfi,321) id(i1:i2),posi
        enddo

        stop

      else ! Species labels already given

        call info0(2,0,0,'SITE')

C       Read through all species, and sites within species
        nbas = 0
        do  i = 1, nspec
          do  j = 1, ibas(i)
            nbas = nbas+1
            read(ifi,*) posi
            pos(1:3,nbas) = posi(1:3)
            call info2(2,0,0,'%5pATOM='//trim(slabl(i))//
     .        '%14p X= %3:2;11,7D',posi,0)
          enddo
        enddo

      endif

      call info2(2,1,0,
     .  " Writing site data to file 'site', basis in %?#n#Cartesian coordinates#multiples of PLAT#",
     .  isw(lcart),0)

      ifi = fopng('site',-1,0)
      rewind ifi
      i = awrite('%% site-data vn=3.0 %?#n##xpos #fast io=15 '//
     .  'nbas=%i alat=%d plat=%9:1;8,1d',' ',180,ifi,isw(lcart),
     .  nbas,alat,plat,px,px,px,px)

      nbas = 0
      do  i = 1, nspec
        do  j = 1, ibas(i)
          nbas = nbas+1
          posi(1:3) = pos(1:3,nbas)
          write(ifi,321) trim(slabl(i)),posi
  321     format(1x,a4,3f17.10)
        enddo
      enddo

      end

      subroutine errm(id)
      character id*(*)

      if (id /= ' ') print *, 'failed at line '//trim(id)
      call rx0('usage:  poscar2init [filename]')

      end
