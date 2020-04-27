      subroutine GWinputprodbas(nl,nat,nphimx,ncoremx,lcuta,
     .  nrphiv,nrphic,noccv,nunoccv,noccc,nunoccc,ncwf,ncwf2,tolbas)
C- Read product basis from GWinput
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1 --- here a dimensioning parameter
Ci   nat :number of sites with augmentation spheres
Ci   nphimx:dimensioning parameter : max number of partial waves in one l channel
Ci  ncoremx:dimensioning parameter : max number of core levels in one site
Co Outputs
Co   nrphiv:number of partial waves for this l and site
Co         :nrphiv named nindxv in old GW code
Co   nrphic:number of core levels for this l and site
Co         :nrphic named nindxc in old GW code
Co   noccv :noccv(l,ib) = 1 if to include as 1st fn in product basis, 0 otherwise
Co  nunoccv:nunoccv(l,ib) = 1 if to include as 2nd fn in product basis, 0 otherwise
Co   noccc :noccc(l,i,ib) = 1 if to include this core as 1st fn in product basis, 0 otherwise
Co  nunoccc:nunoccc(l,i,ib) = 1 if to include this core as 2nd fn in product basis, 0 otherwise
Co   ncwf  :ncwf(l,i,ib) 1 => include this core level in susceptibility
Co   ncwf2 :ncwf2(l,i,ib) 1 => include this core level in sigma
Co   tolbas:product basis tolerance, by l
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Dec 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nat,nphimx,ncoremx,lcuta(nat)
      integer nrphiv(0:nl-1,nat),nrphic(0:nl-1,nat),
     .        noccv(0:nl-1,nphimx,nat),noccc(0:nl-1,ncoremx,nat),
     .        nunoccv(0:nl-1,nphimx,nat),nunoccc(0:nl-1,ncoremx,nat),
     .        ncwf(0:nl-1,ncoremx,nat),ncwf2(0:nl-1,ncoremx,nat)
      double precision tolbas(0:2*(nl-1))
C ... Local parameters
      integer,parameter:: NULLI = -99999
      integer ifi,i,j,k,l,m,n,iv(max(10,2*(nl-1))),ib
      procedure(integer) GWinputcat,fopng,a2vec
      character*512 strn

      call info2(30,1,0,' ... reading product basis specs from GWinput.  nphimx=%i  nat=%i',nphimx,nat)
      ifi = GWinputcat(11,'<PRODUCT_BASIS>',1)
      ifi = fopng('GWinput',-1,0)
      read(ifi,*) ! First line is a commment
      read(ifi,"(a)") strn  ! read tolbas
      m = 0
      k = a2vec(strn,len_trim(strn),m,4,' ',1,2,size(tolbas),iv,tolbas)
      if (k < 0) call rx('GWinputprodbas: failed to read tolbas')
      tolbas(k:) = tolbas(k-1)
      call info2(30,0,0,'%5feigenvalue cutoff to reduce rank (by l): %n:1;3e',k,tolbas)

      call iinit(noccv,size(noccv))
      call iinit(nunoccv,size(nunoccv))
      call iinit(noccc,size(noccc))
      call iinit(nunoccc,size(nunoccc))
      call iinit(ncwf,size(ncwf))
      call iinit(ncwf2,size(ncwf2))

      read(ifi,*) ! This line is a commment
      read(ifi,*) lcuta
      call info2(30,0,0,'%5flmax (by site): %n:1i',nat,lcuta)

C --- Read number of core and valence partial waves ---
      read(ifi,*) ! This line is commment 'atom   l  nnvv  ...'
      do  ib = 1, nat
        do  l = 0, nl-1
          read(ifi,*) i,k,nrphiv(l,ib),nrphic(l,ib)
          if (k /= l) call rx2('GWinputprodbas: wrong l, site %i, l=%i',ib,l)
        end do
      end do

C --- Valence product basis ---
      read(ifi,*) ! This line is commment 'atom   l    n  occ  unocc ...'
      do  ib = 1, nat
        do  l = 0, nl-1
          do  n = 1, nrphiv(l,ib)
            read(ifi,*) j,k,m,noccv(l,n,ib),nunoccv(l,n,ib)
            if (k /= l) call rx('GWinputprodbas: wrong l, site %i, l=%i',ib,l)
            if (m /= n) call rx('GWinputprodbas: wrong n, site %i, l=%i',ib,n)
          end do
        end do
      end do

C --- Core product basis ---
      read(ifi,*) ! This line is commment 'atom   l    n  occ unocc   ForX0'
      do  ib = 1, nat
        do  l = 0, nl-1
          do  n = 1, nrphic(l,ib)
            read(ifi,*) j,k,m,noccc(l,n,ib),nunoccc(l,n,ib),ncwf(l,n,ib),ncwf2(l,n,ib) !ncwf2 is for Sigma calcuation
            if (k /= l) call rx('GWinputprodbas: wrong l, site %i, l=%i',ib,l)
            if (m /= n) call rx('GWinputprodbas: wrong n, site %i, l=%i',ib,n)
          enddo
        enddo
      enddo
      call fclose(ifi)

      end

      integer function GWinputcat(opt,cat_tags,ntags)
C- Move file pointer to start of category
C ----------------------------------------------------------------------
Ci Inputs
Ci  opt    : 1s digit
Ci         : 0 do not rewind GWinput before read
Ci         : 1 rewind GWinput before read
Ci         : 10s digit
Ci         : 0 return regardless of whether category was found.
Ci         : 1 fatal error if category not found
Ci cat_tags: string (or strings labelling categories
Ci  ntags  :Number of categories to seek
Co Outputs
Cu Updates
Cu   12 Dec 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ntags
      character*(15) cat_tags(ntags)
C ... Local parameters
      integer i1,i2,j1,j2,ifi
      character strn*512
      procedure(logical) rdstrn ! cmdopt
      procedure(integer) fopng

      ifi = fopng('GWinput',-1,1)
      if (mod(opt,10) == 1) rewind ifi

      do  while (rdstrn(ifi,strn,len(strn),.false.))
        if (strn(1:1) == '<') then
          call word(strn,1,i1,i2)
          call tokmat(strn(i1:i2),cat_tags,ntags,len(cat_tags(1)),' ',j1,j2,.false.)
          GWinputcat = j1
          if (j1 >= 0) return
        endif
      enddo

      if (mod(opt/10,10) == 1) call rx('GWinputcat cannot find category '//cat_tags(1))
      GWinputcat = -1           ! no category found

      end
