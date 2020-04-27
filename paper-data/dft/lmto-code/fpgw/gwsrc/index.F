      integer function nocc (ek,ef,locc,
     d                       n)
C- Counts the no. occupied/unoccupied states

c ek   = eigenvalues
c ef   = fermi level
c locc = .true.  ==> occupied
c      = .false. ==> unoccupied
c n    = no. states
      implicit none
      logical, intent(in) :: locc
      integer, intent(in) :: n
      real(8), intent(in) :: ek(n), ef
      integer :: i
      nocc       = 0
      if (locc) then
        do       i = 1,n
          if (ek(i) <= ef) nocc = nocc + 1
        end do
      else
        do       i = 1,n
          if (ek(i) > ef) nocc = nocc + 1
        end do
      endif
      return
      end function nocc

      subroutine maxdim (nocc,nunocc,nindx,nl,nn,nclass,
     o nprodx,nlnx,nlnmx,nlnax)
c largest number of product basis, (l,n) and (l,n,m)
      implicit real*8(a-h,o-z)
      dimension nocc(nl*nn,nclass),
     i          nunocc(nl*nn,nclass),
     i          nindx(0:nl-1,nclass)
      nprodx     = 0
      nlnx       = 0
      nlnmx      = 0
      nlnax      = 0
      do      ic = 1,nclass
        nprod      = nallow (nocc(1,ic),nunocc(1,ic),
     .                     nindx(0,ic),nl,nn)
        nln        = nofln(nindx(0,ic),nl)
        nlnm       = noflnm(nindx(0,ic),nl)
        nlna       = nalwln (nocc(1,ic),nunocc(1,ic),
     .                     nindx(0,ic),nl,nn)
        if(nprod > nprodx) nprodx = nprod
        if(nln > nlnx)   nlnx   = nln
        if(nlnm > nlnmx)  nlnmx  = nlnm
        if(nlna > nlnax)  nlnax  = nlna
      end do
      return
      end

      integer function noccx1 (ekt,nk,nt,ef)
c 92.02.27
c find the highest number of occupied states
c ekt   = eigenvalues for all k-points and states
c nk    = number of k-points in the 1st BZ
c nt    = number of states
c ef    = fermi level
c noccx1= highest number of occupied states
      implicit real*8(a-h,o-z)
      dimension ekt(nt,nk)
      logical:: allocc
!! Find the highest number of occupied states ---
!! mar2016. (Seungwoo)
!!   Previous version works,But this is necessary when all bands occuied (for model calculaitons). -->I had bug -->see
!! aug2017 below (okumura test)
!!   bugfix (probably affects little only for metallic case). For single band model, it stops since previous version givs noccx1=0.
!!
      noccx = 0
      allocc=.true.
      do   k  = 1,nk !k index
        do it = 1,nt !band index
          if(ekt(it,k) > ef) then
             allocc=.false.
            if(it-1> noccx) noccx = it-1 !bugfix aug2017takao => previous verion check noccx only at it=nt because of bug.
            exit
          endif
        enddo
      enddo
      if(allocc) noccx=nt   !completely occupied case. aug2017takao
      noccx1 = noccx
      return
      end

      subroutine nolnma  (nindx,nl,nclass,
     o nlnm )
c 92.jan.07
c number of l,n,m for all classes
      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass)
      dimension nlnm(nclass)
      do     ic = 1,nclass
        noflnm    = 0
        do 1    l = 0,nl-1
          noflnm    = noflnm + nindx(l,ic)*(2*l+1)
    1   continue
        nlnm(ic)  = noflnm
      end do
      return
      end

      integer function maxocc2 (nspin,ef, nband, qbz,nqbz)
c maximum no. occupied states
      use m_readeigen, only: readeval
      implicit none
      integer(4):: nspin,nqbz,nband,noccx,is,iq,noccxt,noccx1
      real(8) :: qbz(3,nqbz),ef
      real(8),allocatable :: ekt(:,:)
      allocate( ekt(nband, nqbz ) )
      noccx      = 0
      do  is = 1,nspin
        do iq = 1,nqbz
          call readeval(qbz(:,iq),is, ekt(1,iq))
        enddo
        noccxt  = noccx1 (ekt,nqbz,nband,ef)
        if (noccxt > noccx) noccx = noccxt
      enddo
      maxocc2 = noccx
      deallocate( ekt )
      end
c------------------------------------------------------------------


C      subroutine idxlmto (nindxv,iclass,
C     d                    nl,nnv,nlmto,natom,nclass,
C     o il,in,im,ilmto)
C
Cc 92.03.14
Cc 92.10.13 from idxlmto (see above) modified for all atoms
Cc indexing of LMTO basis functions
Cc follows that in TB-LMTO program
Cc
Cc l    0     1                 2
Cc n    1  2  1        2        1               2
Cc m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
Cc ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
Cc
Cc il(ind),in(ind),im(ind) = l,n,m
Cc ilmto(n,lm,atom) = index of lmto
Cc lm = l*l + l + m + 1
C
C      implicit real*8(a-h,o-z)
C      dimension nindxv(0:nl-1,nclass),iclass(natom),
C     o          ilmto(nnv,nl*nl,natom),in(nlmto),il(nlmto),im(nlmto)
C
C      ind       = 0
C      do 1   ia = 1,natom
C        ic        = iclass(ia)
C      do 1    l = 0,nl-1
C        l2        = l*l
C      do 1    n = 1,nindxv(l,ic)
C      do 1    m = -l,l
C        ind       = ind + 1
C        lm        = l2 + l + m + 1
C        il(ind)   = l
C        in(ind)   = n
C        im(ind)   = m
C        ilmto(n,lm,ia) = ind
C    1 continue
CCstop2rx 2013.08.09 kino      if (ind /= nlmto) stop 'idxlmto: wrong no. l,n,m'
C      if (ind /= nlmto) call rx( 'idxlmto: wrong no. l,n,m')
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine idxlmt1 (nindxv,iclass,
C     d                    nl,nnv,nlmto,natom,nclass,
C     o il,in,im,ilmto,ialmto)
C
Cc 93.08.?? from idxlmto (see above) modified to include atom index
Cc indexing of LMTO basis functions
Cc follows that in TB-LMTO program
Cc
Cc l    0     1                 2
Cc n    1  2  1        2        1               2
Cc m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
Cc ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
Cc
Cc il(ind),in(ind),im(ind) = l,n,m
Cc ilmto(n,lm,atom) = index of lmto
Cc ialmto(ind) = ia
Cc lm = l*l + l + m + 1
C
C      implicit real*8(a-h,o-z)
C      dimension nindxv(0:nl-1,nclass),iclass(natom),
C     o          ilmto(nnv,nl*nl,natom),in(nlmto),il(nlmto),im(nlmto),
C     o          ialmto(nlmto)
C
C      ind       = 0
C      do 1   ia = 1,natom
C        ic        = iclass(ia)
C      do 1    l = 0,nl-1
C        l2        = l*l
C      do 1    n = 1,nindxv(l,ic)
C      do 1    m = -l,l
C        ind       = ind + 1
C        lm        = l2 + l + m + 1
C        il(ind)   = l
C        in(ind)   = n
C        im(ind)   = m
C        ialmto(ind) = ia
C        ilmto(n,lm,ia) = ind
C    1 continue
CCstop2rx 2013.08.09 kino      if (ind /= nlmto) stop 'idxlmto: wrong no. l,n,m'
C      if (ind /= nlmto) call rx( 'idxlmto: wrong no. l,n,m')
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine idxlna (nindx,
C     d                   nl,nn,nclass,nlnx,
C     o                   ilr,inr,ilnr,nlnr)
C
Cc 92.jan.13
Cc index of l,n for radial phi(l,n) for all classes
C
Cc nlnx = max. number of l,n  from maxdim.f
C
Cc ilr(index) = l
Cc inr(index) = n
Cc ilnr(l,n)  = index
Cc nlnr       = no. l,n = max(index)
C
C      implicit real*8(a-h,o-z)
C      dimension nindx(0:nl-1,nclass)
C      dimension ilr(nlnx,nclass),inr(nlnx,nclass),
C     o          ilnr(0:nl-1,nn,nclass),nlnr(nclass)
C
C      do      ic = 1,nclass
C        i          = 0
C        do       l = 0,nl-1
C          do       n = 1,nindx(l,ic)
C            i          = i + 1
C            ilr(i,ic)   = l
C            inr(i,ic)   = n
C            ilnr(l,n,ic)= i
C          end do
C        end do
CCstop2rx 2013.08.09 kino        if(i > nlnx)stop 'idxlna: nln > nlnx'
C        if(i > nlnx)call rx( 'idxlna: nln > nlnx')
C        nlnr(ic)    = i
C      end do
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine idxlnp (nocc,nunocc,nindx,
C     i                   nl,nn,nclass,nlnax,
C     o                   il1a,in1a,il2a,in2a,ilna,nlna)
C
Cc gives an index for the allowed product phi
C
Cc nocc(l,n) = 0,1 ==> unoccupied, occupied
Cc nlna      = number of allowed product phi,
Cc             from integer function nalwln
C
Cc in1(lbas) = n1 etc.
Cc ilna(l1,n1,l2,n2) = index of allowed product phi
Cc ilna(l1,n1,l2,n2) = 0 ==> not allowed
C
C      implicit real*8(a-h,o-z)
C      dimension nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
C     i          nindx(nl,nclass)
C      dimension il1a(nlnax,nclass),in1a(nlnax,nclass),
C     o          il2a(nlnax,nclass),in2a(nlnax,nclass),
C     o          ilna(nl*nn*nl*nn,nclass),nlna(nclass)
C
C      do      ic = 1,nclass
C        nlna(ic)   = nalwln (nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn)
C
C        call indxlna (nocc(1,ic),nunocc(1,ic),nindx(1,ic),
C     d              nl,nn,nlna(ic),
C     o              il1a(1,ic),in1a(1,ic),il2a(1,ic),in2a(1,ic),
C     o              ilna(1,ic) )
C      end do
C
C      return
C      end
c---------------------------------------------------------------
      subroutine incor  (ncwf,nindxc,iclass,
     d                   nl,nnc,nclass,natom,
     o                   icore,ncore,nctot)

c 92.03.18
c sorts out allowed core states and count the number of core states

c ncwf(l,n,cl) = 1 ==> allowed, 0 ==> not allowed
c nindxc(l,cl)  = no. core states/l,class
c nl,nnc = max. no. l,n

c icore(i,cl) = index for allowed core states
c ncore(cl)   = no. allowed core states
c nctot       = total no. allowed core states

      implicit real*8 (a-h,o-z)
      dimension ncwf(0:nl-1,nnc,nclass),nindxc(0:nl-1,nclass),
     i          iclass(natom)
      dimension icore(nl*nl*nnc,nclass),ncore(nclass)

      ncx        = nl*nl*nnc
      do      ic = 1,nclass
        i          = 0
        j          = 0
        do       l = 0,nl-1
          do       n = 1,nindxc(l,ic)
            do       m = -l,l
              j          = j + 1
              if (ncwf(l,n,ic) == 1) then
                i          = i + 1
Cstop2rx 2013.08.09 kino                if (i > ncx) stop 'incore: wrong ncx'
                if (i > ncx) call rx( 'incore: wrong ncx')
                icore(i,ic)= j
              endif
            end do
          end do
        end do
        ncore(ic)  = i
      end do

c total no. allowed core states
      nctot      = 0
      do       i = 1,natom
        ic         = iclass(i)
        nctot      = nctot + ncore(ic)
      end do

      return
      end
c---------------------------------------------------------------
C      subroutine indxbas(nocc,nunocc,nindx,nl,nn,
C     i                   nallow,
C     o                   il1,in1,im1,
C     o                   il2,in2,im2,indxb)
C
Cc gives an index for the allowed product basis
C
Cc nocc(l,n) = 0,1 ==> unoccupied, occupied
Cc nallow    = number of allowed product basis,
Cc             from integer function nallow,
Cc             must be equal to final lbas
C
Cc in1(lbas) = n1 etc.
Cc indxb(n1,lm1,n2,lm2) = index of allowed product basis
Cc indxb(n1,lm1,n2,lm2) = 0 ==> not allowed
C
C      implicit real*8(a-h,o-z)
C      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
C     i          nindx(0:nl-1)
C      dimension indxb(nn,nl*nl,nn,nl*nl),
C     o          in1(nallow),il1(nallow),im1(nallow),
C     o          in2(nallow),il2(nallow),im2(nallow)
C
Cc$$$c write to  file PRODUCT if it exists
Cc$$$      ifprod     = ifile('PRODUCT')
Cc$$$      if (ifprod > 0)
Cc$$$     .write (ifprod,*)'lbas,n1,l1,im1,n2,l2,im2'
C
Cc initialise indxb
Cc      call iinit   (indxb,nn*nn*nl**4)
C      indxb = 0
C      lbas       = 0
C      do 10   l1 = 0,nl-1
C        ll1        = l1*l1
C      do 10   n1 = 1,nindx(l1)
C      do 10   m1 = 1,2*l1+1
C      do 10   l2 = 0,nl-1
C        ll2        = l2*l2
C      do 10   n2 = 1,nindx(l2)
C      do 10   m2 = 1,2*l2+1
Cc     indxb(n1,ll1+m1,n2,ll2+m2) = 0
C        if (nocc(l1,n1) == 0)   goto 10
C        if (nunocc(l2,n2) == 0) goto 10
C        if (l1 /= l2 .or. n1 /= n2) then
C          if (indxb(n2,ll2+m2,n1,ll1+m1) /= 0) goto 10
C        endif
C
Cc temporary
C        if (l1 == l2 .and. n1 == n2) then
C          if (indxb(n2,ll2+m2,n1,ll1+m1) /= 0) goto 10
C        endif
C
C        lbas       = lbas + 1
C        indxb(n1,ll1+m1,n2,ll2+m2) = lbas
C        in1(lbas)  = n1
C        il1(lbas)  = l1
C        im1(lbas)  = m1 - l1 - 1
C        in2(lbas)  = n2
C        il2(lbas)  = l2
C        im2(lbas)  = m2 - l2 - 1
Cc$$$        if (ifprod > 0)
Cc$$$     .write (ifprod,*)lbas,n1,l1,im1(lbas),n2,l2,im2(lbas)
C   10 continue
CCstop2rx 2013.08.09 kino      if(lbas /= nallow)stop 'indxbas: lbas /= nallow'
C      if(lbas /= nallow)call rx( 'indxbas: lbas /= nallow')
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine indxk(n1,n2,n3,qbas,
C     i                 qbz,nqbz,
C     o                 ginv,indk)
C
Cc generates indexing of k-points inside 1st BZ
Cc given k inside BZ indk(ia,ib,ic) = k-point inside 1st BZ
Cc where ia=a*n1, ib=b*n2, ic=c*n3 and
Cc k = a*G1 + b*G2 + c*G3, G1,G2,G3 are the reciprocal base vectors
C
Cc n1,n2,n3 = number of divisions along G1,G2,G3
Cc qbz      = k-points cartesian coordinates in the 1st BZ
Cc nqbz     = number of k-points
Cc qbas     = reciprocal base vectors
C
Cc indk     = s. above
Cc ginv     = inverse of base reciprocal vector
Cc notes:   all vectors are in units of 2*pi/a
C
C      implicit real*8(a-h,o-z)
C      dimension qbz(3,nqbz),qbas(3,3)
C      dimension ginv(3,3),
C     o          indk(-n1:n1,-n2:n2,-n3:n3)
C      dimension work(3)
C
Cc write to unit ifkp
Cc     ifkp       = ifile('KPNT')
Cc     if(ifkp > 0)then
Cc     write(ifkp,*)' k-points in the 1st BZ '
Cc     write(ifkp,*)' nG1,nG2,nG3, qx,qy,qz, index '
Cc     endif
C
Cc initialise indk
Cc      call iinit(indk,(2*n1+1)*(2*n2+1)*(2*n3+1))
C      indk=0
Cc invert qbas
C      call dinv33(qbas,0,ginv,det)
C
C      do       k = 1,nqbz
C        call dmv(3,3,ginv,3,qbz(1,k),work)
C        ia         = idnint(n1*work(1))
C        ib         = idnint(n2*work(2))
C        ic         = idnint(n3*work(3))
C        indk(ia,ib,ic) = k
Cc     if (ifkp > 0)
Cc    .write(ifkp,6000)ia,ib,ic,qbz(1,k),qbz(2,k),qbz(3,k),
Cc    .indk(ia,ib,ic)
Cc6000 format(3i4,3f10.4,i5)
C      end do
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine indxlm(nl,
C     o ila,ima)
C
Cc index of l,m
Cc i = (l+1)**2 + m - l
C
C      implicit real*8(a-h,o-z)
C      dimension ila(nl*nl),ima(nl*nl)
C
C      i          = 0
C      do       l = 0,nl-1
C        do       m = -l,l
C          i          = i + 1
C          ila(i)     = l
C          ima(i)     = m
C        end do
C      end do
CCstop2rx 2013.08.09 kino      if(i /= nl*nl)stop 'indxlm: wrong no. of l,m'
C      if(i /= nl*nl)call rx( 'indxlm: wrong no. of l,m')
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine indxln(nindx,nl,
C     i                  nn,nofln,
C     o il,in,iln)
C
Cc index of l,n
Cc nofln = number of l,n  from integer function nofln
C
C      implicit real*8(a-h,o-z)
C      dimension nindx(0:nl-1)
C      dimension il(nofln),in(nofln),iln(0:nl-1,nn)
C
C      i          = 0
C      do       l = 0,nl-1
C        do       n = 1,nindx(l)
C          i          = i + 1
C          il(i)      = l
C          in(i)      = n
C          iln(l,n)   = i
C        end do
C      end do
CCstop2rx 2013.08.09 kino      if(i /= nofln)stop 'indxln: wrong no. of ln'
C      if(i /= nofln)call rx( 'indxln: wrong no. of ln')
C
C      return
C      end
c---------------------------------------------------------------
C      subroutine indxlna(nocc,nunocc,nindx,nl,nn,
C     i                   nlna,
C     o                   il1a,in1a,il2a,in2a,ilna)
C
Cc gives an index for the allowed product phi
C
Cc nocc(l,n) = 0,1 ==> unoccupied, occupied
Cc nlna      = number of allowed product phi,
Cc             from integer function nalwln
C
Cc in1(lbas) = n1 etc.
Cc ilna(l1,n1,l2,n2) = index of allowed product phi
Cc ilna(l1,n1,l2,n2) = 0 ==> not allowed
C
C      implicit real*8(a-h,o-z)
C      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
C     i          nindx(0:nl-1)
C      dimension il1a(nlna),in1a(nlna),
C     o          il2a(nlna),in2a(nlna),ilna(0:nl-1,nn,0:nl-1,nn)
C
Cc      call iinit   (ilna,nl*nl*nn*nn)
C      ilna=0
C      lbas       = 0
C      do 10   l1 = 0,nl-1
C      do 10   n1 = 1,nindx(l1)
C      do 10   l2 = 0,nl-1
C      do 10   n2 = 1,nindx(l2)
Cc     ilna(l1,n1,l2,n2) = 0
C        if(nocc(l1,n1) == 0)goto 10
C        if(nunocc(l2,n2) == 0)goto 10
C        if((l1 /= l2 .or. n1 /= n2) .and. ilna(l2,n2,l1,n1) /= 0)goto 10
C        lbas       = lbas + 1
C        il1a(lbas) = l1
C        in1a(lbas) = n1
C        il2a(lbas) = l2
C        in2a(lbas) = n2
C        ilna(l1,n1,l2,n2) = lbas
C   10 continue
CCstop2rx 2013.08.09 kino      if(lbas /= nlna)stop 'indxlna: lbas /= nlna'
C      if(lbas /= nlna)call rx( 'indxlna: lbas /= nlna')
C
C      return
C      end
c--------------------------------------------------------------------
C      subroutine indxlnm(nindx,nl,
C     i                   nn,noflnm,
C     o il,in,im,ilnm)
C
Cc indexing of LMTO basis functions for a given class,
Cc follows that in TB-LMTO program
Cc il,in,im = l,n,m
Cc ilnm(n,lm) = index of n,l,m
Cc lm = l*l + l + m + 1
C
C      implicit real*8(a-h,o-z)
C      dimension nindx(0:nl-1)
C      dimension ilnm(nn,nl*nl),in(noflnm),il(noflnm),im(noflnm)
C
C      ind       = 0
C      do 1    l = 0,nl-1
C        l2        = l*l
C      do 1    n = 1,nindx(l)
C      do 1    m = 1,2*l+1
C        ind       = ind + 1
C        lm        = l2 + m
C        il(ind)   = l
C        in(ind)   = n
C        im(ind)   = m - l - 1
C        ilnm(n,lm) = ind
C    1 continue
CCstop2rx 2013.08.09 kino      if (ind /= noflnm) stop 'indxlnm: wrong no. l,n,m'
C      if (ind /= noflnm) call rx( 'indxlnm: wrong no. l,n,m')
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine idxlnma(nindx,nl,nclass,
C     i                   nn,nlnmx,
C     o il,in,im,ilnm)
C
Cc 92.jan.07
Cc indexing of LMTO basis functions for all classes,
Cc follows that in TB-LMTO program
Cc il,in,im = l,n,m
Cc ilnm(n,lm) = index of n,l,m
Cc lm = l*l + l + m + 1
C
C      implicit real*8(a-h,o-z)
C      dimension nindx(0:nl-1,nclass)
C      dimension ilnm(nn,nl*nl,nclass),
C     o          in(nlnmx,nclass),
C     o          il(nlnmx,nclass),
C     o          im(nlnmx,nclass)
C
C      do     ic = 1,nclass
C        ind       = 0
C        do 1    l = 0,nl-1
C          l2        = l*l
C        do 1    n = 1,nindx(l,ic)
C        do 1    m = 1,2*l+1
C          ind       = ind + 1
C          lm        = l2 + m
C          il(ind,ic)= l
C          in(ind,ic)= n
C          im(ind,ic)= m - l - 1
C          ilnm(n,lm,ic) = ind
C    1   continue
C      end do
C
C      return
C      end
c-------------------------------------------------------------------
      subroutine idxlnmc(nindxv,nindxc,
     d                   nl,nn,nnv,nnc,nlnmx,nlnmxv,nlnmxc,nclass,
     o                   il,in,im,ilnm,
     o                   ilv,inv,imv,ilnmv,
     o                   ilc,inc,imc,ilnmc)

c 92.jan.07
c 92.03.17 include core states
c indexing of core states and LMTO basis functions for all classes,
c follows that in TB-LMTO program
c il,in,im = l,n,m
c ilnm(n,lm) = index of n,l,m
c lm = l*l + l + m + 1
c NOTE: the indexing starts with core first and then valence on top
c       of core (not the same as index generated from nindx)
Ci  nindxv :number of partial waves for this l and site, read from GWinput col nnvv
Ci         :Eventually rename to nphiv
Ci  nindxc :number of core levels for this l and site, read from GWinput col nnc
Ci         :Eventually rename to nphic
Ci   nn    :max number of radial product functions for a given l (for dimensioning)
Ci         :Eventually rename to: ndphi
Ci   nnv   :max number of valence states for any l (for dimensioning)
Ci   nnc   :max number of core states for any l (for dimensioning)
Ci  nlnmx  :dimensions in,il,im,ilnm
Ci  nlnmxv :dimensions inv,ilv,imv,ilnmv
Ci  nlnmxc :dimensions inc,ilc,imc,ilnmc
      implicit real*8(a-h,o-z)
      dimension nindxv(0:nl-1,nclass),nindxc(0:nl-1,nclass)
      dimension ilnm(nn,nl*nl,nclass),
     o          ilnmv(nnv,nl*nl,nclass),
     o          ilnmc(nnc,nl*nl,nclass),
     o          in(nlnmx,nclass),il(nlnmx,nclass),im(nlnmx,nclass),
     o         inv(nlnmxv,nclass),ilv(nlnmxv,nclass),imv(nlnmxv,nclass),
     o         inc(nlnmxc,nclass),ilc(nlnmxc,nclass),imc(nlnmxc,nclass)

      do     ic = 1,nclass
        ind       = 0

c core
        do      l = 0,nl-1
          l2        = l*l
          do      n = 1,nindxc(l,ic)
            do      m = 1,2*l+1
              ind       = ind + 1
Cstop2rx 2013.08.09 kino              if (ind > nlnmx) stop 'idxlnmc: ind > nlnmx'
              if (ind > nlnmx) call rx( 'idxlnmc: ind > nlnmx')
              lm        = l2 + m
              il(ind,ic)= l
              in(ind,ic)= n
              im(ind,ic)= m - l - 1
              ilnm(n,lm,ic) = ind
              ilc(ind,ic)= l
              inc(ind,ic)= n
              imc(ind,ic)= m - l - 1
              ilnmc(n,lm,ic)= ind
            end do
          end do
        end do

c valence
        indv      = 0
        do      l = 0,nl-1
          l2        = l*l
          ncore     = nindxc(l,ic)
          do      n = 1,nindxv(l,ic)
Cstop2rx 2013.08.09 kino            if (ncore+n > nn) stop 'idxlnmc: ncore+n > nn'
            if (ncore+n > nn) call rx( 'idxlnmc: ncore+n > nn')
            do      m = 1,2*l+1
              ind       = ind + 1
              indv      = indv + 1
Cstop2rx 2013.08.09 kino              if (ind > nlnmx) stop 'idxlnmc: ind > nlnmx'
              if (ind > nlnmx) call rx( 'idxlnmc: ind > nlnmx')
              lm        = l2 + m
              il(ind,ic)= l
              in(ind,ic)= ncore + n
              im(ind,ic)= m - l - 1
              ilnm(ncore+n,lm,ic) = ind
              ilv(indv,ic)= l
              inv(indv,ic)= n
              imv(indv,ic)= m - l - 1
              ilnmv(n,lm,ic) = indv
            end do
          end do
        end do

      end do

      return
      end
c-------------------------------------------------------------------
      subroutine indxmdm (mdim,nclass,
     i                    iclass,natom,
     o imdim )

c given an atom, imdim points to the starting address of the
c total optimal product basis
c   1 2 3 ...  m1, m1+1, ...  m2, m2+1, ...  m3, m3+1 ...
c   ----atom1----  ----atom2----  ----atom2----  ----atom4----
c imdim(1) = 1
c imdim(2) = m1 + 1  etc.

      implicit real*8(a-h,o-z)
      dimension mdim(nclass),iclass(natom)
      dimension imdim(natom)

      imdim(1)   = 1
      do       i = 2,natom
cmod 93.06.07
c     ic         = iclass(i)
c     imdim(i)   = imdim(i-1) + mdim(i-1)
        ic         = iclass(i-1)
        imdim(i)   = imdim(i-1) + mdim(ic)
c-----
      end do

      return
      end
c-------------------------------------------------------------------
C      integer function indxpp (l1,n1,l2,n2,
C     i                         iln,nl,nn,nln)
C
Cc index for phi(l1,n1) phi(l2,n2), ln1 >= ln2 (lower triangle)
Cc indxpp = (l1,n1 ; l2,n2)
C
C      implicit real*8(a-h,o-z)
C      dimension iln(0:nl-1,nn)
C
C      ln1        = iln(l1,n1)
C      ln2        = iln(l2,n2)
C      irow       = max0(ln1,ln2)
C      icol       = ln1 + ln2 - irow
C      indxpp     = lowtri(nln,irow,icol)
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine indxql (lmx,
C     i                   il,in,nln,
C     i                   nql,
C     o iql,iql1,iqn1,iql2,iqn2,iqlln )
C
Cc indexing for qls
C
C      implicit real*8(a-h,o-z)
C      dimension il(nln),in(nln)
C      dimension iql(nql),iql1(nql),iqn1(nql),iql2(nql),iqn2(nql),
C     o          iqlln(0:lmx,nln,nln)
C
C      do       l = 0,lmx
C        do     ln1 = 1,nln
C          do     ln2 = 1,nln
C            iqlln(l,ln1,ln2) = 0
C          end do
C        end do
C      end do
C
C      i          = 0
C      do       l = 0,lmx
C        do     ln1 = 1,nln
C          l1         = il(ln1)
C          n1         = in(ln1)
C          do     ln2 = ln1,nln
C            l2         = il(ln2)
C            n2         = in(ln2)
C            if(l > l1+l2) goto 1111
C            i          = i + 1
C            iql(i)     = l
C            iql1(i)    = l1
C            iqn1(i)    = n1
C            iql2(i)    = l2
C            iqn2(i)    = n2
C            iqlln(l,ln1,ln2) = i
C 1111       continue
C          end do
C        end do
C      end do
C
C      return
C      end
c--------------------------------------------------------------------
C      subroutine indxrk(n1,n2,n3,qp,nq,g,ng,qbas,
C     o                  ginv,indxr,indxk)
C      implicit real*8(a-h,o-z)
C
Cc generates indexing of symmetry operators and k-point
Cc inside the irreducible BZ (IBZ).
Cc given k inside BZ indxr(ia,ib,ic) = symmetry operator
Cc                   indxk(ia,ib,ic) = k-point inside IBZ
Cc where ia=a*n1, ib=b*n2, ic=c*n3 and
Cc k = a*G1 + b*G2 + c*G3, G1,G2,G3 are the reciprocal base vectors
Cc input:
Cc n1,n2,n3 = number of divisions along G1,G2,G3
Cc qp       = k-points cartesian coordinates in the IBZ
Cc nq       = number of irreducible k-points
Cc g        = symmetry operators
Cc ng       = number of symmetry operators
Cc qbas     = reciprocal base vectors
Cc output:
Cc indxr, indxk
Cc remarks:
Cc all vectors are in units of 2*pi/a
C
C      dimension qp(3,nq),g(3,3,ng),qbas(3,3),ginv(3,3),
C     .          indxr(-n1:n1,-n2:n2,-n3:n3),work(3*3),
C     .          indxk(-n1:n1,-n2:n2,-n3:n3)
C
Cc initialise indxr and indxk
Cc      call iinit(indxr,(2*n1+1)*(2*n2+1)*(2*n3+1))
Cc      call iinit(indxk,(2*n1+1)*(2*n2+1)*(2*n3+1))
C      indxr=0
C      indxk=0
Cc invert qbas
C      call dinv33(qbas,0,ginv,det)
C
Cc rotate each k-point and label it with the symmetry operator
C      do       i = 1,ng
C        do       k = 1,nq
C          call dmv(3,3,ginv,3,qp(1,k),work)
C          ia         = idnint(n1*work(1))
C          ib         = idnint(n2*work(2))
C          ic         = idnint(n3*work(3))
C          indxr(ia,ib,ic) = i
C          indxk(ia,ib,ic) = k
C        end do
C      end do
C
C      return
C      end
c--------------------------------------------------------------------
      integer function maxocc (ifev,ispin,ef,
     d                         nlmto,nqbz)
c kino removed common w
c 92.02.28
c finds the max. no. occupied states from eigenvalues at all k
c NOTE: assumes that the eigenvalues are stored in direct access file

c ifev(2) = direct access unit file for eigenvalues
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level
c nlmto   = no. states
c nqbz    = no. k-points

      implicit real*8 (a-h,o-z)
      integer   ifev(2)

c integer work array
c      integer w
c      common /w/ w(1)
      real*8,allocatable :: w_iekt(:)

c maximum no. occupied states
c      call defdr   (iekt,nqbz*nlmto)
      allocate( w_iekt( nqbz*nlmto ) )
      noccx      = 0
      do      is = 1,ispin
c> read eigenvalues
        call rwdd    (ifev(is),
     d              nlmto,nqbz,
     o              w_iekt )
        noccxt     = noccx1 (w_iekt,nqbz,nlmto,ef)
        if (noccxt > noccx) noccx = noccxt
      end do
      maxocc     = noccx

c      call rlse    (iekt)
      deallocate( w_iekt )
      return
      end
c------------------------------------------------------------------
      integer function nallow (nocc,nunocc,nindx,nl,nn)

c gives the number of allowed product basis
c nocc(n,l) = 0,1 ==> unoccupied, occupied
c nallow    = number of allowed product basis

      implicit real*8(a-h,o-z)
      parameter (lmax=6,nnx=10)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension icheck(0:lmax,nnx,0:lmax,nnx)

Cstop2rx 2013.08.09 kino      if(nl-1 > lmax) stop 'nallow: increase lmax'
      if(nl-1 > lmax) call rx( 'nallow: increase lmax')
Cstop2rx 2013.08.09 kino      if(nn > nnx) stop 'nallow: increase nnx'
      if(nn > nnx) call rx( 'nallow: increase nnx')
c      call iinit   (icheck,nnx*nnx*(lmax+1)*(lmax+1))
      icheck=0
      do      l1 = 0,nl-1
        do      n1 = 1,nindx(l1)
          do      l2 = 0,nl-1
            do      n2 = 1,nindx(l2)
              icheck(l1,n1,l2,n2) = nocc(l1,n1)*nunocc(l2,n2)
              if (l1 /= l2 .or. n1 /= n2) then
                if (icheck(l1,n1,l2,n2)*icheck(l2,n2,l1,n1) /= 0)
     . icheck(l1,n1,l2,n2) = 0
              endif
            end do
          end do
        end do
      end do

      nallow     = 0
      do 10   l1 = 0,nl-1
      do 10   n1 = 1,nindx(l1)
      do 10   m1 = 1,2*l1+1

      do 10   l2 = 0,nl-1
      do 10   n2 = 1,nindx(l2)
      do 10   m2 = 1,2*l2+1

c     if (nocc(l1,n1) == 0)goto 10
c     if (nunocc(l2,n2) == 0)goto 10
        if (icheck(l1,n1,l2,n2) == 0) goto 10

c temporary
        if (l1 == l2 .and. n1 == n2 .and. m1 < m2)goto 10

        nallow     = nallow + 1
   10 continue

      return
      end
c------------------------------------------------------------------
      integer function nalwln (nocc,nunocc,nindx,nl,nn)

c gives the number of allowed product radial phi
c nocc(l,n)   = 0,1 ==> unoccupied, occupied
c nunocc(l,n) = 1,0 ==> unoccupied,occupied
c nalwln    = number of allowed phi(l1,n1) phi(l2,n2)

      implicit real*8(a-h,o-z)
      parameter (lmax=6,nnx=10)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension icheck(0:lmax,nnx,0:lmax,nnx)

Cstop2rx 2013.08.09 kino      if (nl-1 > lmax) stop 'nalwln: increase lmax'
      if (nl-1 > lmax) call rx( 'nalwln: increase lmax')
Cstop2rx 2013.08.09 kino      if (nn > nnx) stop 'nalwln: increase nnx'
      if (nn > nnx) call rx( 'nalwln: increase nnx')
c      call iinit   (icheck,nnx*nnx*(lmax+1)*(lmax+1))
      icheck=0
      nalwln     = 0
      do 10   l1 = 0,nl-1
      do 10   n1 = 1,nindx(l1)
        if(nocc(l1,n1) == 0)goto 10
        do 20   l2 = 0,nl-1
        do 20   n2 = 1,nindx(l2)
          if(nunocc(l2,n2) == 0)goto 20
          if((l1 /= l2 .or. n1 /= n2) .and. icheck(l2,n2,l1,n1) /= 0)
     . goto 20
          nalwln     = nalwln + 1
          icheck(l1,n1,l2,n2) = nalwln
   20   continue
   10 continue

      return
      end
c-------------------------------------------------------------------
C      subroutine ncorex (nindxc,
C     d                   nl,nclass,
C     o                   ncore,ncx)
C
Cc 92.03.18
Cc counts the no. core states/class
C
Cc nindxc(l,cl) = no. core states
C
Cc ncore(cl) = no. core states
Cc ncx       = max. no. core states
C
C      implicit real*8 (a-h,o-z)
C      dimension nindxc(0:nl-1,nclass)
C      dimension ncore(nclass)
C
C      do      ic = 1,nclass
C        ncore(ic)  = ivsum (nindxc(0,ic),nl)
C      end do
C      ncx        = maxim (ncore,nclass)
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine noccmax(ekt,nk,nt,ef,
C     o                   noccx)
C
Cc find the highest number of occupied states
C
Cc ekt   = eigenvalues for all k-points and states
Cc nk    = number of k-points in the 1st BZ
Cc nt    = number of states
Cc ef    = fermi level
C
Cc noccx = highest number of occupied states
C
C      implicit real*8(a-h,o-z)
C      dimension ekt(nk,nt)
C
Cc find the highest number of occupied states
C      noccx      = 0
C      do       k = 1,nk
C        do      it = 1,nt
C          if(ekt(k,it) > ef)goto 1111
C        end do
C 1111   if(it > noccx) noccx = it
C      end do
C
C      return
C      end
c------------------------------------------------------------------
      integer function nofblch(mdim,iclass,nclass,natom)

c number of bloch optimal product basis functions

      implicit real*8(a-h,o-z)
      dimension mdim(nclass),iclass(natom)

      nofblch    = 0
      do      ia = 1,natom
        ic         = iclass(ia)
        nofblch    = nofblch + mdim(ic)
      end do

      return
      end
c-------------------------------------------------------------------
      integer function noflmto(nindx,iclass,nl,nclass,natom)

c total number of LMTO basis functions

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass),iclass(natom)

      noflmto   = 0
      do 1    i = 1,natom
        ic        = iclass(i)
      do 1    l = 0,nl-1
        noflmto   = noflmto + (2*l+1)*nindx(l,ic)
    1 continue

      return
      end
c--------------------------------------------------------------------
      integer function nofln(nindx,nl)

c count the number of l,n

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1)

      nofln      = 0
      do       l = 0,nl-1
        nofln      = nofln + nindx(l)
      end do

      return
      end
c------------------------------------------------------------------
      integer function noflnm(nindx,nl)

c number of l,n,m

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1)

      noflnm    = 0
      do 1    l = 0,nl-1
        noflnm    = noflnm + nindx(l)*(2*l+1)
    1 continue

      return
      end
c-------------------------------------------------------------------
C      integer function noflnmx(nindx,nl,nclass)
C
Cc maximum number of l,n,m
C
C      implicit real*8(a-h,o-z)
C      dimension nindx(0:nl-1,nclass)
C
C      noflnmx   = 0
C      do     ic = 1,nclass
C        nlnm      = noflnm(nindx(0,ic),nl)
C        if(nlnm > noflnmx) noflnmx = nlnm
C      end do
C
C      return
C      end
c-------------------------------------------------------------------
C      subroutine nofpppp(nofln,
C     o npppp,npppd,nppdd)
C
Cc number of
Cc <phi(l1,n1) phi(l2,n2) phi(l3,n3) phi(l4,n4)>
Cc (l1,n1) =< (l2,n2) =< (l3,n3) =< (l4,n4)
C
Cc <phi(l1,n1) phi(l2,n2) phi(l3,n3) phidot(l4,n4)>
Cc (l1,n1) =< (l2,n2) =< (l3,n3)
C
Cc <phi(l1,n1) phi(l2,n2) phidot(l3,n3) phidot(l4,n4)>
Cc (l1,n1) =< (l2,n2)  ;  (l3,n3) =< (l4,n4)
C
C      implicit real*8(a-h,o-z)
C
C      npppp      = 0
C      do 10   i1 = 1, nofln
C      do 10   i2 = i1,nofln
C      do 10   i3 = i2,nofln
C      do 10   i4 = i3,nofln
C        npppp      = npppp + 1
C   10 continue
C
C      npppd      = 0
C      do 20   i1 = 1, nofln
C      do 20   i2 = i1,nofln
C      do 20   i3 = i2,nofln
C      do 20   i4 = 1 ,nofln
C        npppd      = npppd + 1
C   20 continue
C
C      nppdd      = 0
C      do 30   i1 = 1, nofln
C      do 30   i2 = i1,nofln
C      do 30   i3 = 1 ,nofln
C      do 30   i4 = i3,nofln
C        nppdd      = nppdd + 1
C   30 continue
C
C      return
C      end
c-----------------------------------------------------------------
C      integer function nofyyyy(nl)
C
Cc number of ylm(i1)*ylm(i2)*ylm(i3)*ylm(i4)
Cc i1 <= i2 <= i3 <= i4
C
C      implicit real*8(a-h,o-z)
C
C      nofyyyy    = 0
C      ll         = nl*nl
C      do 10   i1 = 1,ll
C      do 10   i2 = i1,ll
C      do 10   i3 = i2,ll
C      do 10   i4 = i3,ll
C        nofyyyy    = nofyyyy + 1
C   10 continue
C
C      return
C      end

C      subroutine nolna (nocc,nunocc,nindx,
C     d                  nl,nn,nclass,
C     o                  nlna)
C
Cc gives the number of allowed product radial phi
Cc nocc(l,n)   = 0,1 ==> unoccupied, occupied
Cc nunocc(l,n) = 1,0 ==> unoccupied,occupied
Cc nalwln    = number of allowed phi(l1,n1) phi(l2,n2)
C
C      implicit real*8(a-h,o-z)
C      dimension nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
C     i          nindx(nl,nclass)
C      dimension nlna(nclass)
C
C      do     ic = 1,nclass
C        nlna(ic)  = nalwln (nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn)
C      end do
C
C      return
C      end
