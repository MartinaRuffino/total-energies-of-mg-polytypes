C#define VEC
      subroutine maks0(npr,plat,bas,indxcg,jcg,cg,cy,alphv,adotv,iax,e,
     .  iop,lmaxw,s0,sd,ndimW)
C- Real space unscreened structure constants centered at one site
C ----------------------------------------------------------------
Ci Inputs
Ci   npr   :number of pairs in this cluster
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   bas   :basis vectors (input)
Ci   cy    :normalization constants for spherical harmonics (sylmnc.f)
Ci   cg,indxcg,jcg:Clebsch Gordan coefficients, and indices (scg.f)
Ci   alphv :vector of alpha's to be added to s0
Ci   adotv :vector of adots's to be added to sd
Ci   iax   :array of parameters containing info about each pair
Ci   e     :energy
Ci   iop   :1s digit
Ci           0 make s0 and sdot
Ci           1 make s0 only
Ci           2 make sd only
Ci           else, make both s0 and sd
Ci          10s digit
Ci           1 add alphv to s0
Ci           2 add adotv to sd
Ci           else both 1 and 2
Ci          100s digit
Ci           1 scale s0 by -1
Ci           2 scale sd by -1
Ci           else both 1 and 2
Ci          1000s digit
Ci           0 Methfessel conventions
Ci           1 Scale strux to conform to to Andersen conventions:
Ci             scale by factor : (2l-1)!! (2l''-1)!! / 2
Ci   lmaxw :maximum l for Watson-sphere (input)
Ci   ndimW :leading dimension of s0,sd
Co Outputs
Co   s0    :alphv^-1 +/- s0, depending on iop
Co   sd    :adotv^-1 +/- s0-dot, depending on iop
Cr Remarks
Cr   Conventions for matrix s0.
Cr   ROWS of matrix s0 are the 'field' or 'augmentation' dimension
Cr   COLS of matrix s0 are the 'source' or 'basis' dimension
Cr   Quantity dr is vector R_source - R_expansion.
Cr   Example: for two atoms at (0,0,0) and (0,0,z) S would be
Cr           AUG  \  SOURCE, iax(1) ---->
Cr          iax(2)       (S(dr=0)    S(dr=z))
Cr            |     S =  (                  )
Cr            |          (S(dr=-z)   S(dr=0))
Cu Updates
Cu   06 Aug 06 Added 1000s digit (1000 was assumed before)
Cu    8 Sep 00 Bug fix handling combination of lmaxw>=0 and sdot
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer npr,iop,lmaxw,ndimW,niax
      parameter (niax=10)
      integer indxcg(*),jcg(*),iax(niax,*)
      double precision cg(*),cy(*),e,plat(3,3),alphv(*),
     .  adotv(*),bas(3,*),s0(ndimW,ndimW),sd(ndimW,ndimW)
C ... Dynamically allocated arrays
      integer, allocatable :: offs(:),lmxb(:)
      real(8), allocatable :: x(:),y(:),z(:),s(:),d(:),hl(:),hd(:)
      integer, pointer :: jj(:)
!      include 'structures.h'
      type(str_strxsu) s_strx
C ... Local parameters
      integer i,j,lm1,lm2,ll,iop0,iop1,loka,nlmi,nlmj,offi,offj,nlmw
      integer nlmbx,nlmp,npow,offh,iprint,k1,k2,nprw
      double precision dr(3),dsqr,drr2,fs,fd
      external mstrx2,mstrx3,drr2,iprint

C ... Setup
      if (npr < 0) return
      call tcn('maks0')
      fs = 1
      iop0 = mod(iop,10)
      if (iop0 > 2) iop0 = 0
      iop1 = mod(iop/10,10)
      i = mod(iop/100,10)
      if (mod(i,2) == 1) fs = -1
      fd = 1
      if (i > 1) fd = -1
      loka = mod(iop/1000,10)
      nprw = npr
      if (lmaxw >= 0) nprw = npr+1
      nlmw = (lmaxw+1)**2

      allocate(offs(0:npr+1),lmxb(0:npr))
      allocate(x(npr),y(npr),z(npr))

C ... Make offsets
      offs(1) = 0
      nlmbx = 0
      do  i = 1, npr
        offs(i+1) = offs(i) + iax(9,i)
        lmxb(i) = ll(iax(9,i))
        nlmbx = max(nlmbx,iax(9,i))
      enddo

      do  i = 1, nprw
C   ... Add diagonal to s0
        offi = offs(i)
        nlmi = iax(9,min(i,npr))
        if (i > npr) nlmi = nlmw
        if (iop1 /= 2) then
          do  lm2 = 1, nlmi
            do  lm1 = 1, nlmi
              s0(offi+lm1,offi+lm2) = 0d0
            enddo
          enddo
          do  lm1 = 1, nlmi
            s0(offi+lm1,offi+lm1) = alphv(offi+lm1)
          enddo
        endif
        if (iop1 /= 1) then
          do  lm2 = 1, nlmi
            do  lm1 = 1, nlmi
              sd(offi+lm1,offi+lm2) = 0d0
            enddo
          enddo
          do  lm1 = 1, nlmi
            sd(offi+lm1,offi+lm1) = adotv(offi+lm1)
          enddo
        endif
        if (i == nprw) cycle

C   --- Put strux into array s0 ---
        if (i == npr) goto 34
C   ... Setup for vectorized strux
C#ifdef VEC
        do  j = i+1, npr
C         dr = pos(j) - pos(i)
          dsqr = drr2(plat,bas(1,iax(2,i)),bas(1,iax(2,j)),
     .      iax(3,j)-iax(3,i),iax(4,j)-iax(4,i),iax(5,j)-iax(5,i),dr)
          k1 = j-i
          call dvset(x,k1,k1,dr(1))
          call dvset(y,k1,k1,dr(2))
          call dvset(z,k1,k1,dr(3))
        enddo
        call strxsu(nlmi,-1,lmxb(i+1),1,loka,npr-i,-1,cg,jcg,indxcg,
     .    nlmbx,nlmp,npow,s_strx)
        allocate(hl(nlmp*npr),hd(nlmp*npr))
        call rstr0(1,lmxb(i+1),e,nlmp,npr-i,x,y,z,
     .    lmxb(i),mod(iop0,2),hl,hd)
C#endif

C   ... For j>i, get strux for (i,j) pair
        if (iop0 == 1) then
          allocate(s(nlmi*ndimW))
          allocate(d(1))
        else
          allocate(s(nlmi*nlmbx))
          allocate(d(nlmi*nlmbx))
        endif
        k2 = 0
        do  j = i+1, npr
          k1 = j-i
          if (k1 <= k2) cycle
          k2 = k1
          offj = offs(j)
          nlmj = iax(9,j)
C#ifndefC VEC
C          dsqr = drr2(plat,bas(1,iax(2,i)),bas(1,iax(2,j)),
C     .         iax(3,j)-iax(3,i),iax(4,j)-iax(4,i),iax(5,j)-iax(5,i),dr)
C          call mstrx2(e,dr,nlmi,nlmj,nlmi,cg,indxcg,jcg,cy,10*loka+iop0,
C     .      s,d)
C          if (iop0 /= 2) call dmscop(s0,ndimW,s,nlmi,1,nlmi,1,
C     .      nlmj,offi+1,offj+1,fs)
C          if (iop0 /= 1) call dmscop(sd,ndimW,d,nlmi,1,nlmi,1,
C     .      nlmj,offi+1,offj+1,fd)
C#else
          jj => s_strx%jkl(lmxb(j))%p

C         For strux only, make entire row in one fell swoop
          if (iop0 == 1) then
   35       continue
            if (k2+i < npr) then
              if (lmxb(k2+i+1) == lmxb(j)) then
                k2 = k2+1
                goto 35
              endif
            endif
c           k2 = k1
            call hstrux(e,nlmi,nlmj,nlmp,npow,k1,k2,s_strx%ikl,jj,
     .        s_strx%ip,s_strx%cf,hl,s)
            call dmscop(s0,ndimW,s,nlmi,1,nlmi,1,nlmj*(k2-k1+1),
     .        offi+1,offj+1,fs)
C           call prmx('s',s0(offi+1,offj+1),ndimW,nlmi,nlmj*(k2-k1+1))
C         For strux+sdot, make row block by block
          else
            offh = nlmp*(j-i-1)
            call hstrud(e,nlmi,nlmj,nlmp,npow,offh,s_strx%ikl,jj,
     .        s_strx%ip,s_strx%cf,hl,hd,s,d)
            if (iop0 /= 2)
     .      call dmscop(s0,ndimW,s,nlmi,1,nlmi,1,nlmj,offi+1,offj+1,fs)
            call dmscop(sd,ndimW,d,nlmi,1,nlmi,1,nlmj,offi+1,offj+1,fd)
          endif
C#endif
        enddo
        deallocate(s_strx%cf,s_strx%ip,s_strx%ikl)
        do  j = 0, ll(nlmbx)
          if (associated(s_strx%jkl(j)%p)) then
            deallocate(s_strx%jkl(j)%p)
          endif
        enddo
        deallocate(hl,hd,s,d)
   34   continue

C   ... Make strux for (i,Watson-sphere) pair
        if (lmaxw >= 0) then
          allocate(s(nlmi*nlmw))
          allocate(d(nlmi*nlmw))
          offj = offs(npr+1)

          nlmj = nlmw
          dsqr = drr2(plat,bas(1,iax(2,i)),bas(1,iax(2,1)),
     .         iax(3,1)-iax(3,i),iax(4,1)-iax(4,i),iax(5,1)-iax(5,i),dr)
          call mstrx3(e,dr,nlmi,nlmj,nlmi,cg,indxcg,jcg,cy,iop0,
     .      s,d)
          if (iop0 /= 2) call dmscop(s0,ndimW,s,nlmi,1,nlmi,1,
     .      nlmj,offi+1,offj+1,fs)
          if (iop0 /= 1) call dmscop(sd,ndimW,d,nlmi,1,nlmi,1,
     .      nlmj,offi+1,offj+1,fd)
          deallocate(s,d)
        endif
      enddo

      deallocate(offs,lmxb,x,y,z)

      if (iop1 /= 2) then
        do  i = 1, ndimW
          do  j = 1, i-1
            s0(i,j) = s0(j,i)
          enddo
        enddo
      endif
      if (iop1 /= 1) then
        do  i = 1, ndimW
          do  j = 1, i-1
            sd(i,j) = sd(j,i)
          enddo
        enddo
      endif

C      if (iop1 /= 2) then
C       call prmx('maks0: s0',s0,ndimW,ndimW,ndimW)
C      endif
C      if (iop1 /= 1) then
C        call prmx('maks0: sd',sd,ndimW,ndimW,ndimW)
C      endif

      call tcx('maks0')
      end
