      subroutine bangs(dist0,nbas,bas,alat,dclabl,ipc)
C- Print table of bond angles and interatomic distances
C- Adapted from pr_bndang.f in nfp (see also mol/bonds.f)
      implicit none
      integer nnmax
      parameter (nnmax=100)
      integer nbas,ipc(*)
      double precision bas(3,nbas),alat,dclabl(*)

      integer nn,nd,na,numa,ipr,iprint,ic,jc,
     .        i,j,ib,jb,i1,i2,j1,j2,k1,k2,in1,in2,jn1,jn2,kn1,kn2,ip
      integer i1mach
      character*8 spid1,spid2,spnum1,spnum2
      character*8 dstr(2,200),astr(2,200),dnstr(2,200),anstr(2,200),
     .            spp(50),sppn(50)
      character*40 s1,s2
      double precision pos(3,nnmax),dd(3),dda,aj,ai,dist(200),
     .                 dist0,sp,angs(200),angle

      ipr=iprint()

c --- start loop over atoms; set up list of neighbors
      call awrit1('%NBond angles and bond lengths up to %d a.u.',' ',
     .            120,i1mach(2),dist0)
      write(6,788)
      do 5 ib = 1,nbas
        ic = ipc(ib)
        nn = 0
        nd = 0
        call r8tos8(dclabl(ic),spid1)
        ip = 0
        spnum1 = ' '
        call bin2a(' ',0,0,ib,2,0,16,spnum1,ip)
        do 1 jb = 1,nbas
          jc = ipc(jb)
          call r8tos8(dclabl(jc),spid2)
          ip = 0
          spnum2 = ' '
          call bin2a(' ',0,0,jb,2,0,16,spnum2,ip)
          dd(1) = alat*(bas(1,jb)-bas(1,ib))
          dd(2) = alat*(bas(2,jb)-bas(2,ib))
          dd(3) = alat*(bas(3,jb)-bas(3,ib))
          dda = dsqrt(dd(1)**2+dd(2)**2+dd(3)**2)
          if(dda <= dist0.and.dda > 0.1d0) then
            nn = nn+1
            if (nn > nnmax) call rx('bonds: increase nnmax')
            pos(1,nn) = dd(1)
            pos(2,nn) = dd(2)
            pos(3,nn) = dd(3)
            spp(nn) = spid2
            sppn(nn) = spnum2
            call anput1(nd,dist,dstr,dnstr,dda,spid1,spid2,
     .                 spnum1,spnum2)
          endif
    1   continue
        call ansrt1(nd,dist,dstr,dnstr)

c --- get the bond angles
        na = 0
        numa = 0
        do 3 i = 1,nn
          do 2 j = 1,i-1
            ai = dsqrt(pos(1,i)**2+pos(2,i)**2+pos(3,i)**2)
            aj = dsqrt(pos(1,j)**2+pos(2,j)**2+pos(3,j)**2)
            sp = pos(1,i)*pos(1,j)+pos(2,i)*pos(2,j)+pos(3,i)*pos(3,j)
            angle = dacos(min(max(sp/(ai*aj),-1d0),1d0))
     .              *360d0/6.2831853d0
            if(angle <= 180.00001d0) then
              numa = numa+1
              call anput1(na,angs,astr,anstr,angle,spp(i),spp(j),
     .                   sppn(i),sppn(j))
            endif
    2     continue
    3   continue
        call ansrt1(na,angs,astr,anstr)

c --- printout
        do 4 i = 1,max0(nd,na)
          call strip(spid1,i1,i2)
          call strip(spnum1,in1,in2)
          s1 = ' '
          if (i <= nd) then
            call strip(dstr(2,i),j1,j2)
            call strip(dnstr(2,i),jn1,jn2)
            write(s1,471) dist(i),
     .         spid1(i1:i2),dstr(2,i)(j1:j2),
     .         spnum1(in1:in2),dnstr(2,i)(jn1:jn2)
  471       format(f12.9,3x,a,'-',a,1x,a,'-',a)
          endif
          s2 = ' '
          if (i <= na) then
            call strip(astr(1,i),j1,j2)
            call strip(astr(2,i),k1,k2)
            call strip(anstr(1,i),jn1,jn2)
            call strip(anstr(2,i),kn1,kn2)
            write(s2,472) angs(i),
     .         astr(1,i)(j1:j2),spid1(i1:i2),astr(2,i)(k1:k2),
     .         anstr(1,i)(jn1:jn2),spnum1(in1:in2),anstr(2,i)(kn1:kn2)
  472       format(f16.8,2x,a,'-',a,'-',a,1x,a,'-',a,'-',a)
          endif

          if (i == 1) write(6,786) ib,spid1,s1,s2
          if (i > 1) write(6,787) s1,s2
  786     format(/i4,2x,a8,a30,a35)
  787     format( 4x,2x,8x,a30,a35)
  788     format(/' atom    spec',9x,'bond length',16x,' bond angle')

    4   continue

    5 continue

      end

c --- sub anput1: add entry to list
      subroutine anput1(na,angs,astr,anstr,angle,sp1,sp2,spn1,spn2)
!       implicit real*8 (a-h,p-z)
      implicit none
      real(8) :: angle
      integer na
      double precision angs(*)
      character*8 astr(2,1),anstr(2,1),a1,a2,n1,n2
      character*(*) sp1,sp2,spn1,spn2
      a1 = sp1
      a2 = sp2
      n1 = spn1
      n2 = spn2
      na = na+1
      angs(na) = angle
      astr(1,na) = a1
      astr(2,na) = a2
      anstr(1,na) = n1
      anstr(2,na) = n2
      end

c --- sub ansrt1: sort by bond length
      subroutine ansrt1(na,angs,astr,anstr)
!       implicit real*8 (a-h,p-z)
      implicit none
      integer :: na, i0, j, i
      real(8) :: bot, angs(1), ang1(200)
      character*8 astr(2,1),astr1(2,200),anstr(2,1),anstr1(2,200)

      do i = 1,na
        i0 = 0
        bot = 10000
        do j = 1,na
          if(angs(j) < bot) then
            i0 = j
            bot = angs(j)
          endif
        enddo
        ang1(i) = angs(i0)
        astr1(1,i) = astr(1,i0)
        astr1(2,i) = astr(2,i0)
        anstr1(1,i) = anstr(1,i0)
        anstr1(2,i) = anstr(2,i0)
        angs(i0) = 10000
      enddo

      do i = 1,na
        angs(i) = ang1(i)
        astr(1,i) = astr1(1,i)
        astr(2,i) = astr1(2,i)
        anstr(1,i) = anstr1(1,i)
        anstr(2,i) = anstr1(2,i)
      enddo

      end

