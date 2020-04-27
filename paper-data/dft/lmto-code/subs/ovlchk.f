      subroutine ovlchk(nbas,nbasp,pos,alat,rmax,hcr,clabl,
     .  ips,mode,plat,fovl,volsph)
C- Check volume and sphere overlaps
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   pos   :basis vectors
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   rmax  :augmentation radius, in a.u., by species (class)
Ci   hcr   :augmentation radius, in a.u., by species (class)
Ci   clabl :species (class) name
Ci   ips   :species (class) table: site ib belongs to species ips(ib)
Ci   mode  :same as mode in pairc.f
Ci   rmax,hcr (see remarks)
Co Outputs
Co   fovl, a function of the positive overlaps, for now set to
Co         sum (ovl)^6
Co   volsph sum of sphere volumes
Cr Remarks
Cr   Checks overlaps for both rmax and hcr.
Cr   hcr(1) <= 0 switches off part with hcr.
Cu Updates
Cu   24 Jul 15 Replace dclabl with clabl
Cu   14 Sep 11 Update printout when hcr is present
Cu   21 Aug 02 Can print out positions as multiples of plat
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nbas,nbasp
      double precision plat(3,3),pos(3,nbasp),rmax(*),hcr(*),alat,fovl,volsph
      character*8 clabl(*)
      integer ips(nbasp),mode(3)
C Local parameters
      double precision dd(3),dd1,dd2,sumrs,summt,ovlprs,ovlpmt,
     .  ctrs,ctmt,ddot,dlat(3),xx,avwsr,vol
      double precision qlat(3,3),volspp
      integer ibas,jbas,ic,jc,kc,m,ipr,i1mach,m1,m2,m3,isw,lgunit,stdo
      character*80 a, ch*1
      logical lterse,cmdopt,lhcr
      character*8 clabli,clablj

      call getpr(ipr)
      lhcr = hcr(1) > 0
      stdo = lgunit(1)
      call mkqlat(plat,qlat,xx)

C --- Determine which linear combination of plat is shortest ---
      dd1 = 9d9
      call dpzero(dlat,3)
      do  m1 = 1, 3
      do  m2 = 1, 3
      do  m3 = 1, 3
        if (mode(m1)*mode(m2)*mode(m3) /= 0) then
          do  ic = -1, 1
          do  jc = -1, 1
          do  kc = -1, 1
            call dpcopy(plat(1,m1),dd,1,3,dble(ic))
            call dpadd(dd,plat(1,m2),1,3,dble(jc))
            call dpadd(dd,plat(1,m3),1,3,dble(kc))
            dd2 = ddot(3,dd,1,dd,1)
            if (dd1 > dd2 .and. dd2 > 1d-6) then
              call dpcopy(dd,dlat,1,3,1d0)
              dd1 = dd2
            endif
          enddo
          enddo
          enddo
        endif
      enddo
      enddo
      enddo

      call info2(10,1,0,'   Site    Class%12fRmax%?;n;%8fHcr ;;'//
     .  '%16fPosition',isw(lhcr),0)
      volsph = 0d0
      volspp = 0d0
      do  ibas = 1, nbasp
        ic = ips(ibas)
        if (ipr <= 10) cycle
        if (ibas == nbas+1) write(stdo,'(''  ... Padding basis'')')
        clabli = clabl(ic)
        if (clabl(1) == ' ') call awrit1('%x%,4i',clabli,8,0,ic)
        if (lhcr) then
          write (stdo,1) ibas,ic,clabli,rmax(ic),hcr(ic),(pos(m,ibas),m=1,3)
    1     format(i5,3x,i4,2x,a8,2F12.6,3F11.5)
        else
          write (stdo,2) ibas,ic,clabli,rmax(ic),(pos(m,ibas),m=1,3)
    2     format(i5,3x,i4,2x,a8,f12.6,3F11.6)
          if (ipr >= 41) then
          call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,ibas),3,0d0,dd,3)
C          do  m = 1, 3
C            dd(m) = pos(1,ibas)*qlat(1,m)+pos(2,ibas)*qlat(2,m)+
C     .        pos(3,ibas)*qlat(3,m)
C          enddo
          write (stdo,3) dd
    3     format(13x,'as multiples of plat:',3F11.6)
          endif
        endif
        if (ibas <= nbas) volsph = volsph + 4.188790205d0*rmax(ic)**3
        volspp = volspp + 4.188790205d0*rmax(ic)**3
      enddo

      xx = avwsr(plat,alat,vol,nbas)
      if (ipr >= 10 .and. volspp == volsph) then
        call awrit3('%N Cell volume= %1,5;5d   Sum of sphere volumes='//
     .  ' %1,5;5d (%1;5d)',' ',80,i1mach(2),vol,volsph,volsph/vol)
      elseif (ipr >= 10) then
        volspp = 2*volspp-volsph
        call info5(0,1,0,
     .    ' Cell volume= %1,5;5d'//
     .    ' Sum of sphere volumes= %1,5;5d + %1,5;5d(2 x pad) '//
     .  ' ratio=%1;5d',vol,volsph,volspp-volsph,volspp/vol,0)
      endif

C --- Check sphere overlaps ---
      fovl = 0
      if (cmdopt('--terse=0',9,0,a)) return
      lterse = cmdopt('-terse',6,0,a) .or. cmdopt('--terse',7,0,a)
      if (lhcr .and. ipr > 20) write (stdo,4)
    4 format(/' ib jb',2x,'cl1     cl2',8x,'Pos(jb)-Pos(ib)',6x,
     .        'Dist  sumrs   Ovlp    %    summt   Ovlp   %')
      if ( .not. lhcr .and. ipr > 20) write (stdo,5)
    5 format(/' ib jb',2x,'cl1     cl2',8x,'Pos(jb)-Pos(ib)',6x,
     .        'Dist  sumrs   Ovlp    %')

      do  ibas = 1, nbasp
        ic = ips(ibas)
        if (ipr > 20) then
          clabli = clabl(ic)
          if (clabl(1) == ' ') call awrit1('%x%,4i',clabl,8,0,ic)
        endif
        do  jbas = ibas, nbasp
        jc = ips(jbas)
        if (rmax(ic) == 0 .and. rmax(jc) == 0) cycle
        if (ipr > 20) then
          clablj = clabl(jc)
          if (clabl(1) == ' ') call awrit1('%x%,4i',clablj,8,0,jc)
        endif
        if (ibas == jbas) then
          if (ddot(3,dlat,1,dlat,1) == 0) cycle
          call dcopy(3,dlat,1,dd,1)
        else
          forall (m = 1:3) dd(m) = pos(m,jbas) - pos(m,ibas)
          dd1 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
          call shorps(1,plat,mode,dd,dd)
          dd2 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
          if (dd2 > dd1+1d-9) call rx('bug in ovlchk')
C ...     test against shorbz
C         call mkqlat(plat,qlat,xx)
C         call shorbz(dd,dd,plat,qlat)
C         print *, dd2,dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
        endif
        forall (m = 1:3) dd(m) = dd(m)*alat
        dd1 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
        sumrs = rmax(ic) + rmax(jc)
        ovlprs = sumrs - dd1
        if (dd1 < 1d-6) then
          write (stdo,8) ibas,jbas,clabli,clablj,dd,dd1
          call rx('ovlchk: positions coincide')
        endif
        ctrs = nint(1000*ovlprs/dd1)/10d0
        if (lhcr) then
          summt = hcr(ic) + hcr(jc)
          ovlpmt = summt - dd1
          ctmt = nint(1000*ovlpmt/dd1)/10d0
        endif
        fovl = fovl + max(ovlprs/dd1,0d0)**6
        if ((lterse .or. ipr <= 40) .and. (ctrs <= -10 .or. ipr <= 20)) cycle
        ch = ' '
        if (ovlprs >= 0d0) ch='*'
        if (lhcr) then

          if (ctmt > -100) then
C           write(stdo,8) ibas,jbas,clabli,clablj,dd,dd1,sumrs,ovlprs,ctrs,ch,summt,ovlpmt,ctmt
            call info8(1,0,0,'%2,3i  '//clabli//clablj//'%3;7,3D%2;7,3D%;7,2D%;6,1D'//ch//'%;7,3D%;7,2D%;6,1D',
     .        [ibas,jbas],dd,[dd1,sumrs],ovlprs,ctrs,summt,ovlpmt,ctmt)
          else
          write(stdo,6) ibas,jbas,clabli,clablj,dd,dd1,
     .      sumrs,ovlprs,ctrs,ch,summt,ovlpmt
          endif
    6     format(2i3,2x,a8,a8,3f7.3,2f7.3,f7.2,f6.1,a1,f7.3,f7.2,'  --')
        else
C         write(stdo,7) ibas,jbas,clabli,clablj,dd,dd1,sumrs,ovlprs,ctrs,ch
          call info8(1,0,0,'%2,3i  '//clabli//clablj//'%3;7,3D%2;7,3D%;7,2D%;6,1D'//ch,
     .      [ibas,jbas],dd,[dd1,sumrs],ovlprs,ctrs,summt,ovlpmt,ctmt)
    7     format(2i3,2x,a8,a8,3f7.3,2f7.3,f7.2,f6.1,a1)

        endif
        enddo
      enddo
    8 format(2i3,2x,a8,a8,3f7.3,2f7.3,f7.2,f6.1,a1,f7.3,f7.2,f6.1)
      end
