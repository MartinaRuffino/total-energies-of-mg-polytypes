
! read sig from 'se' file and projects it. Part extracted from sudmft.f but used at the moment. might be used latter so keep here.


!      if (cmdopt('--rsgw',6,0,strn))then
!         rsgwdiag=.true.        !off diag element not yet implemented
!         nblstf = NULLI
!         fn = 'se'
!         ifi = fopna(fn,-1,i)
!         lunits=1
!         lascii=20
!
!         i = 1000 + lunits*100+lascii+4 !read nblst,assumed unit eV,ascii,read only
!         call ioseh(i,fmode,j,nspse,nband,nqf,nomg,ommin,ommax,chempotgw,nblstf,ix,ix,ifi)
!         allocate(iblstf(nband))
!         do ib=1, nband
!            iblstf(ib) = ib
!         enddo
!         nblstf = nband
!         ioff = iblstf(1)-1
!!     if (chempot /= NULLI .and. lunits == 1) chempotgw = chempotgw/rytoev ! chempot always in Ry units
!         lqsgws = mod(fmode/10,2) >= 1 ! T if have static QSGW sigma
!         lhf = mod(fmode/10,4) >= 2 ! T if have Fock Exchange
!         llxc = mod(fmode/10,10) >= 4 ! T if have LDA exchange-correlation potential
!         allocate(qpf(3,nqf),eigf(nband,nqf,nspse),sigvf(nband,nqf,nspse),sigwf(nomg,nband,nqf,nspse))
!
!         call dpzero(sigvf,size(sigvf))
!         call iose2(lascii,fmode,fmode,nspse,nband,nqf,nomg,qpf,eigf,sigvf,sigwf,sexf,vxcf,ifi)
!         forall (i=1:nomg) sigwf(i,:,:,:)=sigwf(i,:,:,:)-sigvf(:,:,:) !Why ?!
!         multi=300
!         nomgx=(nomg-1)*multi + 1
!         allocate(omgn(nomg),omgx(nomgx),sigwfx(nomgx,nband,nqf,nspse),wki(nomgx),wkii(nomg),ag(nomgx))
!         forall (iomg= 1:nomg) omgn(iomg)=ommin+(ommax-ommin)*dble(iomg-1)/dble(nomg-1)
!         forall (iomg= 1:nomgx) omgx(iomg)=ommin+(ommax-ommin)*dble(iomg-1)/dble(nomgx-1)
!         call info0(10,1,0,' se file is read and interpolate on finer energy mesh')
!         do i=1,nband
!            do iqfbz = 1, nqf
!               do isp = 1,nsp
!                  eigfq=eigf(i,iqfbz,isp)
!                  forall (iomg=1:nomg) wkii(iomg)=sigwf(iomg,i,iqfbz,isp)
!                  call sfitrp(102,nomg,nomgx,ommin,ommax,multi,eigfq,wkii,
!     .                 0.001,omgn,omgx,wki,ag,ag,ag,ag,ib,ib,ib)
!
!                  do iomg= 1, nomgx
!                     sigwfx(iomg,i,iqfbz,isp)=wki(iomg)
!                  enddo
!
!               enddo
!            enddo
!         enddo
!
!         deallocate(sigwf,wki)
!         allocate(sigwf(nomgx,nband,nqf,nspse))
!         sigwf(:,:,:,:)=sigwfx(:,:,:,:)
!         nomg=nomgx
!
!
!
!
!
!c$$$         call mpibc1(nqf,1,2,.false.,'','')
!c$$$         call mpibc1(qpf,1,2,.false.,'','')
!c$$$         call mpibc1(sigwf,1,2,.false.,'','')
!c$$$         call mpibc1(eigf,1,2,.false.,'','')
!!     sigwf is the GW self energy (be careful with the sign!!!!!!)
!!     Now, need to generate the projector for the q point.
!
!         s_ctrl%lwsig = LW5     ! bndfp will write evals, evecs to file
!         s_ctrl%plbnd = 2       ! bndfp will not make output density
!         s_optic%loptic = 0     ! bndfp will not generate optics
!         s_ctrl%lfp = s_ctrl%lfp - bitand(s_ctrl%lfp,2) ! Do not shorten q vectors
!         nkp = s_bz%nkp
!         allocate(s_ham%evals(ndham,nsp,nkp))
!         allocate(s_pot%ppn(1))
!         lrout = 0; lfrce = 0
!         s_ctrl%lfp = s_ctrl%lfp - IAND(s_ctrl%lfp,4) + 4 ! Retain hab and sab
!         s_ctrl%lfp = s_ctrl%lfp - IAND(s_ctrl%lfp,8) + 8 ! Update s_ctrl%elind
!
!C     bndfp will need to replace nkp, s_bz%qp
!         call bndfp(s_strn,s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot,s_bz,
!     .        s_optic,s_gw,nbas,nsp,0,0,0,0,lrout,lfrce,0,xv,1,1,xv,xv,xv)
!         ppnl => s_pot%ppn
!         lrsig = s_ham%lsig
!         nkfbz = s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3)
!         call info0(10,1,0,' Make and renormalize projectors')
!         ifi = fopna('proj',-1,4); rewind ifi ! Ensure proj file is rewound
!         i = 0 ; if (s_dmft%knorm == 0) i = 20
!C     print 777, 'sudmft call lmlproj',procid
!         allocate(s_dmft%Olapp(ldcix,ldcix,nsp,ncix))
!         call lmlproj(11,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,sigdc)
!c$$$
!         allocate(dmftu(nlohi(1):nlohi(2),ldcix,nsp,ncix))
!         ifi = fopna('proj',-1,5); rewind ifi
!         nk1=s_bz%nkabc(1)
!         nk2=s_bz%nkabc(2)
!         nk3=s_bz%nkabc(3)
!
!         nkfbz =nk1*nk2*nk3
!         if(nkfbz .ne. s_gw%nkabc(1)*s_gw%nkabc(1)*s_gw%nkabc(1)) call rx('sudmft: k-point mismatch with file se') ! kmesh needs to be equal to gw one
!         call bzmsh00(s_lat%plat,s_bz%lshft,0,s_bz%nkabc,ifac,qb) ! Needed for iqstar
!
!         allocate(sigwloc(ldcix,nsp,ncix,nomg))
!         allocate(ggwloc(ldcix,nsp,ncix,nomg))
!         allocate(gl(ldcix,nsp,ncix,nomg))
!         call dpzero(sigwloc,2*size(sigwloc))
!         call dpzero(ggwloc,2*size(sigwloc))
!         call dpzero(gl,2*size(sigwloc))
!
!
!         do iqfbz = 1, nkfbz
!            call iodmftu(s_dmft,.true.,ndham,nlohi(1),nlohi(2),ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)
!            call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,iv)
!            do iomg=1, nomg
!               call makesiggwloc(1,s_dmft,ncix,ldcix,nsp,iq,nqf,iomg,nomg,ommin,ommax,nband,eigf,sigwf,dmftu,ggwloc)
!               call makesiggwloc(2,s_dmft,ncix,ldcix,nsp,iq,nqf,iomg,nomg,ommin,ommax,nband,eigf,sigwf,dmftu,sigwloc)
!               call makesiggwloc(11,s_dmft,ncix,ldcix,nsp,iq,nqf,iomg,nomg,ommin,ommax,nband,eigf,sigwf,dmftu,gl)
!            enddo               !iomg
!         enddo
!
!         allocate(siggwl(s_dmft%ndsig,nomg))
!         allocate(ggwl(s_dmft%ndsig,nomg))
!         call dpzero(siggwl,2*size(siggwl))
!         call dpzero(ggwl,2*size(siggwl))
!
!
!         allocate(nicixi(s_dmft%ncix)) ! Loop over inequivalent cix only
!
!         call ineqcix(s_dmft,s_dmft%ncix,nicixi)
!         do iomg=1, nomg
!            do isp=1,nsp
!               do i=1, s_dmft%ncix
!                  if (nicixi(i) >= 0) cycle ! skip equivalent cix
!                  do j = s_dmft%nzsigi(i-1)+1, s_dmft%nzsigi(i)
!                     if (s_dmft%iasig(1,j) == s_dmft%iasig(2,j)) then
!                        siggwl(j,iomg)=sigwloc(s_dmft%iasig(1,j),isp,i,iomg)
!                        ggwl(j,iomg)=ggwloc(s_dmft%iasig(1,j),isp,i,iomg)
!                     endif
!                  enddo
!               enddo
!            enddo
!         enddo
!
!         deallocate(nicixi)     ! Loop over inequivalent cix only
!
!         call info0(30,0,0,' Writing files G^gw_loc ...')
!         ifig = fopna('gloc_GW',-1,2); rewind ifi
!         ifis = fopna('sig_GW',-1,2); rewind ifi
!         ifid = fopna('dos',-1,2); rewind ifi
!         do iomg=1, nomg
!            write(ifig,'(2(x,f14.8))',advance='no')  ommin+(ommax-ommin)*dble(iomg-1)/dble(nomg-1)
!            write(ifis,'(2(x,f14.8))',advance='no')  ommin+(ommax-ommin)*dble(iomg-1)/dble(nomg-1)
!            write(ifid,'(2(x,f14.8))',advance='no')  ommin+(ommax-ommin)*dble(iomg-1)/dble(nomg-1)
!            do i=1, s_dmft%ndsig
!               write(ifis,'(2(x,f14.8))',advance='no') siggwl(i,iomg)/nkfbz
!               write(ifig,'(2(x,f14.8))',advance='no') ggwl(i,iomg)/nkfbz
!            enddo
!            write(ifid,'(2(x,f14.8))',advance='no') gl(1,1,1,iomg)/nkfbz
!            write(ifid,*)
!            write(ifis,*)
!            write(ifig,*)
!         enddo
!         call fclose(ifis)
!         call fclose(ifig)
!         call rx0('end of rsgw')
!      endif                     !rsgw



      subroutine makesiggwloc(mode,s_dmft,ncix,ldcix,nsp,iq,nqf,iomg,nomg,ommin,ommax,nband,eigf,sigw,dmftu,output)
C- Compute the local sigma from GW sum_kn U^+(nk,l)siggw(nk)U(nk,l' )
Cio Structures1
Cio  s_dmft
Ci     Elts read:  lsigim icix l
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lsigim
Cio    Passed to:  *
Ci Inputs
C     i mode : 1digit :
Ci                         1   output Sig_gw_loc
Ci                         2 output G_gw_loc
Ci           10s   digit : 0 U= usual projector
Ci                         1  U=1 ( for debuging or comparaison
Ci ncix   :dimensioning parameter, number of correlated l/site channels
Ci ldcix  :dimensioning parameter, max number of correlated orbitals (lms character)
Ci nsp    :dimensioning parameter,2 for spin-polarized case, otherwise 1
Ci ip     : index of irr q point
Ci nqf    : nb of q pts
Ci nomg,ommin,ommax : parameter of the energy mesh
Ci nband  : nb of band
Ci eigf   : GW e.v
Ci siggw  : GW self energy
Ci dmftu  : projector
ci
      !siggwloc_ll'=
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_dmft)         :: s_dmft
      integer, intent(in)    :: ldcix,nsp,ncix,iq,nqf,iomg,nomg,nband,mode
      real(8), intent(in)    :: ommin,ommax

      complex(8), intent(in)  :: sigw(nomg,nband,nqf,nsp)
      real(8), intent(in)  :: eigf(nband,nqf,nsp)
      complex(8), intent(in)  :: dmftu(s_dmft%nlohi(1):s_dmft%nlohi(2),ldcix,nsp,ncix)
      complex(8), intent(out) :: output(ldcix,nsp,s_dmft%ncix,nomg)
      complex(8) ::new,g
      real(8) ::w
      integer :: mode0,mode1,cix,isp,il,ib
C     local variable
      mode0= mod(mode,10)
      mode1=mod(mode/10,10)

      w=ommin+(ommax-ommin)*dble(iomg-1)/dble(nomg-1)
      do cix=1, ncix
         do isp=1,nsp
            do il=1, ldcix
               do ib=s_dmft%nlohi(1), s_dmft%nlohi(2)
                  if(w>0) then
                     if(mode0 .eq. 1) then
                        g=1./(w-(eigf(ib,iq,isp)+sigw(iomg,ib,iq,isp)))
                     else
                        g=sigw(iomg,ib,iq,isp)
                     endif
                  else
                     if(mode0 .eq. 1) then
                        g=1./(w-(eigf(ib,iq,isp)+dconjg(sigw(iomg,ib,iq,isp))))
                     else
                        g=dconjg(sigw(iomg,ib,iq,isp))
                     endif
                  endif
                  new=conjg(dmftu(ib,il,isp,cix))*g*dmftu(ib,il,isp,cix)
                  if (mode1 .eq. 0) then
                     output(il,isp,cix,iomg)=output(il,isp,cix,iomg)+new
                  else if (mode1 .eq. 1) then
                     output(il,isp,cix,iomg)=output(il,isp,cix,iomg)+g
                  else
                     call rx ('makesiggwloc.f : unkown mode')
                  endif
               enddo            !ib
            enddo               !ldcix
         enddo                  !cix
      enddo



      end subroutine
