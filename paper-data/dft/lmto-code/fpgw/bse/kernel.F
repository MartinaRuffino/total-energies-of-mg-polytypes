      subroutine BSE_kernel(qlat,ginv,mrecl,nl,nband,nlmto,nclass,
     .             natom,iclass,nspin,numkpts,nqibz,qbz,qibz,wbz,
     .             ngrp,symgg,irk,miat,tiat,invgx,shtvg,qsub,nblochpmx,
     .             nemlosp,Hmax2,nxx,nfilosp,ntp0sp,nt0sp,Hmax,nnc,nn,il
     .             ,in,im,nlnm,nlnmx,nlnmc,nlnmv,ncore,nctot,icore,
     .             min_diagval,opt_abs,kernel,kernel2,szk2,Hermit_H,npm,
     .    qsub_indx_ibz,switchr,EELS,nq0,Chi_RIXS,mat_rnk,conv_bnds,nwp,
     .    Chi_dir,alat)

C------------------------------------------------------------------------------
C
C       ISSUES!
C       is the wbz needed when adding V(0) to zw?
C       cgr calculated for all sym ops?
C       Symmetry acting on W(qibz) to get W(qbz)?
C       Units? 2*Rydberg to get kernel to eV?
C       Hermitian transpose of ppovlz in calculation below?
C       il and im etc or ilv etc in ppbafp_v2?
C
C-----------------------------------------------------------------------------

C qlat(or qbas) : reciprocal lattice vectors
C ginv:           inverse of qlat
C mrecl:          for reading WVR file. Read in from WV.d
C nl:             max number of l for all atoms
C nband:          Number of bands
C nlmto:          number of muffin tin orbitals. Depends on natom, the l
C                 value, the number of radial functions per l and the number of m per l
C nclass:         set equal to natom
C natom:          number of atoms
C iclass:         set to 1,2,3,4,....,natom
C nspin:          1 or 2 (needs work for nspin=2)
C numkpts(nqbz):  number of k-points in FBZ
C nqibz:          number of k-points in IBZ
C qbz:            k-points in FBZ
C qibz:           k-points in IBZ(note size of qibz is 3,numkpts)
C wbz:            weights for qbz. Seems that all elements are 1/numkpts
C ibzx:           for i in FBZ ibzx(i) gives the corresponding point in the IBZ
C ngrp:           Number of symmetry operations
C symgg:          Matrix representation of sym ops. For cubic sym there
C                 are 48 operations, e.g. xyz>xyz; xyz>xy-z; xyz>xzy etc
C irk:            irk(nqibz,ngrp) for each q in IBZ, irk for a given symmetry
C                 operation gives the point in the FBZ related by the given sym op
C miat:           Tells mapping between the atoms for the given sym op.
C tiat:           Actual shift in real space. pos(:,miat(atom,ig))+tiat(:,atom,ig)=pos_new
C invgx:          inverse. invg(3)=5 means inv(symgg(:,:,3))=symgg(:,:,5) (Or something to that effect)
C shtvg:          Shift in primitive lattice between the mapped atoms
C qsub:           determine chi(q), currently set to zero (optical limit)
C nblochpmx:      maximum number of Bloch functions and plane waves used/needed
C nemlo:          index to first unoccuppied state
C nfilo:          index to forst valence state used in calculation
C ntp0:           number of conduction states used
C nt0:            number of valence states used
C Hmax:           Size of kernel. number of cond states * num valence
C                 states * numkpts * nspin**2 (or just nspin?)
C nnc:            max number of radial functions (phi and phidot etc)
C                 used for the core
C nn:             Max number of radial functions (core + val) for all l
C il:             l=0,1,2,...,nl (value depends on what in and im are
C                 doing). Atom 1 comes first, then atom 2 etc
C in:             n=1,2,3 (phi, phidot etc) but the valence are added to
C                 the total number of core radial functions. e.g. if the number of core n for
C                 the given l is 2 then the first n for the valence will be 2+1=3 (depends on what im is doing)
C im:             m=-l,...0,...,+l
C nlnm:           size of il, in and im for each atom
C nlnmx:          Max of nlnm
C nlnmc:          number of core contribution to nlnm (for each atom)
C nlnmv:          number of valence contributions to nlnm (for each atom)
C ncore:          number of core states for each atom
C nctot:          total number of core states
C icore:          all elements set to 9999999. Then may stay this value
C                 depending on the value of incwfin passed to genallcf_v3. if incwfin=0 or
C                 -1, they will stay 9999999. If incwfin=-2 (which is used to get core exchange in hsfp0) the first
C                 nore(1) atoms are indexed from 1 up to ncore(1). Then all other
C                 indices remain 9999999 until we reach atom 2 etc and
C                 start from 1 up to ncore(2) etc.

      use m_readqg

      implicit none

      integer :: natom,sind,sindp,valind,valindp,condind,condindp,Hmax
      integer :: ifvcoud,iopen,iclose,ind1,ind2,ind3,ind4,ifrcw,szk2,
     .           mrecl,iatom,nspin,mdimx,nlnmx,Hmax2,n3_fin,n4_ini,
     .           iqdiff,iqdiffibz,sindb,sindpb,n3_num,n4_num,Hmax_new
      integer ::k_ind,numkpts,k_ind_ibz,nqibz,nfilo,nemlo,k_mat,kp_mat,
     .           k_indp,k_ind_ibzp,ngrp,irk(nqibz,ngrp),ig,npm,switchr,
     .           igp,k_indp_ibz,nblochpmx,miat(natom,ngrp),invgx(ngrp)
      integer :: nemlosp(*),nfilosp(*),ntp0sp(*),nt0sp(*),ary,recd
      integer ::ngb,ngc,nbloch,nbloch2,nl,nnc,nclass,nband,invr,invrot,
     .          iclass(*),nlmto,isp1,isp2,il(*),isp1p,nbloch3,
     .          in(*),im(*),nlnm(*),isp2p,numproc,procst,nq0,k_ind2
      integer(4) :: icore(*),numstates
      integer :: ncore(*),nlnmv(*),ngbkt,ngckt,ngbkmax,qsub_indx_ibz,
     .           nlnmc(*),nctot,ntp0,nt0,nxx,ngrpx,nn,qsub_indx,
     .           l2nl,irot,k_indfbz,k_indfbzp,ifidmlx,nq0ix,lxklm,
     .           qsub_indx_ibz2,mat_rnk,n3_ini,n4_fin,nwp
      integer,allocatable :: kdifftableIBZ(:,:),kdifftable(:,:),
     .                  lx(:),ngbk(:),ngck(:),imdim(:),iatomp(:),
     . nblocha(:),igvMPB(:,:),ibzx(:),sind_mat(:,:),sindp_mat(:,:),
     . iatomp_id(:)

      integer, allocatable :: nx(:,:),igv1(:,:),proc(:,:)
      real(8),allocatable::ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),ppb_id(:),
     .                      ppb(:),cgrx(:,:,:,:),ppbrd_id(:,:,:,:,:,:,:)
      real(8),allocatable:: dmlx(:,:),epinvq0i(:,:),wklm(:),epinv(:,:,:)

      character(5):: charnum5
      real(8) :: qvcoud(3),qlat(3,3),qlatinv(3,3),det,qibzk(3),quu(3),
     .        qsub(3),qbz(3,numkpts),ginv(3,3),wbz(numkpts),min_diagval,
     .           oneovernt0,oneoverntp0,fac,alat
      real(8) :: symope(3,3),shftv(3),qibz(3,nqibz),qvcoudk(3),fpi,
     .           symope_id(3,3),shftv_id(3)
      complex(8) :: kernel(Hmax2,Hmax2,mat_rnk),kernel2(szk2,szk2)
      complex(8),allocatable :: ZWBSE(:,:),ZWBSEZ(:,:)
      real(8) :: symgg(3,3,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp),
     .  qbzdiff(3),vc1
      real(8),allocatable :: gvecdiff(:),tr(:,:)
      complex(8) :: VBSE0,WBSE0,kern,kern2

      real(8),allocatable :: vcoud(:),Vcousq(:)
        integer, allocatable :: gvec(:,:)
      complex(8),allocatable :: zcousq(:,:),ppovl(:,:),VBSE(:,:),
     &                          ppovlz(:,:),WBSE(:,:,:),VBSE0_mat1(:),
     &                ppovlzkeep(:,:),zzmeltppovlz(:,:),VBSE0_mat2(:)
      complex(8),allocatable :: cphikc(:,:),cphikqc(:,:),
     .                         cphikv(:,:),cphikqv(:,:),zzmeltc(:,:,:),
     .                          zzmeltk(:,:),tmp1(:),zzmeltk_tmp(:),
     .                          cphikp(:,:),cphikqp(:,:),zw(:,:),
     .                 zzmelto(:,:,:),WBSE0_mat1(:,:),WBSE0_mat2(:,:),
     .                  WBSE0_mat3(:),zzmeltv(:,:,:),zzmeltcv(:,:,:)

      logical :: tiii,tiiiout,iprintx,Hermit_H,EELS,Chi_RIXS
      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),
     .                          igxt(:,:,:),eibzsym(:,:,:)
      integer,allocatable:: nrkip(:,:,:,:),irkip(:,:,:,:)
      logical :: lprt,Chi_dir,Chi_dir_pt
      logical :: opt_abs,conv_bnds

C The material below is for parallelizing over k' and k. Because we only
C calculate the upper triangular; for k'=1 we only need to calculate for
C k=1, for k'=2 we need k=1,2; k'=3, k=1,2,3 etc.  To distribute evenly
C we need to use triangular numbers! As an example, for 16 processors and 1000 k-points
C the distribution would be
C       1       250     31375 pairs to calculate
C       251     354     31460
C       355     433     31126
C       434     500     31289 etc
C       501     559
C       560     613
C       614     662
C       663     707
C       708     750
C       751     791
C       792     830
C       831     866
C       867     902
C       903     936
C       937     969
C       970     1000
C there are 1000*1001/2=500500 pairs of k',k to calculate. over 16
C processors this should be about 31281 pairs of k',k for each processor.
      Chi_dir_pt=.false.
      if(qsub_indx_ibz>nqibz.and.Chi_dir)Chi_dir_pt=.true.
      numproc=1 ! Need to determine in the future
      allocate(proc(numproc,2))
      proc(1,1)=1
      do ind1 = 2, numproc
      procst=0.5d0*(-1.0d0+(1+4.0d0*(ind1-1)/numproc*(numkpts+1.0d0)
     .          *numkpts)**0.5d0)
        proc(ind1,1)=procst+1
        proc(ind1-1,2)=procst
      enddo
      proc(numproc,2)=numkpts*numkpts/2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      lprt = .false.
      fpi = 16.0d0*datan(1.0d0)
      shftv_id=(/0d0,0d0,0d0/)
      symope_id(1:3,1) = (/1d0,0d0,0d0/)
      symope_id(1:3,2) = (/0d0,1d0,0d0/)
      symope_id(1:3,3) = (/0d0,0d0,1d0/)
      allocate(iatomp_id(natom))
      do ind1=1,natom
        iatomp_id(ind1) = ind1
      enddo
!      call dinv33(qlat,0,qlatinv,det)
      call minv33(qlat,qlatinv)

C We need to get wklm: coefficient w_L in Eqns. 28,33 of JPSP 83, 094711 (2014) for V(0)
C      call iqindx2qg(qsub,1,qsub_indx,quu)
      allocate(ibzx(numkpts))
      call invkibzx(irk,nqibz,ngrp,numkpts,ibzx)
C      qsub_indx_ibz2 = ibzx(qsub_indx)

C      write(*,*)qsub_indx_ibz2,qsub_indx_ibz

        if(qsub_indx_ibz==1)then
        ifidmlx = iopen('EPSwklm',0,0,0)
         read(ifidmlx) nq0ix,lxklm
         allocate( dmlx(nq0ix,9))
         allocate( epinvq0i(nq0ix,nq0ix),epinv(3,3,nq0ix) )
         allocate( wklm((lxklm+1)**2))

         read(ifidmlx) dmlx, epinv,epinvq0i
         write(*,*)'epinvq0i,lxklm=',epinvq0i,lxklm
         read(ifidmlx) wklm
        ifidmlx = iclose('EPSwklm')
        deallocate(dmlx,epinvq0i,epinv)

        allocate(zw(nblochpmx,nblochpmx))! This is how zw=W-V is written to file
        allocate(WBSE(nblochpmx,nblochpmx,nqibz)); WBSE=0d0

        endif ! qub_indx_ibz==1

        allocate(ngbk(nqibz+nq0),ngck(nqibz+nq0)) ! number of basis functions and number of plane waves
        ngbk = 0; ngck =0

      do k_ind = 1, nqibz+nq0

      if(qsub_indx_ibz==1)goto 524  ! Only need to calculate W_IJ for q=0, since W is independent of q: W(k'-k)
      if(k_ind/=qsub_indx_ibz)cycle
 524  continue
      if(qsub_indx_ibz <= nqibz.and.k_ind>nqibz)cycle
      k_ind2=k_ind
      if(k_ind>nqibz)k_ind2=1
      ifvcoud = iopen('Vcoud.'//charnum5(k_ind2),0,0,0)
      read(ifvcoud) ngbkt
      read(ifvcoud) qvcoudk
      if(allocated(vcoud)) deallocate(vcoud)
      if(allocated(zcousq)) deallocate(zcousq)
      allocate(zcousq(ngbkt,ngbkt),vcoud(ngbkt))
      read(ifvcoud) vcoud
      read(ifvcoud) zcousq

      if(sum(abs(qibz(:,k_ind2)-qvcoudk))>1.0E-07) call rx('BSE: wrong
     .                  qibz from Vcoud file')

      ngbk(k_ind) = ngbkt
      if (k_ind==qsub_indx_ibz) ngb=ngbkt

C     --- Make overlap matrix (Step A, Remarks) ---
!!        <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!!        zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!!        ppovlz= O Zcousq
!!        (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.

      if(k_ind>nqibz)qvcoudk=qsub
      call getngvec('QGcou',qvcoudk,qlatinv, quu,ngckt)
      ngck(k_ind)=ngckt
      if(k_ind==qsub_indx_ibz)ngc=ngckt

      nbloch = ngbkt - ngckt

C      if (nbloch /= nbloch2 ) call rx('BSE: wrong nbloch,
C     &                  nbloch=',nbloch,'nbloch2=',nbloch2)

C overlap of plane wave basis!
      if(allocated(ppovl))deallocate(ppovl)
      allocate(ppovl(ngckt,ngckt))
        write(*,*)'ngckt',ngckt
      call readppovl0(qvcoudk,ngckt,ppovl)
      if(allocated(ppovlz))deallocate(ppovlz)
      allocate(ppovlz(ngbkt,ngbkt))
      ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)! up to nbloch (MTO's), the overlap is the identity

      call zgemm('N','N',ngckt,ngbkt,ngckt,(1d0,0d0),ppovl,ngckt,
     & zcousq(nbloch+1,1),ngbkt,(0d0,0d0),ppovlz(nbloch+1,1),
     & ngbkt)

      if(qsub_indx_ibz==1)then
      ifrcw = iopen('WVR.'//charnum5(k_ind),0,-1,mrecl)
      recd=1
      write(*,*)nwp
      if(npm==2)recd=nwp ! when T.R.S. off, the w=0 part is not the first!
        write(*,*)'recd=',recd
      read(ifrcw,rec=recd) zw! W(0) - V, i.e., for omega=0! Maybe read in rec as is done in /gwsrc/sxcf.sc.F ?
      do ind1 = 1, ngbkt
C zw is W-V. We are adding V here to get W

        if(ind1==1.and.k_ind==1)then
           zw(ind1,ind1)=zw(ind1,ind1)+wklm(1)*fpi*sqrt(fpi)/wbz(1) ! is the wbz needed here? Is this cancelled with when we divide by N?
C zw for k=0 and \mu=0 is calculated in /gw/fpgw/main/hx0fp0.sc.m.F, but we need to add V.
C V_{\mu=0}(k=0) is taken from /gw/fpgw/gwsrc/sxcf.sc.F and mptauof.F
C getw0l
        else
           zw(ind1,ind1) = zw(ind1,ind1) + vcoud(ind1) ! in basis that diagonalizes V
        endif
      enddo
      endif ! qsub_indx_ibz == 1
      if(switchr==0.and.k_ind==1) vcoud(1) = 0d0 ! We dont include the long range part (G=0) of V for optical absorption
      if(switchr /= 0.and.k_ind==1)Vc1 = wklm(1)*fpi*sqrt(fpi)/wbz(1)

      if(k_ind==qsub_indx_ibz.and.switchr /= 0)then
         allocate(ppovlzkeep(ngbkt,ngbkt))
         ppovlzkeep = ppovlz
         allocate(Vcousq(ngb))
         DO ind1=1,ngb
           Vcousq(ind1)=sqrt(Vcoud(ind1))
         ENDDO
        open(2892,file='XVcousq',form='unformatted')
        rewind 2892
        write(2892)ngb
        if(qsub_indx_ibz==1)write(2892)Vc1
        write(2892)Vcoud(1:ngb)

      endif
C--------The stuff below is probably not ever going to be used!!!
      if(switchr /= 0.and.k_ind==1)vcoud(1)=wklm(1)*fpi*sqrt(fpi)/
     .          wbz(1)
      if(switchr <= 2.and.k_ind==1)vcoud(1)=0.0d0 ! IS THIS THE CASE for calculating eps_{mu=nu=1}(q=0)
      if(switchr==0.and.k_ind==1.and.EELS) vcoud(1)=
     .           wklm(1)*fpi*sqrt(fpi)/wbz(1)
      if(switchr==3.and.k_ind>nqibz)vcoud(1)=
     .          0.0d0!fpi/sum(qsub**2*(0.5d0*fpi/alat)**2)
C      if(.not.opt_abs.and.k_ind==1)Vcousq(1)=sqrt(vcoud(1))
C------------------------------------------------------------------
      if(k_ind==qsub_indx_ibz)then

      allocate(VBSE(ngb,ngb))
      do ind1=1,ngbkt
        do ind2=1,ngbkt

         VBSE0 = 0d0

         do ind3=1,ngbkt
             VBSE0 = VBSE0 + ppovlz(ind1,ind3) * vcoud(ind3) *
     &              dconjg(ppovlz(ind2,ind3))
         enddo

         VBSE(ind1,ind2) = VBSE0

        enddo
      enddo

      endif
C      if(allocated(vcoud)) deallocate(vcoud)

      if(qsub_indx_ibz==1)then
      if(allocated(WBSE0_mat1))deallocate(WBSE0_mat1)
        allocate(WBSE0_mat1(ngbkt,ngbkt)); WBSE0_mat1 = 0.0d0

        call zgemm('N','N',ngbkt,ngbkt,ngbkt,dcmplx(1d0,0d0),ppovlz,
     &      ngbkt,zw,nblochpmx,dcmplx(0d0,0d0),WBSE0_mat1,ngbkt)
        call zgemm('N','C',ngbkt,ngbkt,ngbkt,dcmplx(1d0,0d0),WBSE0_mat1,
     &    ngbkt,ppovlz,ngbkt,dcmplx(0d0,0d0),WBSE(1:ngbkt,1:ngbkt,k_ind)
     &    ,ngbkt)

      endif
      if(qsub_indx_ibz==1)ifrcw = iclose('WVR.'//charnum5(k_ind))
      ifvcoud = iclose('Vcoud.'//charnum5(k_ind2))
      enddo ! End of the loop over IBZ

      if(allocated(vcoud)) deallocate(vcoud)
      if(allocated(zcousq)) deallocate(zcousq)
      if(allocated(ppovl))deallocate(ppovl)
      if(allocated(ppovlz))deallocate(ppovlz)
      if(allocated(wklm))deallocate(wklm)
      if(qsub_indx_ibz==1)deallocate(WBSE0_mat1,zw)!,WBSE0_mat2
C      if (rank == 0)

C Currently I do not us the EIBZ, however, this will probably be
C employed later!

C irkip is such that for each k point, we have a subset from irk.  For
C example imagine the third dimension in irkip runs over the FBZ, and we take index
C 4, which for Si 8*8*8 is 0.375(-1,1,1). If we look at the first
C dimension in irk and take the second index, k_IBZ(2)=0.125(1,1,-1).

C irk tells the mapped point in FBZ from the given k_IBZ and SYMOP.
C irk(2,1)=65 means S1*k_IBZ=k(65), and irk(2,2)=449; and S2*k_IBZ+(1,1,-1)=k(449) etc.
C Notice a G vector to bring back into the 1st BZ may be required.
C When we use the these symops and operate on the third dimension, for example k(4), we end up
C with DIFFERENT rotated points, in the example for k_IBZ=2 and k'=4 for
C Si, we have four DIFFERENT rotated points.  Therefore we need to work out
C quantities for these four different rotations.

C nrkip tells how many points the irkip rotates into.

C       allocate(nwgt(numkpts,1:numkpts),igx(ngrp*2,numkpts,numkpts),
C     .          igxt(ngrp*2,numkpts,numkpts)
C     .          , eibzsym(ngrp,-1:1,numkpts))
C        tiii=.false.    !Enforce no time reversal. time reversal not yet...
C        print *,'NOTE:TimeReversal not yet implemented in BSE'
C         iprintx=.false.
C         call eibzgen(numkpts,symgg,ngrp,qbz,
C     &        1,numkpts,qbz,numkpts,tiii,ginv,iprintx,
C     o        nwgt,igx,igxt,eibzsym,tiiiout)
C
C        allocate(irkip(nqibz,ngrp,numkpts,nspin))
C        allocate(nrkip(nqibz,ngrp,numkpts,nspin))
C         nrkip=0
C         irkip=0
C         isp1=1                   ! not spin dependent
C         do k_ind=1,numkpts
C           do k_indp=1,nqibz
C             do ig=1,ngrp
C               k_indfbz = irk(k_indp,ig)
C               if(k_indfbz==0) cycle
C               if(nwgt(k_indfbz,k_ind)/=0) then
C                irkip(k_indp,ig,k_ind,isp1)= irk(k_indp,ig)
C                nrkip(k_indp,ig,k_ind,isp1)= nwgt(k_indfbz,k_ind)
C               endif
C             enddo !igrp
C           enddo !k_indp
C         enddo !k_ind
C
C----------------------------------------------------------------------

C The stuff below needs rewritten for the cases: 1)T.R.S. off!,
C 2) Going beyond the Tamm-Dancoff approx,
C 3)

C The augmentation part of the quantity <psi|psi M> = <phi1|phi2 B> is
C independent of k, as these are just Bloch summed over later? l2nl is
C an initial guess for the maximum l for MPB. nx is number of radial functions for a given l for the MPB
C and nxx is the maximum nx. lx is max l value for each atom. nblocha is
C number of Bloch functions in MPB for each atom. ppbrd(L on phi1,n index
C on phi1,L on phi2,n on phi2,l from MPB,n from MPB,atom&spin). cgr are
C the Clebsch Gordan coefficients.

      ngrpx = 1
      l2nl=2*(nl-1)
      allocate( nx(0:l2nl,nclass),lx(nclass),
     &     ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:l2nl,nxx, nspin*nclass),
     &     nblocha(natom),cgr(nl**2,nl**2,(2*nl-1)**2,ngrp))
      if(qsub_indx_ibz/=1)then
      allocate(ppbrd_id(0:nl-1,nn, 0:nl-1,nn, 0:l2nl,nxx, nspin*nclass),
     &     cgrx(nl**2,nl**2,(2*nl-1)**2,ngrpx))
C-----------------------------------------------------------------
C
C   Maybe dont need ppbrd_id and just use cgr(:,:,:,1) below?
C
C-----------------------------------------------------------------
       call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope_id,qlat,
     o     nblocha, lx, nx,  ppbrd_id , mdimx, nbloch3, cgrx)
      endif
      call rdpp_v3(nxx, nl,ngrp, nn, nclass, nspin,symgg,qlat,
     o     nblocha, lx, nx,  ppbrd , mdimx, nbloch2, cgr)

      if (mdimx /= maxval(nblocha) ) call rx('BSE: wrong mdimx')

      allocate(imdim(natom),iatomp(natom))

      call indxmdm(nblocha,nclass,iclass,natom,imdim) ! Index to first Bloch state for the particular atom

      WRITE(*,*)'If the code crashes here, maybe turn KeepPPOVL off in
     .           GWinput?'
      ngbkmax=nblochpmx!maxval(ngbk)

      allocate(ppb(nlnmx*nlnmx*mdimx*nclass))
      allocate(tr(3,natom))
        ary=1
        if(mat_rnk>1)ary=2

      allocate(zzmeltk(ngb,Hmax*ary))
      zzmeltk=0d0
      Hmax_new=Hmax
      if(switchr /= 0)open(2891,file='ZWBSE',
     .          form='unformatted')
      if(switchr /= 0)rewind 2891
      do isp1p = 1,nspin
      write(*,*)'spin=',isp1p
      nfilo = nfilosp(isp1p)
      nt0 = nt0sp(isp1p)
      nemlo = nemlosp(isp1p)
      ntp0 = ntp0sp(isp1p)
C      n3_fin = nfilosp(isp1p)+nt0sp(isp1p)-1
C      n4_ini = nemlosp(isp1p)
C      n3_num = nt0sp(isp1p)
C      n4_num = ntp0sp(isp1p)
      n3_ini=nfilo
      n3_fin=nfilo+nt0-1+ntp0*(ary-1)
      n4_ini=nemlo-nt0*(ary-1)
      n4_fin=nemlo+ntp0-1
      n3_num=nt0+ntp0*(ary-1)
      n4_num=ntp0+nt0*(ary-1)

      allocate(cphikqv(nlmto,n3_ini:n3_fin),
     .          cphikv(nlmto,n3_ini:n3_fin))
      allocate(cphikqc(nlmto,n4_ini:n4_fin),
     .          cphikc(nlmto,n4_ini:n4_fin))

      if(qsub_indx_ibz==1)
     .   allocate(zzmeltv(ngbkmax,1:n3_num,1:n3_num))
      if(qsub_indx_ibz==1)zzmeltv=0d0
      allocate(zzmeltc(ngbkmax,1:n4_num,1:n4_num));zzmeltc=0d0
      allocate(ZWBSE(ngbkmax,n3_num*n3_num));ZWBSE=0d0
      allocate(ZWBSEZ(n3_num*n3_num,n4_num*n4_num));ZWBSEZ=0d0
      allocate(zzmeltcv(1:ngbkmax,n4_ini:n4_fin,n3_ini:n3_fin))
      zzmeltcv = 0d0

C ----------------------------------------------------
C For W we have <v'k'|vk M_I>, which is indexed over I and then over
C v'v, with the same for <c'k'|ck M_J>^*.  The kernel is indexed over
C cv.  The matrices below give the transformation, i.e., if you know the index
C for vv' and cc', the matrices below give the indices for the kernel.
      oneovernt0=1d0/n3_num
      oneoverntp0=1d0/n4_num

      allocate(sind_mat(n3_num*n3_num,n4_num*n4_num),
     .         sindp_mat(n3_num*n3_num,n4_num*n4_num))
C-------------------------------------------------------
C
C
C
C------------------------------------------------------
      DO sindpb=1,n4_num*n4_num
        DO ind2=1,n4_num
          if(mod(sindpb-ind2,n4_num)==0)exit
        ENDDO
        condindp=ind2
        condind=(sindpb-ind2)*oneoverntp0+1
        DO sindb=1,n3_num*n3_num
        DO ind1=1,n3_num
          if(mod(sindb-ind1,n3_num)==0)exit
        ENDDO
        valindp=ind1
        valind=(sindb-ind1)*oneovernt0+1
        sindp = (valindp-1)*ntp0+condindp
        sind = (valind-1)*ntp0+condind
        if(ary>1)then
        sind=0;sindp=0
        if(valind <= nt0.and.valindp <= nt0.and.condind>nt0.and.
     .          condindp>nt0)then
        sindp = (valindp-1)*ntp0 +condindp-nt0
        sind = (valind-1)*ntp0 + condind-nt0
        elseif(valind>nt0.and.valindp>nt0.and.condind <= nt0.and.
     .          condindp <= nt0)then
        sindp = (condindp-1)*ntp0+valindp-nt0+Hmax ! The index valind(condind) is now conduction(valence) states
        sind = (condind-1)*ntp0 + valind-nt0+Hmax
        endif
        endif
        sind_mat(sindb,sindpb)=sind
        sindp_mat(sindb,sindpb)=sindp
      ENDDO
      ENDDO
C-----------------------------------------------------------
C
C
C
C-----------------------------------------------------------
      do k_indp = 1, numkpts
        if(mod(k_indp-1,20)==0)write(*,*)'k=',k_indp,'of',numkpts

          call getcphi(qbz(:,k_indp),nlmto,nfilo,n3_fin,
     .                  isp1p,quu,cphikqv)
          isp2p=isp1p
CC Spin-flip?
CC       DO isp2p = 1, nspin
          call getcphi(qbz(:,k_indp)+qsub,nlmto,n4_ini,n4_fin,
     .                  isp2p,quu,cphikqc)
         kp_mat = (isp1p-1)*numkpts*nt0sp(1)*ntp0sp(1)+(k_indp-1)*(
     .                 nt0*ntp0)

          do k_ind = 1,k_indp!numkpts!k_indp ! Since W is Hermitian. For k' k we have
CC                             ! all states, then when we swap we already have what we need!
          if(conv_bnds.and.k_ind /= k_indp)cycle
          do isp1 = 1, nspin
          if(isp1/=isp1p)cycle ! W is diagonal in the spins

            qbzdiff = qbz(:,k_indp)-qbz(:,k_ind)

            call iqindx2qg(qbzdiff,2,iqdiff,quu)

            iqdiffibz = ibzx(iqdiff)
            call getngvec('QGcou',qibz(:,iqdiffibz),qlatinv, quu,ngckt)
            ngbkt = nbloch + ngckt
         do irot = 1, ngrp
           if(irk(iqdiffibz,irot)==iqdiff)then; ig=irot;exit;endif
         enddo
CC --- The stuff below has been taken from hsfp0.sc.m.F
         call ppbafp_v2(ig,ngrp,isp1p,nspin,il,in,im, nlnm,
     .    nl,nn,nclass,nlnmx,mdimx,lx,nx,nxx,  !Bloch wave
     .    cgr,nl-1,         !rotated CG
     .    ppbrd,            !radial integrals
     o    ppb)

            invr   = invrot(ig,invgx,ngrp) ! index to the inverse of rotation ig above
            tr     = tiat(:,:,invr)
            iatomp = miat(:,invr)
            symope = symgg(:,:,ig)
            shftv  = matmul(symope,shtvg(:,invr))
C ----------------------------------------------------

         IF(qsub_indx_ibz==1.and.switchr <= 3)then
          call getcphi(qbz(:,k_ind),nlmto,nfilo,
     .                  n3_fin,isp1,quu,cphikv)
         call mkmelt(1,isp1p,isp1,qlat,qbz(:,k_indp),qbzdiff,! For the lower right quadrant of A (Eqn B.20 Onida Rev. Mod. Phys) do we need k'+q here instead of just k'?
     .            qibz(:,iqdiffibz),cphikqv(1:nlmto,nfilo:n3_fin),
     .            cphikv(1:nlmto,nfilo:n3_fin),ppb,icore,ncore,
     .            symope,shftv,tr,nl,natom,nclass,iclass,iatomp,
     .            nlnmx,nlnmv,nlnmc,imdim,nblocha,mdimx,nlmto,
     .         ngckt,nbloch,nfilo,n3_num,[0],npm,nfilo,n3_num,1,ngbkt,1,
     .      n3_num,1,n3_num,nnc,lprt,zzmeltv(1:ngbkt,1:n3_num,1:n3_num))

       call zgemm('T','N',ngbkt,n3_num*n3_num,ngbkt,dcmplx(1d0,0d0),(
     . WBSE(1:ngbkt,1:ngbkt,iqdiffibz:iqdiffibz)),ngbkt,
     . zzmeltv(1:ngbkt,1:n3_num,1:n3_num),ngbkt,dcmplx(0d0,0d0),
     . ZWBSE(1:ngbkt,1:n3_num*n3_num),ngbkt) ! zzmeltv, combines dimensions 2 and 3 into one, and increments over the 2nd first

          if(switchr /= 0)write(2891)ZWBSE(1:ngbkt,1:n3_num*n3_num)
        ELSE
        read(2891)ZWBSE(1:ngbkt,1:n3_num*n3_num)
        ENDIF

        isp2 = isp1
CC Spin-flip?
CC       DO isp2 = 1, nspin
          call getcphi(qbz(:,k_ind)+qsub,nlmto,n4_ini,
     .                  n4_fin,isp2,quu,cphikc)

          call mkmelt(1,isp2p,isp2,qlat,qbz(:,k_indp)+qsub,qbzdiff,! For the lower right quadrant of A (Eqn B.20 Onida Rev. Mod. Phys) do we need just k' here instead of k'+q?
     .          qibz(:,iqdiffibz),cphikqc(1:nlmto,n4_ini:nemlo+ntp0-1)
     .          ,cphikc(1:nlmto,n4_ini:nemlo+ntp0-1),ppb,icore,ncore,
     .          symope,shftv,tr,nl,natom,nclass,iclass,iatomp,
     .          nlnmx,nlnmv,nlnmc,imdim,nblocha,mdimx,nlmto,
     .       ngckt,nbloch,n4_ini,n4_num,[0],npm,n4_ini,n4_num,1,ngbkt,1,
     .      n4_num,1,n4_num,nnc,lprt,zzmeltc(1:ngbkt,1:n4_num,1:n4_num))

        call zgemm('T','N',n3_num*n3_num,n4_num*n4_num,ngbkt,
     . dcmplx(1d0,0d0),ZWBSE(1:ngbkt,1:n3_num*n3_num),ngbkt,
     . dconjg(zzmeltc(1:ngbkt,1:n4_num,1:n4_num)),ngbkt,dcmplx(0d0,0d0),
     . ZWBSEZ(1:n3_num*n3_num,1:n4_num*n4_num),n3_num*n3_num)

      k_mat = (isp1-1)*numkpts*nt0sp(1)*ntp0sp(1)+(k_ind-1)*nt0*ntp0
CC Spin-flip
      DO sindpb = 1, n4_num*n4_num; DO sindb = 1, n3_num*n3_num

      if(sind_mat(sindb,sindpb)>0.and.sindp_mat(sindb,sindpb)>0)then
      if(sind_mat(sindb,sindpb) <= Hmax.and.sindp_mat(sindb,sindpb)
     . <= Hmax)then
      kernel(sind_mat(sindb,sindpb)+k_mat,sindp_mat(sindb,sindpb)+
     .          kp_mat,1) = ZWBSEZ(sindb,sindpb)
      elseif(sind_mat(sindb,sindpb)>Hmax.and.sindp_mat(sindb,sindpb)>
     .          Hmax)then
      kernel(sind_mat(sindb,sindpb)-Hmax+k_mat,sindp_mat(sindb,sindpb)-
     .          Hmax+kp_mat,2) = ZWBSEZ(sindb,sindpb)
      endif
      endif
      ENDDO; ENDDO

           if (k_indp==k_ind.and.isp1p==isp1)then
CC Spin-flip?
CC          if(k_indp==k_ind.and.isp1p==isp1.and.isp2p==isp2)then
        if(qsub_indx_ibz/=1)then
         call ppbafp_v2(1,ngrpx,isp1p,nspin,il,in,im, nlnm,
     .    nl,nn,nclass,nlnmx,mdimx,lx,nx,nxx,  !Bloch wave
     .    cgrx,nl-1,         !rotated CG
     .    ppbrd_id,            !radial integrals
     o    ppb)
        endif

           call mkmelt(0,isp2p,isp1p,qlat,qbz(:,k_indp)+qsub,qsub,
     .         qsub,cphikqc(1:nlmto,n4_ini:nemlo+ntp0-1)
     .         ,cphikqv(1:nlmto,nfilo:n3_fin),ppb,icore,ncore,
     .       symope_id,shftv_id,[0d0],nl,natom,nclass,iclass,iatomp_id,
     .         nlnmx,nlnmv,nlnmc,imdim,nblocha,mdimx,nlmto,ngc,nbloch,
     .         n4_ini,n4_num,[0],npm,nfilo,n3_num,1,ngb,1,n3_num,1,
     .         n4_num,nnc,lprt,
     .         zzmeltcv(1:ngb,n4_ini:n4_fin,n3_ini:n3_fin))

           do valindp = nfilo, nfilo+nt0-1
           do condindp = nemlo, nemlo+ntp0-1

            sindp = (isp1p-1)*numkpts*nt0sp(1)*ntp0sp(1)+
     .                  (k_indp-1)*nt0*ntp0+(valindp-nfilo)*ntp0+
     .                  condindp-nemlo+1
          zzmeltk(1:ngb,sindp) = zzmeltcv(1:ngb,condindp,valindp)
          if(ary>1)
     .         zzmeltk(1:ngb,sindp+Hmax) = zzmeltcv(1:ngb,valindp,condindp)

             enddo;enddo
          ENDIF

CC Spin-flip?
CC          enddo ! isp2
          enddo ! isp1
          enddo ! k
C Spin-flip?
C        enddo ! isp2p
      enddo ! k'
      if(qsub_indx_ibz==1)deallocate(zzmeltv)
      deallocate(zzmeltc)
      deallocate(cphikc)
      deallocate(cphikv)
      deallocate(sind_mat)
      deallocate(sindp_mat)
      deallocate(ZWBSE)
      deallocate(ZWBSEZ)
      deallocate(zzmeltcv)
      deallocate(cphikqc)
      deallocate(cphikqv)
      enddo ! isp1p
      deallocate(nx)
      deallocate(lx)
      deallocate(ppbrd)
      deallocate(cgr)
      deallocate(ppb)
      if(qsub_indx_ibz==1)deallocate(WBSE)
      if(qsub_indx_ibz/=1)deallocate(ppbrd_id,cgrx)
      write(*,*)'W calculated'

      if(switchr /= 0)CLOSE(2891)

        if(nspin==1.and.switchr==0)then        ! THIS WILL NEED CHECKING?
        fac = -0.5d0
        else
        if(nspin==1)then
           fac = -0.5d0
        else;fac=-1.0d0
        endif
       endif
      kern=0.0d0;kern2=0.0d0
      ALLOCATE(VBSE0_mat1(1:ngb),VBSE0_mat2(1:ngb))
      kernel =fac*kernel
      if(switchr /= 3.or.Chi_dir_pt)then
      do sindp = 1, Hmax_new
C        if(switchr <= 2.or.Chi_RIXS)then
        if(mod(sindp,300)==0)write(*,*)'s= ',sindp,' of ',Hmax_new
          call zgemv('N',ngb,ngb,dcmplx(1d0,0d0),VBSE,ngb,
     &      dconjg(zzmeltk(:,sindp)),1,dcmplx(0d0,0d0),VBSE0_mat1,1)
        if(ary>1)call zgemv('N',ngb,ngb,dcmplx(1d0,0d0),VBSE,ngb,
     &     dconjg(zzmeltk(:,sindp+Hmax)),1,dcmplx(0d0,0d0),VBSE0_mat2,1)
C        endif
       do sind = 1,sindp
C          if(switchr <= 2.or.Chi_RIXS)then
             call zgemm('T','N',1,1,ngb,dcmplx(1d0,0d0),zzmeltk(:,sind),
     &                  ngb,VBSE0_mat1,ngb,dcmplx(0d0,0d0),kern,1)
          if(ary>1)call zgemm('T','N',1,1,ngb,dcmplx(1d0,0d0),
     &  zzmeltk(:,sind+Hmax),ngb,VBSE0_mat2,ngb,dcmplx(0d0,0d0),kern2,1)
C          endif
          kernel(sind,sindp,1) = kernel(sind,sindp,1)  + kern
          if(ary>1)kernel(sind,sindp,2)= kernel(sind,sindp,2)+kern2
C        if(sind /= sindp)kernel(sindp,sind) = dconjg(kernel(sind,sindp))
        enddo
      enddo
      write(*,*)'V calculated'
      endif
      DEALLOCATE(VBSE0_mat1,VBSE0_mat2)
      IF(switchr /= 0)then
        allocate(zzmeltppovlz(Hmax*ary,ngb))
       call zgemm('T','N',Hmax*ary,ngb,ngb,dcmplx(1d0,0d0),zzmeltk(
     .1:ngb,1:Hmax*ary),ngb,ppovlzkeep(1:ngb,1:ngb),ngb,dcmplx(0d0,0d0),
     .  zzmeltppovlz,Hmax*ary)
        deallocate(ppovlzkeep)
      DO ind1 = 1, ngb; DO ind2 = 1, Hmax*ary
        zzmeltppovlz(ind2,ind1)=zzmeltppovlz(ind2,ind1)*
     .          Vcousq(ind1)  ! eps_(mu,nu) = 1 - sqrt(V_mu) * X_(mu,nu) * sqrt(V_nu)
      ENDDO; ENDDO
       write(2892)zzmeltppovlz
       CLOSE(2892)

      deallocate(zzmeltppovlz,Vcousq)
      endif

      deallocate(ibzx)
      deallocate(tr)
      deallocate(imdim,iatomp,iatomp_id,nblocha,zzmeltk)
      deallocate(ngck,ngbk)
      deallocate(VBSE,proc)

      end
