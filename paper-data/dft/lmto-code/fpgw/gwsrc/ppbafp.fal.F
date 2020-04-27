C#define OPENMP
      subroutine ppbafp_v2 (ig,ng,isp,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd,           !radial integrals
     o  ppb)
C- Calculates for all R (class) <Phi(RLn) Phi(RL'n') B(R,i)>
c iln            = s. indxln.f
c nn,nl          = number of n and l
c
c mdim           = dimension of orthonormal basis
c mnl            = total number of partial waves (each with given (m,n,l)) for a given class
c mnlx           = dimensions il,in,im
c in,il,im      = index for n,l,m s. indxlnm.f
c
c ppb            = <Phi(RLn) Phi(RL'n') B(R,i)>
c
c  Parameter n differentiates core phi phidot localOrbital.
c
c B(R,i)   = orthonormal basis, i=1,mdim
c
!     use mpi_mod, rank => mpi_rank

      implicit none
      integer(4):: ic,nclass,nl,nn,mnlx,mdimx,
     &             isp, nxx,lmxax,
     &             i,lb,nb,mb,lmb,i1,ibas,i2,ig,
     &             np,lp,mp,lmp,n,l,m,lm,nspin,ng
      integer(4)::
     i          il(mnlx,nclass),in(mnlx,nclass),im(mnlx,nclass),
     i          mnl(nclass), lx(nclass),nx(0: 2*(nl-1),nclass)
c

      real(8), intent(out) :: ppb(mnlx,mnlx,mdimx,nclass)
      real(8) ::
     &   cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng),
     &   ppbrd( 0:nl-1,nn, 0:nl-1,nn, 0:2*(nl-1), nxx, nclass*nspin)
      integer(4):: zvztest

C     integer omp_get_thread_num
ccccccccccccccccccc
c      write(6,*)'    ppbafp: sum(ppbrd)',sum(ppbrd)
c      ppb=0d0;pdb=0d0;dpb=0d0;ddb=0d0
c      ppb=1d99;pdb=1d99;dpb=1d99;ddb=1d99
cccccccccccccccccccc


!ML      write(0,'("[",I4,"] ppbafp_v2 called.")') rank
!ML      write(6,'("[",I4,"] ppbafp_v2 called.")') rank
!ML      write(6,'("[",I4,"] ppbafp_v2: sum(pprd) = ",E11.3)') rank,sum(ppbrd)

c loop over class  ! nclass=nbas
!$OMP parallel do  private (ic,ibas,i,lb,nb,mb,lmb,i2,np,lp,mp,lmp,i1,n,l,m,lm)
      do ic  = 1, nclass

        ibas = ic
c        do ip1=1,2
c        do ip2=1,2
c          ioff(ip1,ip2) =  ip1+ 2*(ip2-1)+ 4*(ic-1) + 4*nclass*(isp-1)
c        enddo
c        enddo
c     print *, ' ic lx=',ic, lx(ibas)

!ML tests:

!        write(6,'("[",I4,"] ppbafp: ibas = ",I3," lx(ibas) = ",I3,", nx = ",20I4)')
!     &       rank, ibas, lx(ibas), (nx(lb,ibas),lb=1,lx(ibas))

!        write(6,'("[",I4,"] ppbafp: nclass = ",I3," mnl = ",20I4)')
!     &       rank, nclass,  (mnl(i2),i2=1,nclass)

!        write(6,'("[",I4,"] ppbafp: ic = ",I3,", in = ",30I3)')
!     &       rank, ic,  (in(i2,ic),i2=1,min(mnl(ic),30))

!        write(6,'("[",I4,"] ppbafp: ic = ",I3,", il = ",30I3)')
!     &       rank, ic,  (il(i2,ic),i2=1,min(mnl(ic),30))

!        write(6,'("[",I4,"] ppbafp: ic = ",I3,", im = ",30I3)')
!     &       rank, ic,  (im(i2,ic),i2=1,min(mnl(ic),30))

!ML tests end


!! loop over product basis
        i = 0 !i = product basis index.
        do lb  = 0, lx (ibas)
          do nb  = 1, nx (lb,ibas)
            do mb  = -lb, lb
              i    = i+1  !The number of product basis is  =(i at the end of loop).
              lmb  = lb*lb + lb + mb + 1

cccccccccccccccccc
c              write(6,"(' === product basis index: i ibas l n m',10i4)")i,ic,lb,nb,mb
cccccccccccccccccc

c loop over LMTO basis
              do  i2 = 1,mnl(ic)
                np   = in(i2,ic)
                lp   = il(i2,ic)
                mp   = im(i2,ic)
                lmp  = lp*lp + lp + mp + 1

                do  i1 = 1,mnl(ic)
                  n    = in(i1,ic)
                  l    = il(i1,ic)
                  m    = im(i1,ic)
                  lm   = l*l + l + m + 1
cccccccccccccccccccccccccccccccc
c        do j1=1,2;doj2=1,2
c        if(abs(ppbrd(l,n, lp,np, lb,nb, ioff(j1,j2))) >1d10) then
c          print *, abs(ppbrd( l,n, lp,np, lb,nb, ioff(j1,j2)))
c          print *, l,n, lp,np, lb,nb, ioff(j1,j2)
c          stop
c        endif
c       enddo; enddo
ccccccccccccccccccccccccccccccccccccc
                  ppb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
     &                  * ppbrd( l,n, lp,np, lb,nb, isp+nspin*(ic-1))

cccccccccccccccccccccccccccccccccccccccccccccccccc
                  if(zvztest()==2) then
                    if(abs(ppb(i1,i2,i,ic))/=0d0  )then
                      write(6,"('i1 i2 i ic=',4i4,' ppb=',f10.6)") i1,i2,i,ic, ppb(i1,i2,i,ic)
                    endif
                  endif

!                  if(abs(ppb(i1,i2,i,ic))/=0d0  )then
!                     write (*,'("[",I4,"] ppbafp: ppb(",I3,",",I3,",",I3,",",I3,") = ",E11.3)')
!     &                    rank, i1,i2,i,ic,ppb(i1,i2,i,ic)
!                  endif

cccccccccccccccccccccccccccccccccccccccccccccccccc

c        pdb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(1,2))
c        dpb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(2,1))
c        ddb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(2,2))
ccccccccccccccccccc
c         xxx = ppb(5+1,1,i,1)
c         ppb(i1,i2,i,ic) = 0d0
c         ppb(5+1,1,i,1) = xxx
c         if(i1==6.and.i2==1.and.ic==1.and.ig==1) then
c          write(6,'(4i4,2d15.6,7i2)')lm,lmp,i,lmb, cgr(lm,lmp, lmb,ig),
c     &       ppbrd( l,n, lp,np, lb,nb, ioff(1,1))
c     &       ,l,n, lp,np, lb,nb, ioff(1,1)
c         endif
ccccccccccccccccccc
!      if(ppb(i1,i2,i,ic)/=0d0) then
!         write(6,'("[",I4,"] ppb(",4i4,") = ",d15.6)'),rank, i1,i2,i,ic,ppb(i1,i2,i,ic)
!      endif
ccccccccccccccccccc
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
cccccccccccccccccccccccccccccccccccccccccccc
c      pdb =0d0
c     dpb =0d0
c     ddb =0d0
ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc
c     print *,' ppbafp x sum ppb',sum(ppb)
c     stop
cccccccccccccccccccccccccccccccccccccccccccc

!ML      write(6,'("[",I4,"] ppbafp_v2: sum(ppb) = ",E11.3)') rank,sum(ppb)

      end

c--------------------------------------------------------------------
      subroutine rotcg(lmxax,symops,ng,cgr)
C- Rotate Clebsch-Gordan coefficients.
      implicit none
      integer(4) :: lmxax, ng, nlmxa,
     &              lnjcg, lnxcg,
     &              ilma,la,ilmb,lh,ii,indx,icg1,icg2,icg,
     &              ig,lm1,lm2,lm,l,m,lmxcg,ilm ,lx
      real(8) :: cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng),symops(9,ng)
      real(8),allocatable:: cg(:),dlmm(:,:,:,:),cgn(:,:,:)
      integer(4),allocatable :: jcg(:),indxcg(:)
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4):: ll
#endif

c ----------- MPI remark: -------------------------------------
c
c this routine assumes that ll is already set up.
c -> rotcg does not necessarily have to be called on all nodes.
c
c -------------------------------------------------------------

c --- CG coefficients. <LM3|lm1 lm2>
c inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
c Injcg = indxcg(inxcg) to indxcg(inxcg)-1
c cg(inxcg)  : = <lm3|lm1 lm2>
c jcg(lnjcg) : = lm3
c      print *,' rotcg:'
c      do ig=1,ng
c      write(6,*)' transposed symope ig  =',ig
c     write(6,'(3f12.6)') symops(1:3,ig)
c     write(6,'(3f12.6)') symops(4:6,ig)
c     write(6,'(3f12.6)') symops(7:9,ig)
c      enddo
      allocate(cgn((2*lmxax+1)**2,(lmxax+1)**2,(lmxax+1)**2) )
      cgn = 0d0
      nlmxa = (lmxax+1)**2
      lmxcg = lmxax
      if (lmxcg <= 6) then
        lnjcg = 6500
        lnxcg = 1300
      else if (lmxcg <= 8) then
        lnjcg = 22700
        lnxcg = 3400
      else if (lmxcg <= 10) then
        lnjcg = 62200
        lnxcg = 7400
      else
        write(6,*) 'rotcg: cannot handle lmxcg=',lmxcg
        call rx( 'rotcg: cannot handle lmxcg')
      endif
      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
      call scg(lmxcg,cg,indxcg,jcg)
c----
      do ilma = 1, nlmxa
        la = ll(ilma)
        do ilmb = 1, nlmxa
          lh = ll(ilmb)
          ii = max0(ilma,ilmb)
          indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do icg = icg1, icg2
            ilm  = jcg(icg)
            cgn(ilm, ilma,ilmb)  = cg(icg) ! ilm is move to 1st argument.!
          enddo
        enddo
      enddo
ccccccccccccccccccccccccccccccccc
c      print *,' cgn=',cgn(1,1,1),1/sqrt(4*3.1415926) !=1d0
c      stop
ccccccccccccccccccccccccccccccccc

C --- Rotation matrix ---
      allocate(dlmm( -2*lmxax:2*lmxax, -2*lmxax:2*lmxax, 0:2*lmxax,ng))
      call rotdlmm(symops, ng, 2*lmxax+1,dlmm)

c --- Rotate Clebsch Gordan coefficients ---
      do  lm =  1, (2*lmxax+1)**2
        l = ll(lm)
        m = lm - l**2 - l - 1
        lx = l**2 + l + 1
        do  ig  = 1, ng
          do  lm2 = 1, nlmxa
            do  lm1 = 1, nlmxa
              cgr(lm1,lm2,lm,ig) = sum(cgn(lx-l:lx+l,lm1,lm2)*dlmm(-l:l,m,l,ig))
C              sumr = 0d0
C              do md = -l,l
C                lmd = l**2 +l +1 + md
C                sumr = sumr + cgn(lm1,lm2,lmd)*dlmm(md,m,l,ig)
C              enddo
C              cgr(lm1,lm2,lm,ig) = sumr
            enddo
          enddo
        enddo
      enddo

C      do  ig = 1, ng
C      print *, ig; call prmx('cgr',cgr(:,:,:,ig),(lmxax+1)**4,(lmxax+1)**4,(2*lmxax+1)**2)
C      enddo

      deallocate( cg,dlmm,cgn,jcg,indxcg)
      end

      subroutine rotgvec(symops,nqibz,ngcmx,ngcn,qbas,ngveccB,
     o  ngveccBr)
C- Rotate ngveccB
C    G' = R(G), where R denotes rotation
Cr   Determine ngveccBr so that
Cr   matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq))) = matmul(qbas, ngveccBr(1:3, igc, iq)).
      implicit none
      integer :: nqibz,ngcmx,ngcn(nqibz),ngveccB(3,ngcmx,nqibz),ngveccBr(3,ngcmx,nqibz),iq,igc
      real(8) :: symops(3,3),qbas(3,3),qbasinv(3,3),rotnvg(3,3),vec(3),det,sumchk
      logical:: debug=.false.

C     call dinv33(qbas,0,qbasinv,det)
      call minv33(qbas,qbasinv)
      if(debug) write(6,'(3f12.6)') symops(1,1:3)
      if(debug) write(6,'(3f12.6)') symops(2,1:3)
      if(debug) write(6,'(3f12.6)') symops(3,1:3)

      rotnvg = matmul(qbasinv,matmul(symops,qbas))
      sumchk  =0d0
      do  iq = 1, nqibz
        if(debug) print *,' iq=',iq, ' sum ngveccB=', sum(abs(ngveccB(1:3,1:ngcn(iq), iq)))
        do igc = 1, ngcn(iq)
          vec  = matmul( rotnvg, ngveccB(1:3,igc, iq))
! vec should be almost integer and ngveccBr = vec.
! This expression rounds to the correct integer value.
          ngveccBr(1:3, igc, iq) = idint( vec + dsign(.5d0,vec))
        enddo
        do igc= 1, ngcn(iq)
          sumchk = sumchk +
     &    sum(abs(
     &       matmul(qbas, ngveccBr(1:3, igc, iq))
     &     - matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))))
        enddo

C        if(debug) print *," rotgvec: nmin nmax="
C     & ,minval(ngveccBr(1:3, 1:ngcn(iq), iq))
C     & ,maxval(ngveccBr(1:3, 1:ngcn(iq), iq))

      enddo
      if(abs(sumchk)/nqibz/minval(ngcn)>1d-4) then
        print *," rotgvec: sum chk error sumchk=",sumchk
        call rx( "rotgvec: sum chk error >1d-4")
      endif
      end

      subroutine q0iwgt(symops,ngrp,q0i,nq0i,
     o                  wgt0)
c Get weight for each k-points near 0.
c wgt0(irreducible-k, irotation)
      implicit none
      integer(4) :: ixx,i, ngrp,ig,nq0i
      real(8)     :: q0i(3,6),symops(3,3,ngrp),sym(3,3),
     &   qt(3), q0in(3,ngrp*nq0i), wgt0(nq0i,ngrp)
c      print *,' q0iwgt:'
      wgt0 = 0d0
      ixx = 0
      do i = 1,nq0i
        qt = q0i(:,i)
c equivalence check
        do ig = 1,ngrp
          sym = symops(:,:,ig)
c          write(66,'(2i3,3f12.5)')i,ig,matmul(sym,qt)
ccccccccccccccccccccccccccccccccccccccc
c          do ix = 1,ixx
c            if(sum(abs(q0in(:,ix)-matmul(sym,qt)))<1d-10) then
c              goto 1111
c            endif
c          enddo
ccccccccccccccccccccccccccccccccccccccc
          ixx = ixx+1
          q0in (:,ixx) = matmul(sym,qt)
          wgt0 (i,ig)  = 1d0
c          write(6,'(2i3,3f12.5)')i,ig,q0in (:,ixx)
 1111     continue
        enddo
      enddo
      wgt0  = wgt0 /ixx
ccccccccccccccccccccc
c     do ig=1,ngrp
c     write(6,*)' sym ig  =',ig
c     write(6,'(3f12.6)') symops(1,1:3,ig)
c     write(6,'(3f12.6)') symops(2,1:3,ig)
c     write(6,'(3f12.6)') symops(3,1:3,ig)
c      enddo
ccccccccccccccccccccc
      end

c--------------------------------------------------
      subroutine melpln
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc,ngvecc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied
     i       shtv,q, qbas,
     o       zmelp)
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
co rmelp, cmelp
c-----------------------------------------------------
      implicit none
      integer(4):: ngp1, ngvecp1(3,ngp1), ngp2, ngvecp2(3,ngp2),
     &     itp,it, ngc,ngvecc(3,ngc),ig1,ig2,igc,ntp0,nt0,
     &     nadd(3)
      complex(8):: geigq1(ngp1,ntp0),
     &             geigq2(ngp2,nt0), zmelp(ngc,nt0,ntp0),phase(ngc),
     &             img=(0d0,1d0)
      real(8)::  q(3),qg(3),qbas(3,3),shtv(3)
     & ,pi=3.1415926535897932D0
c      inum =0
c
c      print *,' melpln: sht nadd=',shtv,nadd
      do igc = 1,ngc
        qg(1:3) = 2*pi * (q(1:3)+ matmul(qbas, ngvecc(1:3,igc)))
        phase(igc) = exp(img*sum(qg*shtv))
      enddo

      zmelp=0d0
      do ig1=1,ngp1
        do ig2=1,ngp2
          do igc=1,ngc
            if(sum(abs( ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2)
     &             -ngvecc(1:3,igc) -nadd))  ==0) then
c          inum = inum+1
! igdelta(3,inum) means the delta function delta_G1,G2+G.
c          igdelta(1:3,inum) = (/ig1,ig2,igc/)
              do itp=1,ntp0
                do it =1,nt0
                  zmelp(igc,it,itp) = zmelp(igc,it,itp)
     &        + dconjg(geigq1(ig1,itp))*geigq2(ig2,it) * phase(igc)

cccccccccccccccccccccccccccccccccccccccccccc
c      if(sum(abs( ngvecp1(1:3,ig1)))<3.and.
c     &   igc==4.and.it==32.and.itp==1) then
c         write(6,'("  g1 g2 zmeladd=",3i3,1x,3i3,1x,3i3,12d12.3)')
c     &   ngvecp1(1:3,ig1),ngvecp2(1:3,ig2),ngvecc(1:3,igc),
c     &   dconjg(geigq1(ig1,itp))*geigq2(ig2,it) * phase(igc),
c     &    geigq1(ig1,itp),geigq2(ig2,it) ,phase(igc)
c      endif
cccccccccccccccccccccccccccccccccccccccccc
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
c      nnum=inum
c      do inum=1,nnum
c        (/ig1,ig2,igc/) = igdelta(1:3,inum)
c      enddo
ccccccccccccccccccccccccccccccccccccc
c      igc=4
c      it = 32
c      itp=1
c      write(6,'("  igc it itp zmel=",3i4,2d15.6)')
c     &  igc,it,itp,zmelp(igc,it,itp)
c      stop 'TEST end melpln:'
ccccccccccccccccccccccccccccccccccccc
      end

      subroutine melpln4(ngp1,ngvecp1,ngp2,ngvecp2,nbloch,ngc,nadd,geigq1,ntp0,geigq2,
     .  nt0,nctot,ncc,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,shftv,q1,q2,symope,qlat,zpsimb)
C- Interstitial plane wave contributions to matrix elements <M psi |psi>
C ----------------------------------------------------------------------
Ci Inputs
Ci   ngp1   : number of G-vectors representing eigenfunction at q1
Ci          : NEW: ngq1
Ci   ngvecp1: The G-vectors at q1 (as integer multiples of qlat)
Ci          : NEW: igv1
Ci   ngp2   : number of G-vectors representing eigenfunction at q2
Ci          : NEW: ngq2
Ci   ngvecp2: The ngq2 G-vectors at q2 (integer multiples of qlat)
Ci          : NEW: igv2
Ci   nbloch : total number of local product basis functions within MTs
Ci          : NEW: nrprd
Ci   ngc    : number of PWs in interstitial product basis
Ci          : NEW: ngprd
Ci   nadd   : G-vector shift, normally 0 (integer multiples of qlat)
Ci   geigq1 : Coefficients to PW representation eigenvectors associated with q1
Ci   ntp0   : Polarizability context: no. unoccupied states at q+k
Ci          : Self-energy context:    no. states at q for which sigm calculated (n or m in Eq 34)
Ci          : NEW: n1
Ci   geigq2 : Coefficients to PW representation eigenvectors associated with q2
Ci   nt0    : Polarizability context: no. occupied states at k, n in Eq. 32, PRB76, 165106
Ci          : Self-energy context:    no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Ci          : NEW: n2
Ci   nctot  : Offset in 2nd index zpsimb array: first element is zpsimb(:,off2+1,:)
Ci          : Typically off2 = number of core states.
Ci          : NEW off2
Ci   ncc    : Offset in 3rd index zpsimb array: first element is zpsimb(:,:,off3+1)
Ci          : Typically off3 = number of core states when there is no time reversal symmetry.
Ci          : NEW off3
Ci   ndpsi1 : parameters dimensioning zpsimb
Ci   ndpsi2 :
Ci   ndpsi3 :
Ci   shftv  : shift entering into phase for q1+G
Ci   q1     : Polarizability context: q+k, calculating <psi(q+k,t') PW | psi(k,t) >
Ci          : used in integration over k to make Pi(q). Eq. 32, PRB76, 165106
Ci          : Self-energy context: q at which sigma is calculated. Eq. 34, PRB76, 165106
Ci   q2     : Polarizability context: k, for ME <psi(q+k,t') PW | psi(k,t) > .. Eq. 32, PRB76, 165106
Ci          : Self-energy context: k for integration dk psi*(q) psi(q-k)    ... Eq. 34, PRB76, 165106
Ci          : q2 = q - R*k (R is rotation)
Ci   symope : Rotation R
Ci   qlat   : reciprocal lattice vectors
Ci          : NEW: qlat
Co Outputs
Co   zpsimb :The IPW part of matrix elements with the product basis, <psi |psi Plane>
Cr Remarks
Cr   use PPOVL through getppx in m_read_ppovl
Cu Updates
Cu   07 May 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngp1,ngp2,nbloch,ngc,ntp0,nctot,ncc,nt0,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi
      integer ngvecp1(3,ngp1),ngvecp2(3,ngp2),nadd(3)
      real(8):: shftv(3),q1(3),q2(3),qlat(3,3),symope(3,3)
C     complex(8):: geigq1(ngp1,ntp0),geigq2(ngp2,nt0),zpsimb(nbloch+ngc,nctot+nt0,ncc+ntp0)
      complex(8):: geigq1(ngp1,ntp0),geigq2(ngp2,nt0),zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
C ... Local parameters
      integer ig,itp
      complex(8),allocatable :: zwk(:,:)

C     call zprm('geigq2',geigq2,ngp2,ngp2,nt0)

      if (ipblo /= 1 .or. it1lo /= 1 .or. it2lo /= 1) call rx('fix melpln4')
      call melpln3(ngp1,ngvecp1,ngp2,ngvecp2,nbloch,ngc,nadd,geigq1,ntp0,geigq2,nt0,
     .  ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,shftv,q1,q2,symope,qlat,zpsimb(nbloch+1,1,ncc+1))

C ... Undo 2nd column misalignment: shift columns 1:nt in zpsimb to 1+nctot:nt:nctot
      if (nctot > 0) then
        allocate(zwk(nt0,ntp0))
        do  ig = 1, ngc
          call zcopy(nt0*ntp0,zpsimb(nbloch+ig,1,ncc+1),ipbhi-ipblo+1,zwk,1)
          do  itp = 1, ntp0
            zpsimb(nbloch+ig, 1:nctot, ncc+itp) = 0d0
            zpsimb(nbloch+ig, nctot+1:nctot+nt0, ncc+itp) = zwk(1:nt0,itp)
          enddo
        enddo
        deallocate(zwk)
      endif

      end

      subroutine melpln3(ngp1,ngvecp1,ngp2,ngvecp2,nbloch,ngc,nadd,geigq1,ntp0,geigq2,nt0,
     .  ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,shftv,q1,q2,symope,qlat,zpsimb)
C- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ngp1   : number of G-vectors representing eigenfunction at q1
Ci          : NEW: ngq1
Ci   ngvecp1: The G-vectors at q1 (as integer multiples of qlat)
Ci          : NEW: igv1
Ci   ngp2   : number of G-vectors representing eigenfunction at q2
Ci          : NEW: ngq2
Ci   ngvecp2: The ngq2 G-vectors at q2 (integer multiples of qlat)
Ci          : NEW: igv2
Ci   nbloch :Not used REMOVE
Ci   ngc    : number of PWs in interstitial product basis
Ci          : NEW: ngprd
Ci   nadd   : G-vector shift, normally 0 (integer multiples of qlat)
Ci   geigq1 : Coefficients to PW representation eigenvectors associated with q1
Ci   ntp0   : (polarizability context) no. unoccupied states at q+k (polarizability context)
Ci          : (self-energy context) no. states at q for which sigm calculated (n or m in Eq 34)
Ci          : NEW: n1
Ci   geigq2 : Coefficients to PW representation eigenvectors associated with q2
Ci   nt0    : (polarizability context) no. occupied states at k, n in Eq. 32, PRB76, 165106
Ci          : (self-energy context) no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Ci          : NEW: n2
Ci   ndpsi1 : parameters dimensioning zpsimb
Ci   ndpsi2 :
Ci   ndpsi3 :
Ci   shftv  : shift entering into phase for q1+G
Ci   q1     : Polarizability context: q+k, calculating <psi(q+k,t') PW | psi(k,t) >
Ci          : used in integration over k to make Pi(q). Eq. 32, PRB76, 165106
Ci          : Self-energy context: q at which sigma is calculated. Eq. 34, PRB76, 165106
Ci   q2     : Polarizability context: k, for ME <psi(q+k,t') PW | psi(k,t) > .. Eq. 32, PRB76, 165106
Ci          : Self-energy context: k for integration dk psi*(q) psi(q-k)    ... Eq. 34, PRB76, 165106
Ci          : q2 = q - R*k (R is rotation)
Ci   symope : Rotation R
Ci   qlat   :reciprocal lattice vectors
Ci          :NEW: qlat
Ci   ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi : dimensions of zpsimb
Co Outputs
Co   zpsimb :The IPW part of matrix elements with the product basis, <psi |psi Plane>
Cr Remarks
Cr   Reads from PPOVL, IPW matrix elements ppx = <Gc|Gc>^-1 <Gc|Gc2> (Gc2 spans G1-G2)
Cu Updates
Cu   17 May 13 Adapted from melpln2.  Requires Sep12 version of rppovl.F
C ----------------------------------------------------------------------
      use m_read_ppovl,only:getppx, ppx,ngc2,ngvecc2
      implicit none
C ... Passed parameters
      integer nbloch,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,ngc,ngp1,ngp2,nt0,ntp0
      integer nadd(3),ngvecp1(3,ngp1),ngvecp2(3,ngp2)
      real(8):: q1(3),q2(3),qlat(3,3),shftv(3),symope(3,3)
c     complex(8):: geigq1(ngp1,ntp0),geigq2(ngp2,nt0),zpsimb(nbloch+ngc,nt0*ntp0)
      complex(8):: geigq1(ngp1,ntp0),geigq2(ngp2,nt0),zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
C ... Local parameters
      integer itp,it,ig1,ig2,igc,igc2,ngcx,igmn,igmx,nblk,nblki,
     .  nxmin,nymin,nzmin,nxmax,nymax,nzmax,verbose,idalloc,iverbose
      integer omppid,nthreads,tid,iblk1,iblk2
      integer,allocatable:: firstb(:),lastb(:),ipp(:,:,:)
      integer nn(1:3)
      real(8):: qg(3),pi=3.1415926535897932D0
      complex(8),allocatable :: gg(:,:,:),zwk(:)
C     complex(8),allocatable :: ggl(:,:,:)
      complex(8):: img=(0d0,1d0),phas,zfac

C ... Setup
      if (nt0 == 0) return
      iverbose = (verbose()/50)*10
C     getppx returns ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2)
C     ppx = <Gc|Gc>^-1 <Gc|Gc2> (Gc2 spans G1-G2) made in in rdata4gw_v2 (Sep12 version)
      call getppx(qlat,ngcx,q2,symope)
      if (ngc/=ngcx) call rx('melpln3: file mismatch with ppovl, no G vectors')
      nxmin = minval(ngvecc2(1,1:ngc2)); nxmax = maxval(ngvecc2(1,1:ngc2))
      nymin = minval(ngvecc2(2,1:ngc2)); nymax = maxval(ngvecc2(2,1:ngc2))
      nzmin = minval(ngvecc2(3,1:ngc2)); nzmax = maxval(ngvecc2(3,1:ngc2))
C     Sanity check
      do  ig1 = 1, ngp1
      do  ig2 = 1, ngp2
        nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
        if (nn(1)<nxmin .or. nn(1)>nxmax) stop "melpln3:nn(1) out of range"
        if (nn(2)<nymin .or. nn(2)>nymax) stop "melpln3:nn(2) out of range"
        if (nn(3)<nzmin .or. nn(3)>nzmax) stop "melpln3:nn(3) out of range"
      enddo
      enddo
      allocate(ipp(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
      ipp = -10000000
      do  igc2 = 1, ngc2
        ipp(ngvecc2(1,igc2),ngvecc2(2,igc2),ngvecc2(3,igc2)) = igc2
      enddo

C --- Assemble gg in slices to conserve memory.  Slice ranges from igmn to igmx ---
C     nblk is max number of G2 vectors.  Should be adjusted depending on total size of gg
      nblk = 200
C$OMP parallel
      nthreads = omppid(0)
C$OMP end parallel
      nthreads = max(nthreads,1)
      if (nthreads > 1) nblk = nblk - mod(nblk+nthreads,nthreads)
      igc2 = idalloc('gg',iverbose+2,nt0*ntp0,nblk*2)
      allocate(gg(nt0,ntp0,nblk)) ! order gg indices to be compatible with zpsimb

C     print *, 'starting time',nblk,ngc2,ngp2,delwc()

      igmx = 0
      do  igmn = 1, ngc2, nblk
        igmx = min(igmx+nblk,ngc2)
        nblki = igmx-igmn+1  ! Size of this block
        call dinit(gg,2*nt0*ntp0*nblk)

C#ifdef OPENMP
        allocate(firstb(0:nthreads-1),lastb(0:nthreads-1))
        call thdist1(0,0d0,0d0,nthreads,1,nblki,firstb,lastb)
C$OMP   parallel private(tid,nn,ig1,ig2,igc2,igc,qg,phas,it,zwk,itp,iblk1,iblk2,zfac) shared(gg)
        allocate(zwk(nt0))
        tid = omppid(1)

        iblk1=firstb(tid); iblk2=lastb(tid)
C   ... Loop over G1
        do  ig1 = 1, ngp1
C     ... Loop over G2
          do  ig2 = 1, ngp2
            nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
            igc2    = ipp(nn(1),nn(2),nn(3))
            igc     = igc2-igmn+1
            if (igc < iblk1 .or. igc > iblk2) cycle
            qg(1:3) = 2d0*pi * (q1(1:3) + matmul(qlat, nn(1:3)))
            phas = exp(-img*sum(qg*shftv))
            do  it = 1, nt0
              zwk(it) = dconjg(geigq2(ig2,it))*phas
            enddo
C           Loop over states at q1
            do  itp = 1, ntp0
              zfac = geigq1(ig1,itp)
C             Loop over states at q2
              do  it = 1, nt0
                gg(it,itp,igc) = gg(it,itp,igc) +  zfac*zwk(it)
              enddo
            enddo
          enddo
        enddo
        deallocate(zwk)
C$OMP   end parallel
        deallocate(firstb,lastb)

C#elseC
CC$OMP   parallel private(ig1,ig2,nn,igc2,igc,qg,phas,it,zwk,itp,zfac) shared(gg)
C        allocate(zwk(nt0))
CC   ... Loop over G1
C        do  ig1 = 1, ngp1
CC     ... Loop over G2
CC$OMP     do
C          do  ig2 = 1, ngp2
C            nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
C            igc2    = ipp(nn(1),nn(2),nn(3))
C            igc     = igc2-igmn+1
C            if (igc < 1 .or. igc > nblki) cycle
C            qg(1:3) = 2d0*pi * (q1(1:3) + matmul(qlat, nn(1:3)))
C            phas = exp(-img*sum(qg*shftv))
C            do  it = 1, nt0
C              zwk(it) = dconjg(geigq2(ig2,it))*phas
C            enddo
CC           Loop over states at q1
C            do  itp = 1, ntp0
C              zfac = geigq1(ig1,itp)
CC             Loop over states at q2
C              do  it = 1, nt0
C                gg(it,itp,igc) = gg(it,itp,igc) +  zfac*zwk(it)
C              enddo
C            enddo
C          enddo
CC$OMP     end do
C        enddo
C        deallocate(zwk)
CC$OMP   end parallel
C#endif
        zfac = (0d0,0d0)
        if (igmn > 1) zfac = (1d0,0d0)
        call zgemm('N','T',ngc,ntp0*nt0,nblki,(1d0,0d0),ppx(1,igmn),ngc,
     .    gg,ntp0*nt0,zfac,zpsimb,ipbhi-ipblo+1)
      enddo
C     call zprm('zpsimb',zpsimb,ipbhi-ipblo+1,ngc,nt0*ntp0)

      igc2 = idalloc('gg',iverbose+4,nt0*ntp0,nblk*2)
      deallocate(ppx,gg,ngvecc2,ipp)

C     <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
      call ztoccz(zpsimb,ipbhi-ipblo+1,ngc,nt0*ntp0)

      end

C      subroutine melpln3t
C     .  (ngp1,ngvecp1,ngp2,ngvecp2,nbloch,ngc,nadd,
C     .   geigq1, ntp0,           ! q1=q    ---> iq 1:ntp0 q-point (unoccupied)
C     .   geigq2, nt0,            ! q2=q-rk ---> kp 1:nt0 (occupied)
C     .   shftv,q,qi,symope,qbas,
C     .   zpsimb)
CC- Contributions to matrix elements <Plan psi |psi> from interstitial plane wave.
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ngp1   :number of G-vectors representing eigenfunction at q1
CCi   ngvecp1:vector of ngp1 G-vectors at q1 as multiples of qlat
CCi   ngp2   :number of G-vectors representing eigenfunction at q2
CCi   ngvecp2:vector of ngp2 G-vectors at q2 as multiples of qlat
CCi   nbloch :total number of Bloch product basis functions, needed
CCi          :for dimensioning and offsets
CCi   ngc    :number of IPWs in product basis
CCi   nadd   :G-vector shift (why?)
CCi   geigq1 :PW coefficients for eigenvectors associated with q1
CCi   ntp0   :number of eigenvectors at q1 (unoccupied bands)
CCi   geigq2 :PW coefficients for eigenvectors associated with q2
CCi   nt0    :number of eigenvectors at q2 (occupied bands)
CCi   shftv   :shift entering into phase for q1+G = g(delta_{g^-1}) (?)
CCi   q1     :q
CCi   q2     :q-k
CCi   symope :Unit symmetry operation
CCi   qbas   :reciprocal lattice vectors
CCo Outputs
CCo   zpsimb :The IPW part of matrix elments with the product basis
CCl Local variables
CCl  ngcgp2 :number of G vectors Gc+G2
CCl  nggg   :number of G vectors Gc+G2-G1  (module, read from file PPOVLGG)
CCl  ggg    :<Gc G2|G1>, equivalent to <-Gc+G1-G2> (module, read from file PPOVLGG)
CCl  nvggg  :vector of Gc+G2-G1 as multiples of qlat (module, read from file PPOVLGG)
CCl  ngcgp  :number of Gc (?) (module, read from file PPOVLGG)
CCl  ngcgp  :number of G vectors spanning IPWs + basis (?) (module, read from file PPOVLGG)
CCl  nvgcgp2:vector of symope(Gc)+G2 as multiples of qlat (?) (module, read from file PPOVLGG)
CCl  ngcread:number of IPWs in product basis as in: (module, read from file PPOVLGG)
CCl ppovlinv:inverse of <Gc|Gc> overlap matrix (module, read from file PPOVLG.iq)
CCr Remarks
CCr   use PPOVL through getppx in m_read_ppovl
CCu Updates
CCu   17 May 13 Adapted from melpln2
CC ----------------------------------------------------------------------
CC     use m_read_ppovl,only:getppx, ppx,ngc2,ngvecc2
C      use m_read_ppovl,only: getppx2,ggg,ppovlinv,
C     .  nvggg,nvgcgp2,ngvecc,nggg,ngcgp,ngcread
C      implicit none
CC ... Passed parameters
C      integer ngp1,ngp2,nbloch,ngc,ntp0,nt0,nadd(3),
C     .        ngvecp1(3,ngp1),ngvecp2(3,ngp2)
C      real(8):: q(3),qi(3),qbas(3,3),shftv(3),symope(3,3)
C      complex(8):: geigq1(ngp1,ntp0),geigq2(ngp2,nt0)
C      complex(8):: zpsimb(nbloch+ngc,nt0*ntp0)
CC ... Dynamically allocated local arrays
C      complex(8),allocatable::ggitp(:,:),gp2phas2(:),phase(:)
C      integer,allocatable:: firstb(:),lastb(:)
C      integer,allocatable:: igggi(:,:,:),igcgp2i(:,:,:),ngveccR(:,:)
C      complex(8),allocatable :: gg(:,:,:),zwk(:)
CC     complex(8),allocatable :: ggl(:,:,:)
C      complex(8),allocatable:: mmm(:,:),zmelp0(:,:,:),ggitp_(:,:)
CC ... Local parameters
C      integer itp,it,ig1,ig2,igc,igc2,ngcx,igmn,igmx,nblk,nblki,iggg,igp1,igp2,igcgp2,
C     .  nxmin,nymin,nzmin,nxmax,nymax,nzmax,verbose,idalloc,iverbose
C      integer ngcs(1),nxminggg,nyminggg,nzminggg,nxmaxggg,nymaxggg,nzmaxggg
C      integer ngcgp2,nxmingcgp2,nymingcgp2,nzmingcgp2,nxmaxgcgp2,nymaxgcgp2,nzmaxgcgp2
C      integer omppid,nthreads,tid,iblk1,iblk2
C      integer nn(1:3)
C      integer,allocatable :: ipp(:,:,:)
C      real(8):: qg(3),pi=3.1415926535897932D0
C      real(8):: delwc
C      complex(8):: img=(0d0,1d0),phas,zfac
C
CC ... Setup
C      if (nt0 == 0) return
C      iverbose = (verbose()/50)*10
C      if (verbose()>=90) print *,' melpln3t: '
CC     call getppx(qbas,  ngcx,qi,symope) ! rppovl.f
CC     Read from file PPOVLG : nggg,ngcread,ngcgp,ggg(nggg),nvggg(1:nggg),nvgcgp2(3,ngcgp2)
C      call getppx2(qbas,qi) ! rppovl.f
C      if (ngc/=ngcread) then
C         print *,qi,ngcread,ngc
C         print *, 'qi',qi, 'ngcread',ngcread, 'ngc',ngc
C         call rx( 'melpln3t: ngc different from ngcread read from PPOVLG by getppx')
C      endif
C      ngcgp2=ngcgp
C
CC ... Make igggi ... maps G vectors nvggg to multiples of qlat
C      nxminggg =minval(nvggg(1,1:nggg)) ;nxmaxggg=maxval(nvggg(1,1:nggg))
C      nyminggg =minval(nvggg(2,1:nggg)) ;nymaxggg=maxval(nvggg(2,1:nggg))
C      nzminggg =minval(nvggg(3,1:nggg)) ;nzmaxggg=maxval(nvggg(3,1:nggg))
C      allocate( igggi(nxminggg:nxmaxggg,nyminggg:nymaxggg,nzminggg:nzmaxggg) )
C      igggi = -100000
C      do  iggg = 1, nggg
C        igggi(nvggg(1,iggg),nvggg(2,iggg),nvggg(3,iggg)) = iggg
C      enddo
C
CC ... Make ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shftv)*exp(-i(q-Gadd)*shftv)
C      allocate(ngveccR(1:3,1:ngc))
C      ngcs(1) = ngc
C      call rotgvec(symope, 1, ngc, ngcs, qbas, ngvecc, ngveccR)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C!! simple case commented out. This worked (but slow) nov2013takao.
Cc$$$!!      print *,'--- goto simple case for debug---'
Cc$$$      print *,'fffffff: ngc ngp1 ngp2=',ngc,ngp1,ngp2
Cc$$$      zmelp0 = 0d0
Cc$$$      do igc=1,ngc
Cc$$$        phas =exp( img*2d0*pi*sum((q+matmul(qbas,ngveccR(:,igc)))*shftv) ) !conjg
Cc$$$      do igp1=1,ngp1     !for ngp1
Cc$$$      do igp2=1,ngp2
Cc$$$c        nn = ngvecp1(:,igp1)-ngvecp2(:,igp2)-ngveccR(:,igc) - nadd
Cc$$$        nn = -ngvecp1(:,igp1) + ngvecp2(:,igp2) + ngveccR(:,igc) + nadd
Cc$$$        iggg = igggi(nn(1),nn(2),nn(3))
Cc$$$        do it=1,nt0
Cc$$$        do itp=1,ntp0
Cc$$$          zmelp0(igc,it,itp) = zmelp0(igc,it,itp)
Cc$$$     &      +  phas*ggg(iggg)* dconjg(geigq1(igp1,itp))* geigq2(igp2,it)
Cc$$$        enddo
Cc$$$        enddo
Cc$$$      enddo
Cc$$$      enddo
Cc$$$      enddo
Cc$$$      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
Cc$$$      return
Cccccccccccccccccccccccccccccccccccccccccccccccccc
C
CC ... Accumulate for each unocc state ggitp(1:ntp0,igcg2) : ggg(iggg)*geigq1(igp1,1:ntp0)
CC     where ggg = <Gc G2|G1> = <qt+G1 -(rk+Gc) -(q_rk+G2)>
CC     Skips over G vectors for which ggg is is not available
C      allocate(ggitp(ntp0,ngcgp2)) ! element for each (unocc state,Gc+G2)
C      ggitp = 0d0
C      do igp1   = 1, ngp1   !for ngp1
C      do igcgp2 = 1, ngcgp2 !for ngc+ngp2 : a large number but most are missing
C         nn = ngvecp1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd  where -Gadd= -rk + qt -q_rk
C         if(nn(1)<nxminggg .or. nxmaxggg<nn(1) ) cycle
C         if(nn(2)<nyminggg .or. nymaxggg<nn(2) ) cycle
C         if(nn(3)<nzminggg .or. nzmaxggg<nn(3) ) cycle
C         iggg = igggi(nn(1),nn(2),nn(3))
C         if(iggg<0) cycle
CC        Accumulate into each unocc state ggitp(1:ntp0,igcgp2) += ggg(iggg)*geigq1+(1:ntp0,igp1)
C         call zaxpy(ntp0,ggg(iggg), geigq1(igp1,1),ngp1, ggitp(1,igcgp2),1)
C      enddo
C      enddo
C      deallocate(igggi)
C      ggitp = dconjg(ggitp)     ! Is this necessary?
C
CC ... Make igcgp2i ... maps G vectors nvgcgp2 to multiples of qlat
C      nxmingcgp2 =minval(nvgcgp2(1,1:ngcgp2)) ;nxmaxgcgp2=maxval(nvgcgp2(1,1:ngcgp2))
C      nymingcgp2 =minval(nvgcgp2(2,1:ngcgp2)) ;nymaxgcgp2=maxval(nvgcgp2(2,1:ngcgp2))
C      nzmingcgp2 =minval(nvgcgp2(3,1:ngcgp2)) ;nzmaxgcgp2=maxval(nvgcgp2(3,1:ngcgp2))
C      allocate(igcgp2i(nxmingcgp2:nxmaxgcgp2,nymingcgp2:nymaxgcgp2,nzmingcgp2:nzmaxgcgp2))
C      igcgp2i = -100000
C      do  igcgp2 = 1, ngcgp2
C        igcgp2i(nvgcgp2(1,igcgp2),nvgcgp2(2,igcgp2),nvgcgp2(3,igcgp2)) = igcgp2
C      enddo
C
CC ... zmelp0 <=  \sum_G2 ggitp(Gc+G2) geigqg2(G2))
CC     zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
CC     ? note \bfr'= g (\bfr) +\delta_g  (\bfr= {\bf r})
CC     ? mapping of function g[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
C      allocate(gp2phas2(nt0),phase(ngc))
C      do  igc = 1, ngc
C        phase(igc) = exp(img*2d0*pi*sum((q+matmul(qbas,ngveccR(:,igc)))*shftv))
C      enddo
C      allocate(ggitp_(ngp2,ngcgp2))
C      allocate(zmelp0(ngc,nt0,ntp0))  ! element for each IPW, occ, unocc state
C      zmelp0 = 0d0
C      do  igc = 1, ngc
C        do  itp = 1, ntp0
C          do  igp2 = 1, ngp2  ! For each G component of occupied state
C            nn = ngveccR(:,igc) + ngvecp2(:,igp2)  ! Gc+G2
C            igcgp2 = igcgp2i(nn(1),nn(2),nn(3))    ! corresponding G vector
C            ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
C          enddo
C!     ! This do block is time-consuming when we use large QpGcut_psi _cou.
C!     ! But this version (nov2013) is not so bad, maybe.
C          do  it = 1, nt0
C            zmelp0(igc,it,itp) = zmelp0(igc,it,itp) + sum( geigq2(:,it)*ggitp_(:,igc) )
C          enddo
C        enddo
C      enddo
C
C      call rx('melpln3t not finished')
C
C      end

      subroutine ztoccz(zmat,ndz,n1,n2)
C- Turn a subblock of a matrix into its complex conjugate
      implicit none
      integer ndz,n1,n2
      complex(8):: zmat(ndz,n2)
      integer i1,i2
!$omp parallel do collapse(2)
      do  i2 = 1, n2
        do i1 = 1, n1
          zmat(i1,i2) = dconjg(zmat(i1,i2))
        enddo
      enddo
!$omp end parallel do
      end
ccccccccccccccccccccccccccccccccc
c      write(6, '(" *** q=",(3f10.4,2x))'),q(1:3)
c      do ig1 = 1,ngp1
c        write(6, '(" ig geig =",i2,2x,3i3,4(2f10.4,2x))')
c     &   ig1,ngvecp1(1:3,ig1), ( geigq1(ig1,itp),itp=1,4)
c      enddo
ccccccccccccccccccccccccccccccccc
cccccccccccccccccccc
c      write(6,"( ' --- sum check gg=',10d13.6)")
c     &  sum(gg(1:ngc2,1,1)),sum(abs(gg(1:ngc2,1,1)))
cccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
c      allocate(xxx(ngc,ngc2))
c      call matm(ppovlinv,ppovl,xxx,ngc,ngc,ngc2)
c     print *,'invtest'
c      ix=1
c     do iy=1,ngc2
c       write(6,'("  ix iy xxx=",2i6,2f12.6)')  ix,iy,xxx(ix,iy)
c      enddo
c     stop ' invtest end'
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccc
c      do itp=1,ntp0
c      do it =1,nt0
c      do igc = 1,ngc
c      do igc2= 1,ngc2
c        if ( sum(abs(ngvecc(:,igc)-ngvecc2(:,igc2)))==0 ) then
c          zmelp(igc, it,itp) = gg(igc2,it,itp)
c        endif
c      enddo
c      enddo
c      enddo
c      enddo
c      zmelp = dconjg(zmelp) ! <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
c      deallocate(ppx,gg,ngvecc2)
c      return
ccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccc
c      write(6, '(" *** q=",(3f10.4,2x))'),q(1:3)
c      do it = 1,nt0
c      do igc = 1,ngc
c        write(6, '(" it igc nvec gg =",2i4,2x,3i3,4(2f10.4,2x))')
c     &   it ,igc, ngvecc(1:3,igc), ( zmelp(igc,it,itp),itp=1,4)
c      enddo
c      enddo
ccccccccccccccccccccccccccccccccc

