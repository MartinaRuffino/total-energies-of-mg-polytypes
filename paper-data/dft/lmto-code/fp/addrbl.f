      subroutine addrbl(s_site,s_spec,s_lat,s_ham,s_bz,s_optic,isp,nsp,
     .  nspc,q,ndham,ndimh,lmet,lrout,lwtkb,ldmatk,ltso,wtkb,lswtk,swtk,
     .  iq,lfrce,ldos,lekkl,k1,k2,k3,smpot,vconst,lcplxp,numq,qval,evec,
     .  evl,nevl,ef0,def,esmear,emin,emax,ndos,dos,smrho,sumqv,sumev,f,
     .  dmatk,tso)
C- Adds to the smooth and local output density and to eigval sum
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos qkkl qhkl qhhl eqkkl eqhkl eqhhl pikk sigkk
Ci                pihk sighk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  fsmbl fsmbpw rsibl rsibl1 rlocbl bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp ngcut kmxt rsma
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  fsmbl uspecb fsmbpw rsibl tbhsi rsibl1 rlocbl bstrux
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc alat qlat vol plat awald tol nkd nkq gmax ng kv
Ci                kv2
Co     Stored:     *
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:cg indxcg jcg cy qlv dlv igv gv kv igv2
Cio    Passed to:  fsmbl hhigbl phhigb hklbl gklbl fklbl hsmbl rsibl
Cio               sugvec rlocbl bstrux hxpbl ghibl hxpgbl ghigbl hklgbl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Ci Inputs
Ci   isp   :current spin channel.  Always 1 in the noncollinear case
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   q     :Bloch vector
Ci   ndham :leading dimension of evl
Ci   ndimh :dimension of hamiltonian
Ci   lmet  :Not used unless lwtkb is nonzero (see also subzi)
Ci         :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file, if they are needed
Ci         :3 always make two band passes; weights never needed a priori
Ci         :4 BZ integration with 3-point scheme
Ci   lwtkb :0 set of weights is not given; use 3-point interpolation
Ci         :1 or 2 given a set of weights
Ci         :-1 needs weights, but not yet available
Ci         :4 Special interacting hamiltonian mode with the following generalization
Ci         :  1. eigenvectors may not be eigenstates of a 1-body hamiltonian
Ci         :  2. eigenvalues are not given; properties depending on eigenvalues
Ci         :     are not generated.
Ci         :  3. wtkb is supplied: wtkb and eigenvectors are synchronized,
Ci         :     but they are not necessarily ordered by energy
Ci         :  wtkb and evecs are synchronized, but their order is not prescribed
Ci   lrout :0 no eigenvectors generated; exit after accumulating
Ci         :  eigenvalue weights and dos
Ci   ldmatk:0 no density matrix
Ci         :1 density matrix generated for this q, eigenfunction basis
Ci   wtkb  :integration weights, needed if lwtkb is 1 or 2
Ci   lswtk :<1 do nothing
Ci         :1 given a set of weights, make 'spin weights' swtk
Ci   iq    :index to current k-point
Ci   lfrce :if nonzero, accumulate contribution to force
Ci   ldos  :if nonzero, accumulate density-of-states
Ci   k1,k2,k3 dimensions of smpot,smrho
Ci   smpot :smooth potential on uniform mesh (mkpot.f), for forces
Ci   vconst:additional constant potential
Co   numq  :number of Fermi levels. Usu. 1, but numq=3 if lmet=4
Ci   qval  :total valence charge
Ci   evec  :eigenvectors
Ci   evl   :eigenvalues
Ci   nev
Ci   ef0   :estimate for fermi level
Ci   def   :When lmet=4, charge also accmulated for ef0+def and ef0-def
Ci   esmear:(sampling integration) gaussian broadening
Ci         :sign and integer part are extensions; see mkewgt.f
Ci   emin  :energy lower bound when adding to sampling dos (not used unless ldos=1)
Ci   emax  :energy upper bound when adding to sampling dos (not used unless ldos=1)
Ci   ndos  :number of energy mesh points
ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   lekkl :0 do not accumulate eqkkl; 1 do accumulate eqkkl
Co Outputs
Co   sumqv :integrated charge, resolved by spin
Co   sumev :sum of eigenvalues
Co   dos   :sampling density of states, if ldos=1
Co   smrho :smooth density on uniform mesh
Co   f     :eigenvalue contribution to forces
Co   dmatk:density matrix in MTO basis set
Ci   swtk  :'spin weights' to determine global magnetic moment, nspc=2
Ci         : swtk = diagonal part of  (z)^-1 sigmz z
Ci         : where sigmz is the z component of Paul spin matrices
Cl Local variables
Cl   napw  :number of PWs in APW part of basis
Cl   igapwl:PWs in units of reciprocal lattice vectors,
Cl         :possibly modified if q is shortened.
Cl   wtkp  :q-point weights from symmetry operations
Cr Remarks
Cu Updates
Cu   13 Jan 14 (Ben Kaube) start on optical matrix elements of envelope fns
Cu   09 Aug 13 New option to resolve SO energy by site
Cu   12 Apr 13 First cut at generating density matrix
Cu   08 Feb 13 Internally shortens q vector
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   09 Jun 07 Makes spin weights (noncollinear case)
Cu   02 Jan 06 sumqv resolved by spin
Cu   17 Jan 05 Extension of esmear to Fermi distribution
Cu   23 Dec 04 Extended to spin-coupled case
Cu   18 Nov 04 Sampling integration properly handles Fermi distribtion
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   23 Jan 01 Added lrout switch
Cu   17 Jun 00 spin polarized
Cu   22 May 00 Adapted from nfp add_densw.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nsp,nspc,iq,k1,k2,k3,ldos,lmet,lrout,lwtkb,lswtk,
     .  lfrce,ldmatk,ltso,ndham,ndimh,ndos,numq,lekkl,lcplxp,nevl
      double precision def,ef0,emax,emin,esmear,qval,vconst
      double precision q(3),dos(ndos,2,isp),
     .  wtkb(ndham*nspc,nsp/nspc,*),swtk(ndham,nsp,*),
     .  evl(ndham,nsp),sumev(2,3),sumqv(3,2),f(3,*)
      double precision dmatk(ndham*nspc),tso(*)
      double complex evec(ndimh,nspc,ndimh,nspc),smrho(k1,k2,k3,isp),
     .  smpot(k1,k2,k3,isp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_bz)::    s_bz
      type(str_optic):: s_optic
C ... Dynamically allocated arrays
      complex(8),allocatable:: evecc(:,:,:,:),work(:,:,:,:)
      real(8),allocatable:: qpgv(:,:),qpg2v(:),ylv(:,:)
      integer, pointer :: igapwl(:,:)
      real(8), pointer :: wtkp(:)
C ... Local parameters
      integer i,k,nevec,nbas,nglob,ngabc(3),n1,n2,n3,ndimhx,
     .  lmxax,lmxa,nlmax,nlmto,ig,napw,lfrce0
      double precision vavg,wgt,alat,qlat(3,3),tpiba,vol
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer ipiv(ndimh*2)
      double precision ewgt(ndimh*nspc,max(numq,3))
      double precision qs(3),dlength,tol
      parameter (tol=1d-8)

      if (lwtkb < 0) return
      call tcn('addrbl')
      nbas = nglob('nbas')
      napw = s_lat%napw
      nlmto = ndimh-napw
      ndimhx = ndimh*nspc
      ngabc = s_lat%nabc
      wtkp => s_bz%wtkp
      call dpzero(ewgt,ndimh*nspc*max(numq,3))

C     Debugging : rotate eigenvector
CC     call zprm('z',2,evec,ndimhx,ndimhx,ndimhx)
C      allocate(evecc(ndimh,nspc,ndimh,nspc),work(ndimh,nspc,ndimh,nspc))
C      call snot(0,ndimh,evecc)
CC     call zprm('u',2,evecc,ndimhx,ndimhx,ndimhx)
C      call zgemm('N','N',ndimhx,ndimhx,ndimhx,(1d0,0d0),evecc,ndimhx,
C     .  evec,ndimhx,(0d0,0d0),work,ndimhx)
C      call dcopy(ndimhx**2*2,work,1,evec,1)
CC     call zprm('z-rot',2,evec,ndimhx,ndimhx,ndimhx)
C      deallocate(evecc,work)

C ... Shorten q; shift APW G vectors to correspond
      igapwl => s_lat%igv2
      call shorps(1,s_lat%qlat,(/72,2,2/),q,qs)
      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        allocate(igapwl(3,napw))
        call shorigv(q,qs,s_lat%plat,napw,s_lat%igv2,igapwl)
      endif

C ... Setup for PW part of basis
      alat = s_lat%alat
      qlat = s_lat%qlat
C     Find largest lmxa ... should be made elsewhere
      lmxax = -1
      do  i = 1, nbas
        k = s_site(i)%spec
        lmxa = s_spec(k)%lmxa
        lmxax = max(lmxax,lmxa)
      enddo
      nlmax=(lmxax+1)**2
      if (napw > 0) then
        allocate(ylv(napw,nlmax),qpgv(3,napw),qpg2v(napw))
        tpiba = 2d0*4d0*datan(1d0)/alat
        do  ig = 1, napw
          qpgv(:,ig) = tpiba * ( qs + matmul(qlat,igapwl(:,ig)) )
        enddo
        call ropyln(napw,qpgv(1,1:napw),qpgv(2,1:napw),qpgv(3,1:napw),
     .    lmxax,napw,ylv,qpg2v)
      else
        allocate(ylv(1,1),qpgv(1,1),qpg2v(1))
      endif

C --- Decide how many states to include and make their weights ---
C ... Case band weights not passed: make sampling weights
      if (lwtkb == 0) then
        call rxx(nspc /= 1,'lwtkb=0 not implemented in noncoll case')
        wgt = abs(wtkp(iq))/nsp
        call mkewgt(lmet,wgt,qval,ndimh,evl(1,isp),ef0,def,esmear,numq,
     .    nevec,ewgt,sumev,sumqv(1,isp))
        call dscal(nevec*numq,wgt,ewgt,1)

C ... Case band weights are passed
      else
        if (numq /= 1) call rx('addbrl: incompatible numq')
        call dcopy(nevl,wtkb(1,isp,iq),1,ewgt,1)
        k = s_bz%nkp
        do  i = 2, s_bz%nef
          call dcopy(nevl,wtkb(1,isp,iq+(i-1)*k),1,ewgt(1,i),1)
        enddo
        do  i = nevl, 1, -1
          nevec = i
          if (abs(ewgt(i,1)) > 1d-6) exit
        enddo
      endif
      if (s_optic%loptic > 0) then
        nevec = min(nevl,max(nevec,s_optic%unrng(2)))
      endif

C ... Add to sampling dos
      if (ldos == 1) call addsds(nevl,evl(1,isp),abs(wtkp(iq))/nsp,
     .  emin,emax,esmear,ndos,dos(1,1,isp))

C ... Force from smooth analytic hamiltonian and overlap
      if (lfrce > 0 .and. lrout > 0) then
        call rxx(nspc /= 1,'forces not implemented in noncoll case')
        vavg = vconst
c        print *, '!! avg=1=smpot=1'; vavg=1; smpot=1 !; nevec=1
C        print *, '!! smpot=vavg'; smpot=vavg
C       call zprm('evecs',2,evec,ndimh,ndimh,nevec)

        if (nlmto > 0) then
          call fsmbl(nbas,s_site,s_spec,s_lat,vavg,qs,ndimh,nlmto,
     .      s_ham%iprmb,numq,nevec,evl(1,isp),evec,ewgt,f)
        endif
C        print *, 'after fsmbl'
C        do  i = 1, 3
C          print 543, f(:,i)
C  543     format(1p3e15.7)
C        enddo

        if (napw > 0) then
          vol = s_lat%vol
          call fsmbpw(nbas,s_site,s_spec,vavg,ndimh,nlmto,s_ham%iprmb,
     .      numq,nevec,evl(1,isp),evec,ewgt,napw,qpgv,qpg2v,ylv,nlmax,
     .      lmxax,alat,dsqrt(vol),f)
C        print *, 'after fsmblpw'
C        do  i = 1, 3
C          print 543, f(:,i)
C        enddo
        endif
      endif

C     print *, '!!'; call dpzero(f,9)

C ... Add to smooth density, forces, optical matrix elements
      if (lrout > 0) then
        lfrce0 = lfrce
        call rsibl(s_site,s_spec,s_lat,s_optic,lfrce0,nbas,isp,qs,iq,
     .    ndimh,nspc,napw,igapwl,s_ham%iprmb,numq,nevec,evec,ewgt,k1,k2,
     .    k3,smpot,smrho,f)
C       call zprm3('smrho after rsibl',0,smrho(1,1,1,isp),k1,k2,k3)
      endif
C      print *, 'after rsibl'
C      do  i = 1, 3
C        print 543, f(:,i)
C      enddo

C ... Add to local density coefficients
      if (lrout > 0) then
        call rlocbl(s_site,s_spec,s_lat,lfrce,nbas,isp,qs,ndham,ndimh,
     .    nspc,napw,igapwl,s_ham%iprmb,numq,nevec,evec,ewgt,evl,
     .    lcplxp,lekkl,f)
C       call zprm3('smrho after rlocbl',0,smrho(1,1,1,isp),k1,k2,k3)
      endif

C      print *, 'after rlocbl'
C      do  i = 1, 3
C        print 543, f(:,i)
C      enddo


C ... Weights for spin moments
      if (lswtk > 0 .and. nspc == 2) then
        if (ndimhx /= nevl) then
          call info0(30,0,0,' addrbl: eigenvector matrix not square'
     .      //' ... spin weights not evaluated')
          lswtk = -2
        else
        allocate(evecc(ndimh,2,ndimh,2),work(ndimh,2,ndimh,2))
        call zcopy(ndimhx**2,evec,1,evecc,1)
        call zgetrf(nevl,nevl,evecc,ndimhx,ipiv,i)
        if (i /= 0) call rx('addrbl: failed to generate overlap')
        call zgetri(nevl,evecc,ndimhx,ipiv,work,ndimhx**2,i)
        do  i = 1, ndimh
        do  k = 1, ndimh
          swtk(i,1,iq) = swtk(i,1,iq) + evecc(i,1,k,1)*evec(k,1,i,1)
     .                                - evecc(i,1,k,2)*evec(k,2,i,1)
          swtk(i,2,iq) = swtk(i,2,iq) + evecc(i,2,k,1)*evec(k,1,i,2)
     .                                - evecc(i,2,k,2)*evec(k,2,i,2)
        enddo
C        print 345,i,swtk(i,1,iq),swtk(i,2,iq), swtk(i,1,iq)+swtk(i,2,iq)
C  345   format(i4,3f14.8)
        enddo
        deallocate(evecc,work)
        endif
      endif

C ... Site resolution of the SO coupling ...
      if (ltso /= 0 .and. nspc == 2) then
C       Overwrite ewgt(:,2) with wgt'(e) and ewgt(:,3) with wgt''(e)
C        def = s_bz%def
C        do  i = 1, nevl
C          dw1 = (ewgt(i,2) - ewgt(i,3))/(2*def)
C          dw2 = (ewgt(i,2) + ewgt(i,3) - 2*ewgt(i,1))/(def**2)
C          if (abs(dw1) > 1d-8) then
C            print 333, i,evl(i,1),dw1,dw2,
C     .      ewgt(i,1)+dw1*def+dw2/2*def**2-ewgt(i,2),
C     .      ewgt(i,1)-dw1*def+dw2/2*def**2-ewgt(i,3)
C  333     format(i4,3f14.8,2x,3f14.8)
C          endif
C        enddo

C        allocate(wtfitp(6,nevl))
C        if (iq == 68) then
C        call sodwgt(s_bz%ef,s_bz%def,nevl,evl,wtkp(iq)/2,ewgt,wtfitp)
C        endif
        call sosite(0,s_site,s_spec,s_lat,s_ham,s_bz,q,ndimh,
     .    nevl,evl,evec,wtkp(iq)/nsp,ewgt,tso)
      endif

C ... Density matrix ...
C     Note: inner product can be shortened for states ewgt=0.
      if (ldmatk == 1) then
C       call dpzero(dmatk,ndimhx)
        call dscal(ndimhx,1/(abs(wtkp(iq))/nsp),ewgt,1)
        call dcopy(ndimhx,ewgt,1,dmatk,1)
      endif

      deallocate(qpgv,qpg2v,ylv)
      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
C       if (.not. associated(igapwl,s_lat%igv2)) deallocate(igapwl)
        deallocate(igapwl)
      endif

      call tcx('addrbl')
      end

      subroutine addsds(ndimh,evl,wgt,emin,emax,esmear,ndos,dos)
C- Add to sampling dos
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimh :hamiltonian dimension
Ci   evl   :eigenvalues
Ci   wgt   :eigenvalue weights
Ci   emin  :lower bound for DOS
Ci   emax  :upper bound for DOS
Ci   esmear:Parameter that describes gaussian broadening.
Ci         :Integer part >0 for for generalized gaussian broadening
Ci         :and is the the Methfessel-Paxton integration order
Ci         :Fractional part is the broadening width.
Ci         :Integer part <0 => Fermi-Dirac broadening used
Ci         :Fractional part is the temperature
Ci         :(see delstp.f)
Ci         :integer part above 100's digit is stripped.
Ci   ndos  :dimensions dos
Co Outputs
Co   dos   :DOS accumulated for these eigenvalues
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Jan 05 Extension of esmear to Fermi distribution
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ndimh,ndos
      double precision evl(*),dos(ndos,2),wgt,emin,emax,esmear
C ... Local parameters
      integer nord,ie,i1,i2,i
      double precision width,de,eigval,ei,sn,dn,fn,x

      width = dabs(esmear) - int(dabs(esmear))
      nord = dsign(1d0,esmear) * mod(int(dabs(esmear)),100)
      de = (emax-emin)/(ndos-1)
      do  ie = 1, ndimh
        eigval = evl(ie)
        i1 = (eigval-emin-width*5d0)/de + 1
        i2 = (eigval-emin+width*5d0)/de + 2
        i1 = max0(i1,1)
        i1 = min0(i1,ndos)
        i2 = max0(i2,1)
        i2 = min0(i2,ndos)
        do  i = i1, i2
          ei = emin + (i-1)*de
          x = (eigval-ei)/width
          call delstp(nord,x,dn,fn,sn)
          dos(i,2) = dos(i,2) + wgt*fn
          dos(i,1) = dos(i,1) + (wgt/width)*dn
        enddo
        do i = i2+1,ndos
          dos(i,2) = dos(i,2) + wgt
        enddo
      enddo
      end

C      subroutine snot(linv,n,u)
C      integer linv,n
C      double complex u(n,2,n,2)
C
C      complex(8) :: ul0(2,2)
C      double precision eula(3)
C
C      eula(1) = .9d0; eula(2) = .5d0 ; eula(3) = .2d0
C      call rotspu(0,1,1,1,1,eula,1,ul0)
C      if (linv == 1) call zinv22(ul,ul)
C
C      call dpzero(u,n*n*4*2)
C      do  i = 1, n
C        do  i1 = 1, 2
C        do  i2 = 1, 2
C          u(i,i1,i,i2) = ul0(i1,i2)
C        enddo
C        enddo
C      enddo
C      end
C
C      subroutine sodwgt(ef,def,nevl,evl,wtkp,ewgt,wtfitp)
CC- Build interpolating function for change in BZ weights w/ change in Ef
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ef    :Fermi energy
CCi   def   :weights are given for Ef+def and Ef-def
CCi   nevl  :number of eigenvalues
CCi   evl   :eigenvalues
CCi   wtkp  :weight of this k-point, without spin degeneracy
CCi   ewgt  :BZ weights for each eigenvalue, with the property that
CCi         :ewgt/wtkp = 1 if evl far below Ef
CCi         :ewgt/wtkp = 0 if evl far above Ef
CCo Outputs
CCo  wtfitp :Coefficients to interpolating function, for each evl
CCo          index    meaning
CCo            1       wtkp
CCo            2       sig, or 0 if wt is constant wtkp
CCo            3       mu
CCo            4       x0
CCo            5       A
CCo            6       B
CCl Local variables
CCl         :
CCr Remarks
CCr   This routine returns coefficients wtfitp to a smooth function
CCr   that interpolates changes in weights with changes in Ef.
CCr   call wtitrp below to evaluate the shift of a weight for
CCr   one eigenvalue, given its new value.
CCr
CCr  This routine builds an interpolating function for changes in BZ
CCr  integration weights arising from changes in Ef.  This is a proxy to
CCr  estimate changes in weights cause by a perturbation to the eigenvalue.
CCr  For any eigvalue 3 energies E1,E2,E3 are given, E1 = evl, E2=E1+def, E3=E1-def
CCr  with weights w1,w2,w3 that would obtain if the Fermi level were at
CCr  the actual Fermi level Ef, or Ef+def, or at Ef-def,
CCr  and an interpolating function is build to smoothly interpolate weights.
CCr  Interpolation is set up only for states near Ef, which means that
CCr  at least one of w1,w2,w3 must be sufficiently far from the asymptote (0 or 1)
CCr
CCr  Function uses Q-function to smoothly interpolate between 1 and 0
CCr    Q[(m-E)/sig] = 1-phi[(E-m)/sig],  phi(x) = cumulative normal distribution
CCr  sig is the width around Ef. To estimate, use the following ansatz:
CCr     Choose k = 2 or 3, whichever makes |w(k)-w(1)| largest,
CCr     i.e. biggest excursion in w from central point.
CCr     Make Q(x) where x=0 is midway between these energies.  Then:
CCr     Q(def/2sig) - Q(-def/2sig) = -|w(k)-w(1)|
CCr     2*Q(def/2sig) - 1 = -|w(k)-w(1)| => Q(def/2sig) = (1-|w(k)-w(1)|)/2
CCr     sig = def/[2*Q^-1((1-|w(k)-w(1)|)/2)] where Q^-1 is inverse of Q
CCr   This is reasonable provided |w(k)-w(1)| is not too small.
CCr   If it is small, Q^-1 -> 0 making sig large.  So, we limit sig as:
CCr      sig = min(sig_max,def/[2*Q^-1((1-|w(k)-w(1)|)/2)])
CCr   For now, sig_max is hardwired to 10 mRy.
CCr
CCr  Given sig, fit wl(E) = ewgt/wtkp to this form:
CCr    wl(E) = [1 + (Ax + Bx^2)exp(-x^2/2)] phi((mu-E)/sig)
CCr  where x=(x0-E)/sig and A and B are coefficients to be determined.
CCr  x0 is whichever of E1,E2,E3 have a weight closest to 1/2.
CCr  mu is the center of the normal distribution.
CCr  At E=x0, [..]=1 so wl=phi((mu-E)/sig).  This determines mu.
CCr  A and B are determined by fitting wl to the remaining two weights.
CCu Updates
CCu   30 Jan 14 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nevl
C      double precision ef,def,wtkp,ewgt(nevl,3),evl(nevl),wtfitp(6,nevl)
CC ... Local parameters
C      integer i,k,m,nk123,k123(3),k0
C      double precision dw1,dw2,rs,wl(3),sig,cdfi,cdf,xl(3),pl(3),x0,mu
C      double precision cdfix,A,B,w123(3),x123(3),rhs(2),m22(2,2),x
C      real(8), parameter :: sig_max = .01d0
C
C      do  i = 1, nevl
C        wl(1:3) = ewgt(i,1:3)/wtkp
C        dw1 = (wl(2) - wl(3))/(2*def)  ! dw/dEf
C        wtfitp(1,i) = wtkp; wtfitp(2,i) = 0 ! default: weight is constant
C        if (dabs(dw1) < 1d-4) cycle
C
CC       Construct sigma (see Remarks)
C        if (1-abs(wl(2)) < abs(wl(3))) then
C          k = 3  ! Ef -> Ef-def is farthest from limiting case
C        else
C          k = 2  ! Ef -> Ef+def is farthest from limiting case
C        endif
C        if (wl(k)-wl(1) == 0) then
C          sig = sig_max
C        else
C          x = min(-cdfi((1-dabs(wl(k)-wl(1)))/2)*2/def,1/sig_max)
C          sig = 1/x
C        endif
C        wtfitp(2,i) = sig
C
CC       Debugging: should be zero
C        print *, cdf(def/sig/2)-cdf(-def/sig/2)-dabs(wl(k)-wl(1))
C
CC   ... Identify which of E1,E2,E3 has wt closest to 1/2
C        if (dabs(wl(1)-.5d0) < dabs(wl(2)-.5d0) .and.
C     .      dabs(wl(1)-.5d0) < dabs(wl(3)-.5d0)) then
C          k0 = 1
C          x0 = evl(i)
C        elseif (dabs(wl(2)-.5d0) < dabs(wl(3)-.5d0)) then
C          k0 = 2
C          x0 = evl(i)-def
C        else
C          k0 = 3
C          x0 = evl(i)+def
C        endif
C        w123(1) = wl(k0)
C        mu = cdfix(wl(k0))*sig + x0  ! Must be stored
C        wtfitp(3,i) = mu
C        wtfitp(4,i) = x0
CC       Debugging: should be zero
C        print *, cdf((mu-x0)/sig)-w123(1)
C
CC       Count the number of energies to include in the fit
CC       Exclude points too close to asymptote
C        nk123 = 1 !  number of qualifying points
C        k123(1) = k0
C        w123(1) = wl(k0)
C        if (dabs((1-wl(k0))*wl(k0)) < 1d-4)
C     .    call rx('sodwgt: problem with weights')
C        do  k = 1, 3
CC         Include them all?
CC         if (dabs((1-wl(k))*wl(k)) < 1d-4 .or. k == k0) cycle ! necessary?
C          if (k == k0) cycle
C          nk123 = nk123+1
C          k123(nk123) = k
C          w123(nk123) = wl(k)
C
C          if (k == 1) x123(nk123) = evl(i)
C          if (k == 2) x123(nk123) = evl(i)-def
C          if (k == 3) x123(nk123) = evl(i)+def
C        enddo
C
CC   ... Fit to the form wl = [1 + (Ax + Bx^2)exp(-x^2/2)] phi((mu-E)/sig)
CC       where x=(x0-E)/sig.  mu,x0, and sig are already known: make A and B
C        A = 0; B = 0; if (nk123 == 1) goto 10 ! only 1 qualifying point: A = B = 0
C
CCC         w123(2) = [1 + A*x*exp(-x^2/2)] phi((mu-x123(2))/sig)
CCC         A = [w123(2)/phi((mu-x123(2))/sig) - 1]/[x*exp(-x^2/2)]
CCC         wl([1 + A*x*exp(-x^2/2)] phi((mu-E)/sig)
CC          x = (x123(1)-x123(2))/sig
CC          A = (w123(2)/cdf((mu-x123(2))/sig) - 1)/(x*exp(-x^2/2))
C
CC       w123(2)/phi((mu-x123(2))/sig) - 1 = x*exp(-x^2/2) A + x*exp(-x^2/2) B
CC       w123(3)/phi((mu-x123(3))/sig) - 1 = y*exp(-y^2/2) A + y*exp(-y^2/2) B
CC       If only one additional point, set B = 0 and use top equation only
C        x = (x123(1)-x123(2))/sig
C        m22(1,1) = x*exp(-x**2/2); m22(1,2) = x*m22(1,1)
C        rhs(1) = w123(2)/cdf((mu-x123(2))/sig) - 1
C        if (nk123 == 2) then ! 2 qualifying points: B = 0
C          A = rhs(1)/m22(1,1)
C          goto 10
C        endif
CC       Solve two simultaneous equations
C        x = (x123(1)-x123(3))/sig
C        m22(2,1) = x*exp(-x**2/2); m22(2,2) = x*m22(1,1)
C        rhs(2) = w123(3)/cdf((mu-x123(3))/sig) - 1
C        call dinv22(m22,m22)
C        A = m22(1,1)*rhs(1) + m22(1,2)*rhs(2)
C        B = m22(2,1)*rhs(1) + m22(2,2)*rhs(2)
C
CC       At this point mu,sig,A,B,x0 all have been determined
C   10   continue
C
C
C        xl(1) = 0; xl(2) = -def; xl(3) = +def
C        do  k = 1, 3
C          xl(k) = -(evl(i)-ef+xl(k))/sig
C          pl(k) = cdf(xl(k))
C        enddo
C        print *, evl(i),ef
C        print *, wl
C        print *, xl
C        print *, pl
C
C        stop
C
C        print *, exp(-(evl(i)-ef)**2/(rs**2))
C        print *, derfc(0d0)/2
C        print *
C        print *, def,rs
C        print *,
C     .    derfc((evl(i)-def-ef)/rs)/2,
C     .    derfc((evl(i)-ef)/rs)/2,
C     .    derfc((evl(i)+def-ef)/rs)/2
C        print *, sngl(wl)
C        stop
C
C        dw2 = (wl(2) + wl(3) - 2*wl(1))/(def**2)
C        if (abs(dw1) > 1d-8) then
C          print 333, i,evl(i),dw1,dw2,
C     .      wl(1)+dw1*def+dw2/2*def**2-wl(2),
C     .      wl(1)-dw1*def+dw2/2*def**2-wl(3)
C  333     format(i4,3f14.8,2x,3f14.8)
C        endif
C      enddo
C
C      end
C
C      subroutine wtitrp(evl,wtfitp,wtfit)
CC- Return energy interpolated BZ weight
CC ----------------------------------------------------------------------
CCi Inputs
CCi   evl   :eigenvalue
CCi   wtfitp:coefficients to fit; see sodwgt
CCo Outputs
CCo   wtfit :BZ weigh for evl
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr   Functional form developed in sodwgt above
CCu Updates
CCu   31 Jan 14 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      double precision evl,wtfitp(5),wtfit
CC ... Local parameters
C      double precision wtkp,sig,mu,x0,A,B,x,cdf
C
C      wtkp = wtfitp(1)
C      sig  = wtfitp(2)
C      mu   = wtfitp(3)
C      x0   = wtfitp(4)
C      A    = wtfitp(5)
C      B    = wtfitp(6)
C
C      x = (x0-evl)/sig
C      wtfit = (1 + (A*x + B*x*x)*exp(-x*x/2))* cdf((mu-evl)/sig)
C
C      end
