      subroutine pvgfe1(ng,ips0,gv,cJ,csym,tau,rcut,Jsum,ib,jb,ifi,scl,
     .  tol,lwj00)
C- Print out J(T) for one pair of atoms and translation vectors T
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of lattice vectors
Ci   ips0  :pointer to first vector in star of this vector; see sgvsym
Ci   gv    :list of lattice vectors
Ci   cJ    :unsymmetrized Heisenberg parameters J
Ci   csym  :symmetrized Heisenberg parameters J
Ci   tau   :connecting vector for between
Ci   rcut  :largest distance for which to print.
Ci         :rcut<0 => use sign to zero out any cJ for r>rcut
Ci   ib,jb :site indices for this pair
Ci         :ifi>0 ib,jb written to file containing J
Ci   ... the next four tokens deal with writing Ham to file
Ci   ifi   :if >0, write r.s. hamiltonian to file unit ifi
Ci   scl   :factor with which to scale J
Ci   tol   :tolerance; pairs for which J<tol are not written
Ci   lwj00 :(writing to rsj file)
Ci         :1 write J00 into rsj file
Ci         :NB: If tol=0, the sum of rsj should equal Jsum(3,ib)
Co Outputs
Co   jsum  :jsum(3,ib) is added into for pairs within range rcut
Cl Local variables
Cl         :
Cr Remarks
Cr   Prints out J for only one element in the star
Cu Updates
Cu   22 Jan 13 Modifications for CPA case
Cu   12 Apr 03 Better treatment of rcut
Cu   22 Nov 02 Added ability to write r.s. H to file ifi
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,ips0(ng),ib,jb,ifi,lwj00
      double precision gv(ng,3),rcut,tau(3),scl,tol,Jsum(3)
      double complex cJ(ng), csym(ng)
C ... Local parameters
      integer nstar,i0,kstar,i,ndeg,ipr,stdo,nglob
      double precision vv,rms,fac,v(3),sumJ,sumJR2

      stdo = nglob('stdo')
      call getpr(ipr)
      nstar = 0
      fac = 1000
      sumJ = 0
      sumJR2 = 0
      if (ipr >= 30) write(stdo,1)
    1 format(/'   Star   nstar    |R|          -- R(jb)-R(ib) -- ',
     .'       J(mRy)   sum J')
C  250 format(/'   Star   nstar    |R|          -----   R   ----- ',
C     .'      J(mRy)   sum J  -sum J*R**2/2')

      do  i0 = 1, ng
        if (ips0(i0) == i0) then
          nstar = nstar+1

C     ... ndeg <- members of the star; need for writing
          kstar = 0
          do  i = i0, ng
            if (ips0(i) == i0) then
              kstar = kstar+1
            endif
          enddo
          ndeg = kstar

C     ... Check that J is symmetrized; optional write to file ifi
          rms = 0
          kstar = 0
          do  i = i0, ng
            if (ips0(i) == i0) then
              kstar = kstar+1
              rms = rms + dble(dconjg(csym(i)-cJ(i))*(csym(i)-cJ(i)))
C         ... Make connecting vector = rtab = (pos(jb=src)-pos(ib=field))
C             v(1) = gv(i,1) - tau(1)
C             v(2) = gv(i,2) - tau(2)
C             v(3) = gv(i,3) - tau(3)
              v(1) = tau(1) - gv(i,1)
              v(2) = tau(2) - gv(i,2)
              v(3) = tau(3) - gv(i,3)
              vv = v(1)**2+v(2)**2+v(3)**2
              if (vv < rcut**2) then
                Jsum(3) = Jsum(3) + dble(cJ(i))
              endif
              if (ifi /= 0) then
                if (abs(fac*scl*dble(cJ(i))) > tol .and.
     .              vv < rcut**2) then
                  if (vv /= 0 .or. mod(lwj00,2) /= 0) then
                    if (kstar == 1) then
                      write(ifi,2) ib,jb,v,fac*scl*dble(cJ(i)),ndeg
                    else
                      write(ifi,2) ib,jb,v,fac*scl*dble(cJ(i))
                    endif
                  endif
                endif
              endif

            endif
          enddo

          v(1) = gv(i0,1) - tau(1)
          v(2) = gv(i0,2) - tau(2)
          v(3) = gv(i0,3) - tau(3)
          v(1) = tau(1) - gv(i0,1)
          v(2) = tau(2) - gv(i0,2)
          v(3) = tau(3) - gv(i0,3)
          vv = dsqrt(v(1)**2+v(2)**2+v(3)**2)
          if (vv > abs(rcut)) then
            if (rcut < 0) then
              call dpzero(cJ(i0),(ng-i0+1)*2)
              call dpzero(csym(i0),(ng-i0+1)*2)
            endif
            return
          endif
          sumJR2 = sumJR2 - kstar*vv**2*dble(csym(i0))/2
          if (vv /= 0) sumJ = sumJ + kstar*dble(csym(i0))
          if (ipr >= 30) then
          if (rms > 1d-16) then
C           write(stdo,3) nstar,kstar,vv,v,fac*dble(csym(i0)),fac*sumJ,fac*sqrt(rms)
            call info8(30,0,0,
     .        '%,5i%,8i%;12,6D%3;9,3D%;12,5D%;9,3D%;14,3D',
     .        nstar,kstar,vv,v,fac*dble(csym(i0)),fac*sumJ,fac*sqrt(rms),0)
          else
C           write(stdo,3) nstar,kstar,vv,v,fac*dble(csym(i0)),fac*sumJ
            call info8(30,0,0,
C    .        '%,5i%,8i%;12,6D%3;9,3D%;12,5D%;9,3D%;14,3D',
     .        '%,5i%,8i%;12,6D%3;9,3D%;12,5D%;9,3D',
     .        nstar,kstar,vv,v,fac*dble(csym(i0)),fac*sumJ,0,0)
          endif

          if (rms > 1d-8) then
            write(stdo,*) '(warning) symmetrized J not equal to J'
          endif
          endif
        endif
      enddo
    2 format(2I4,3F12.7,2x,f12.6:i6)
    3 format(i5,i8,f12.6,3f9.3,f12.5,f9.3,5x,f9.3)

      end
      subroutine pvgfe2(isp,nsp,nq,offi,offj,ni,nj,npf,
     .  dpf,ddpf,ldi,ldj,ldiag,gij,gji)
C- Convert subblock of proper G to dimensionless g (aka T-matrix)
C ----------------------------------------------------------------------
Ci Inputs
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nq    :number of k-points on a single symmetry line for bands plotting
Ci   offi  :offset in dpf, pdf corresponding to i=1 element in G
Ci   offj  :offset in dpf, pdf corresponding to j=1 element in G
Ci   ni    :size of row subblock to convert
Ci   nj    :size of col subblock to convert
Ci   npf   :dimensions dpf,ddpf
Ci   dpf   : dot-P, where P is potential function
Ci   ddpf  : dot-dot-P, where P is potential function
Ci   ldi   : dimensions gij,gji
Ci   ldj   : dimensions gij,gji
Ci   ldiag : if T, just add diagonal part
Cio Inputs/Outputs
Cio  gij  : On input, subblock of Green's function Gij.  On output, Tij
Cio  gji  : On input, subblock of Green's function Gji.  On output, Tji
Cl Local variables
Cl         :
Cr Remarks
Cr  No downfolding now.
Cu Updates
Cu   06 Apr 16 dpf is now dot-P (change from old sqrt(dot-P))
C ----------------------------------------------------------------------
      implicit none
      logical ldiag
      integer isp,nsp,nq,ni,nj,npf,ldi,ldj
      integer offi,offj
      double complex dpf(npf,nsp),ddpf(npf,nsp)
      double complex gij(nq,ldi,ldj,nsp),gji(nq,ldi,ldj,nsp)
C Local variables
      integer id,jd,iq
      double complex xi,xj,di

      call tcn('pvgfe2')

C --- Loop over orbital pairs ---
      do  id = 1, ni
      do  jd = 1, nj

C        xi = 1/dpf(id+offi,isp)
C        xj = 1/dpf(jd+offj,isp)

        xi = 1/sqrt(dpf(id+offi,isp))
        xj = 1/sqrt(dpf(jd+offj,isp))

        do  iq = 1, nq
          gij(iq,id,jd,isp) = xi*gij(iq,id,jd,isp)*xj
        enddo
        if (id+offi == jd+offj) then
          di = ddpf(id+offi,isp)*xi**2
          do  iq = 1, nq
            gij(iq,id,jd,isp) = gij(iq,id,jd,isp) - di
          enddo
        endif

        if (.not. ldiag) then
          do  iq = 1, nq
            gji(iq,jd,id,isp) = xj*gji(iq,jd,id,isp)*xi
          enddo
        if (id+offi == jd+offj) then
          do  iq = 1, nq
            gji(iq,jd,id,isp) = gji(iq,jd,id,isp) - di
          enddo
        endif
        endif

      enddo
      enddo

      call tcx('pvgfe2')

      end

      subroutine pvgfe3(plat,n1,n2,n3,Jq,D1,D2)
C- Extract q->0 limit of J
C  Returns D1 = 4-point estimate for A as Jq = A * q**2
C          D2 = diff between 2pt and 4-point estimate
C  Here q is in units of (2 pi/alat)
      implicit none
      integer n1,n2,n3
      double precision plat(3,3),D1,D2
      double complex Jq(n1,n2,n3)
C Local variables
      integer is(3),ifac(3),k,jj1,jj2,jj3
      double precision qb(3,3),rb(3,3),qk,q2,D12
      logical llshft(3)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./

      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr
C     tpi = 8*datan(1d0)

C ... 1st order estimate for D: use first point
      q2 = qk(1,2,1,1)**2 + qk(2,2,1,1)**2  + qk(3,2,1,1)**2
      D1 = dble(Jq(2,1,1)-Jq(1,1,1)) / q2
C ... 2nd order estimate for D1 uses first two points
      D12 = dble(Jq(3,1,1)-Jq(1,1,1)) / q2
      D12 = 4d0/3*D1 - 1/12d0*D12
C ... D1 <- 2nd order estimate for q**2 coff; D2 <- q**4 coff
C      print *, D1, D12, q2
      D2 = D1 - D12
      D1 = D12

      end

      subroutine pvgfe4(plat,n1,n2,n3,wq,tcfac,opqres,wint)
C- Integral omega(q), 1/omega(q) for MF and RPA estimates of Tc
C ----------------------------------------------------------------------
Ci Inputs
Ci   plat  :primitive lattice vectors, in units of alat
Ci  n1..n3:number of k-divisions, 3 lattice vectors
Ci   wq    :SW frequencies
Ci   tcfac:scale factor to convert q sum to Tc; see Remarks
Ci   opqres: string specifying qres options; See remarks
Co Outputs
Co   wint :MF and RPA estimtes of Tc.
Cr Remarks
Cr   Note: for many-atom case, see pvgfek.
Cr
Cr   Tyablikov formula for Tc is, see e.g. DOI 10.1007/s10948-009-0569-3
Cr    kTc = 2/3 [1/N sum_q [J(q)-J(0)]^-1]^-1
Cr   Can write as (Kudrnovsky PRB 64, 174402, Eq. 10)
Cr     (kB Tc)^-1 = 6/M [1/N sum_q 1/omega(q)]
Cr   The Tyablikov, RPA, and spherical model formulas are identical.
Cr   Compare to mean-field formula
Cr     kTc = 2/3 1/N sum_q [J(q)-J(0)]
Cr   Here we calculate integrals
Cr     [1/N sum_q 1/w(q)]^-1
Cr     [1/N sum_q w(q)]
Cr   The caller must supply the appropriate scale factor to convert
Cr   these integrals to an estimate of Tc.  For example for 1 atom/cell,
Cr   using lmgf conventions where magnetic moments are folded into
Cr   the exchange, the relation  between omega and J is
Cr         omega = M/4 [J(q)-J(0)]; see pvgfe8.
Cr
Cr   The integral of 1/w(q) is calculated as the gam->0 limit of:
Cr    1/N sum_q 1/w(q) = 1/N lim(gam->0) sum_q w(q)/ (w^2(q) + gam^2)
Cr   The integral is calculated for several values of gam, and
Cr   extrapolated to gam=0 with a quadratic least-squares fit.
Cr
Cr  optres resolves q-point contributions to Tc to locate for what q
Cr  the contributions to Tc are strongest
Cr  If optres starts with --tcqres these integrals are performed:
Cr     I^MF(qcut)  = 1/N sum_q d((q-qcut)/gausw) * w(q)
Cr     I^RPA(qcut) = 1/N sum_q d((q-qcut)/gausw) / w(q)
Cr     where d a normalized gaussian.
Cr  Sample confirmation that the q-resolved integrates to the total:
Cr    mc tcqres.fe -int 0 1.4
Cu Updates
Cu   12 Jul 10 Tyablikov formula for Tc, 1 atom/cell
C ----------------------------------------------------------------------
      implicit none
      double precision plat(3,3),tcfac
      integer n1,n2,n3
      double precision wq(n1,n2,n3),wint(2)
      character opqres*128
C Local variables
      integer i,n,ierr,m,j,iq
      double precision f(4),gam(4),xij(3,3),yi(3),w(3,3),fmf,fac,fgeo
      integer is(3),ifac(3),k,jj1,jj2,jj3,i3,i2,i1,ifi,fopna
      double precision qb(3,3),rb(3,3),qk,q1(3),q,ddot,qmax,x,qcut,qwsr
      double precision wqcut(100,4),qlat(9),
     .  vol0,gausw0,gausw,srpi,d
C     double precision s,derfc
      parameter (gausw0=0.04d0)
      logical llshft(3),lqres
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./
      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr
      call mkqlat(plat,qlat,vol0)
      qwsr = (3d0/vol0/16/datan(1d0))**(1d0/3d0)
      qmax = 1.7d0*qwsr
      gausw = qwsr*gausw0
      srpi = dsqrt(4d0*datan(1d0))
C     stdo = nglob('stdo')
      lqres = opqres(1:8) == '--tcqres'
      if (lqres) then
        ifi = fopna('tcqres',-1,0)
        rewind ifi
        write(ifi,1) qwsr
    1   format('#      q',9x,'Tc_RPA',6x,'Tc_MF',7x,'fracbz',
     .    '   qwsr =',f9.5)
      endif

C      dc = opqres(1:1)
C      if (dc /= ' ') then
C      endif

C ... RPA summation for four values of gam; take limit for result
      gam(1) = (wq(2,1,1) - wq(1,1,1))/2
      gam(2) = 2*gam(1)
      gam(3) = 3*gam(1)
      gam(4) = 4*gam(1)
      n = n1*n2*n3
      do  iq = 0, 100
      qcut = 1.0d0*dble(iq-10)/100 * qmax
C      x = qcut/qmax
C      print *, x, qcut
C      enddo
C      x=(qcut-qmax)*5/qmax
C      print *, x, 1/(exp(x)+1)


      call dpzero(f,4)

      fmf = 0
      fac = 1
      fgeo = 1
C ... Do the integral for this qcut
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        fac = 1
        d = 1
C       s = 1

C          q1(1) = qk(1,i1,i2,i3)
C          q1(2) = qk(2,i1,i2,i3)
C          q1(3) = qk(3,i1,i2,i3)
C          print "(3i4,3f12.5)", i1,i2,i3, q1


        if (iq /= 0) then
          q1(1) = qk(1,i1,i2,i3)
          q1(2) = qk(2,i1,i2,i3)
          q1(3) = qk(3,i1,i2,i3)
          call shorbz(q1,q1,qlat,plat)
          q = dsqrt(ddot(3,q1,1,q1,1))
          x=(q-qcut)/gausw
C         d = sm-delta-function; s is integral of d
          d = dexp(-x*x)/srpi/gausw
C         s = .5d0 * derfc(x)
C          fac = 1/(exp(x)+1)
C          fac = -exp(x)*fac*fac * 25/qmax
C         if (q == 0 .or. i1*i2*i3 == n) print *, iq,sngl(q),x,fac
C         fac = 1
        endif
        fgeo  = fgeo + d
        fmf   = fmf + d*wq(i1,i2,i3)
        f(1) = f(1) + d*wq(i1,i2,i3)/(wq(i1,i2,i3)**2 + gam(1)**2)
        f(2) = f(2) + d*wq(i1,i2,i3)/(wq(i1,i2,i3)**2 + gam(2)**2)
        f(3) = f(3) + d*wq(i1,i2,i3)/(wq(i1,i2,i3)**2 + gam(3)**2)
        f(4) = f(4) + d*wq(i1,i2,i3)/(wq(i1,i2,i3)**2 + gam(4)**2)
      enddo
      enddo
      enddo
      call dscal(4,1/dble(n),f,1)
      fmf = fmf/dble(n)
      fgeo = fgeo/dble(n)

C       print *, fmf,
C     .  fmf*tcfac, fmf*tcfac / 8.617d-5
CC       print *, 2d0/3d0 * 13.6d0 * 2.263322/4, tcfac
C       stop

C ... Extrapolate to the origin
      call dpzero(xij,3*3)
      call dpzero(yi,3)
      do  i = 1, 3
        do  m = 1, 4
          do  j = 1, 3
            xij(i,j) = xij(i,j) + gam(m)**(i-1+j-1)
          enddo
          yi(i) = yi(i) + gam(m)**(i-1) * f(m) / 1d0
        enddo
      enddo
C      call prmx('xij',xij,3,3,3)
C      print *, yi
      call dqinvb('s',xij,3,0,3,1,w,3,w,yi,3,ierr)
      if (ierr /= 0) call rx('pvgfe4: bug at dqinvb')

      if (iq == 0) then
        wint(1) = tcfac/yi(1)
        wint(2) = tcfac*fmf
      else
        wqcut(iq,1) = qcut
        wqcut(iq,2) = tcfac/yi(1)
        wqcut(iq,2) = yi(1)/tcfac * wint(1)**2
        wqcut(iq,3) = tcfac*fmf
        wqcut(iq,4) = fgeo

        if (qcut > 0)
     .  write(ifi,"(f12.6,3f12.6)") (wqcut(iq,i), i=1,4)
      endif

C      call info0(20,0,0,'%N pvgfe4: extrapola I(gam) = '//
C     .  'int d^3q wq/(wq^2+gam^2): to gam=0: %N    gam         I')
C      do  40  i = 1, 4
C   40 print '(1pe12.4,0p,f9.3)', gam(i),f(i)
C      print "(' extrapolate',f9.3)", yi(1)
C      exit

      if (.not. lqres) exit
      enddo
      if (lqres) call fclose(ifi)

      end
      subroutine pvgfe5(isp,nsp,iq,ldim,h1,evec,s)
C- Write  s*z into sz
      implicit none
      integer ldim,isp,nsp,iq
      double precision h1(ldim,ldim,2),err
      double precision evec(ldim,ldim,2,nsp,*),s(ldim,ldim,2,nsp,*)
      integer stdo,lgunit,i,j

C     call yprm('z',2,evec(1,1,1,isp,iq),ldim*ldim,ldim,ldim,ldim)
C     call yprm('ovlp',2,h1,ldim*ldim,ldim,ldim,ldim)

      call yygemm('N','N',ldim,ldim,ldim,1d0,h1,h1(1,1,2),ldim,
     .  evec(1,1,1,isp,iq),evec(1,1,2,isp,iq),ldim,0d0,
     .  s(1,1,1,isp,iq),s(1,1,2,isp,iq),ldim)

C     call yprm('s z',2,s(1,1,1,isp,iq),ldim*ldim,ldim,ldim,ldim)

C ... Debugging: check whether z+ s z is unity
      stdo = lgunit(1)
      call yygemm('C','N',ldim,ldim,ldim,1d0,
     .  evec(1,1,1,isp,iq),evec(1,1,2,isp,iq),ldim,
     .  s(1,1,1,isp,iq),s(1,1,2,isp,iq),ldim,0d0,
     .  h1(1,1,1),h1(1,1,2),ldim)
      err = 0
      do  i = 1, ldim
        err = err + (h1(i,i,1)-1)**2 + h1(i,i,2)**2
        do  j = i+1, ldim
          err = err+h1(i,j,1)**2+h1(i,j,2)**2+h1(j,i,1)**2+h1(j,i,2)**2
        enddo
      enddo
      if (err > 1d-24)
     .  call awrit1(' pvgfe5 (warning): deviation in z+ s z '//
     .  'from unity = %,3g',' ',80,stdo,dsqrt(err))

      end
      subroutine pvgfe6(n1,n2,n3,amom,J)
C- Make omega from Jq, (1 atom/cell)
      implicit none
      integer n1,n2,n3
      double precision amom,J(n1,n2,n3)
      integer i1,i23

      call info2(20,0,0,' pvgfe6: converting J to omega: amom=%,4d',
     .  amom,0)

      do  i23 = 1, n2*n3
        do  i1 = 1, n1
          J(i1,i23,1) = -4/abs(amom) * J(i1,i23,1)
        enddo
      enddo
      end

      integer function pvgfe7(s,lio,ifi,ib,icomp,jb,jcomp,n1,n2,n3)
C- I/O for exchange integrals J
C ----------------------------------------------------------------------
Ci  lio  0  read, return s in real space (see Remarks)
Ci       1  write s in real space
Ci       2  read, find start of array containing 1st (ib,jb) pair
Ci          Notes: array contents not read.  No check on components
Ci      10  read s, return s in recip space
Ci      11  write s in recip space; input s is real space
Ci      21  write s in recip space; input s is recip space
Ci     111  write s in real space; input s is real space
Ci     121  write s in real space; input s is recip space
Ci  ifi     file logical unit number
Ci  ib,jb,icomp,jcomp
Ci           pair (including components) for which to read/write J
Ci  n1..3   dimensions s
Cio Inputs/Outputs
Cio   s    array containing exchange integrals
Cr  Remarks
Cr    It is the caller's responsibility to rewind file.
Cr
Cr    File read (1s digit mode 0):
Cr    pvgfe7 assumes that (ib,jb) pairs are ordered by increasing ib,jb,
Cr    with jb=fast index.  If file ib>passed ib, or if
Cr    file jb>passed jb and file ib>passed ib, pvgfe7 returns with -1
Cr    If file end-of-file encountered, pvgfe7 returns with -2
Cr
Cr    File read (1s digit mode 2):
Cr    Looks in file for start of array that contains (ib,jb) pair.
Cr    No assumption is made about file order; no attempt is made
Cr    to read the array if it is found.  pvgfe7 returns 0 if start
Cr    is found, or -2 if end-of-file encountered.
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   16 Jan 13 Reworked site list for CPA compatibility
Cu   05 Mar 01 Added read mode 2
Cu   22 Dec 00 turned into a function call
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,ib,jb,icomp,jcomp,lio,ifi
      double precision s(n1,n2,n3)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: h(:)
C ... Local parameters
      logical a2bin  !,matchib,matchic,matchjb,matchjc
      integer i,j,ib0,jb0,icmp0,jcmp0,k,nfbz,ip,recl,iprint,lio0
      integer nglob,stdo
      parameter (recl=80)
      character*(20) fmt, space*2, recrd*(recl), dir(4)*7
      data dir /'rows','cols','rs','qs'/

      stdo = nglob('stdo')
      lio0 = mod(lio,10)
      pvgfe7 = 0
C --- Write ---
      if (lio0 == 1) then
        fmt = '(9f15.10)'
        space = 'rs'
        if (lio/100 == 0) space = 'qs'
        if (lio == 11 .or. lio == 121) then
          nfbz = n1*n2*n3
          allocate(h(nfbz)); call dpzero(h,2*nfbz)
          call dcopy(nfbz,s,1,h,2)
          i = -1
          if (lio == 11) i = 1
          call fftz3(h,n1,n2,n3,n1,n2,n3,1,0,i)
          call dcopy(nfbz,h,2,s,1)
          deallocate(h)
        endif
        call info5(41,0,0,' pvgfe7:  writing '//space//' J for '//
     .    'ib=%i%?#n# icomp=%-1j%i##, jb=%i%?#n# jcomp=%-1j%i##',
     .    ib,icomp,jb,jcomp,0)
        call awrit6('%% rows %i cols %i real  '//space//
     .    '  ib=%i%?#n#  icomp=%-1j%i##  jb=%i%?#n#  jcomp=%-1j%i##',
     .    ' ',80,ifi,n1*n2,n3,ib,icomp,jb,jcomp)
        do  i = 1, n1
          do  j = 1, n2
            write (ifi,fmt) (s(i,j,k),k=1,n3)
          enddo
        enddo
C --- Read ---
      else
C   ... Entry point for new array to read, modes 0 and 2
    5   continue
        read (ifi, '(a)',END=99) recrd
        space = ' '
        if (recrd(1:1) == '%') then
          ip = 1
          ib0 = 0; jb0 = 0; icmp0 = 0; jcmp0 = 0
    6     continue
          call skipbl(recrd,recl,ip)
C   ...   If header has been parsed, read array contents
          if (ip >= recl) then
            if (space == ' ') then
              call rx('pvgfe7: file jr not in standard format')
            endif
            do  i = 1, n1
              do  j = 1, n2
                read (ifi,*) (s(i,j,k),k=1,n3)
              enddo
            enddo
C           We are done if the indices match
            if (((icomp >= 1.or.icmp0 /= 0) .and. icomp /= icmp0) .or.
     .          ((jcomp >= 1.or.jcmp0 /= 0) .and. jcomp /= jcmp0) .or.
     .          (jb0 /= 0 .and. jb /= jb0) .or.
     .          (ib0 /= 0 .and. ib /= ib0)) goto 5
            recrd = '%N pvgfe7: read J('//space//'), returning J(  )'
            if (ib0 /= 0) call awrit1('%a, ib=%i',recrd,80,0,ib)
            if (icmp0 /= 0) call awrit1('%a, icomp=%i',recrd,80,0,icomp)
            if (jb0 /= 0) call awrit1('%a, jb=%i',recrd,80,0,jb)
            if (jcmp0 /= 0) call awrit1('%a, jcomp=%i',recrd,80,0,jcomp)
            if (lio == 0 .and. space == 'qs' .or.
     .        lio /= 0 .and. space == 'rs') then
              nfbz = n1*n2*n3
              allocate(h(nfbz)); call dpzero(h,2*nfbz)
              call dcopy(nfbz,s,1,h,2)
              i = -1
              space = 'rs'
              if (lio /= 0) i = 1
              if (lio /= 0) space = 'qs'
              call fftz3(h,n1,n2,n3,n1,n2,n3,1,0,i)

C              call zprm3('jq',2,h,n1,n2,n3)
C              call fftz3(h,n1,n2,n3,n1,n2,n3,1,0,-1)
C              call zprm3('jr',2,h,n1,n2,n3)
C              call fftz3(h,n1,n2,n3,n1,n2,n3,1,0,i)
C              call zprm3('jq',2,h,n1,n2,n3)
C
              call dcopy(nfbz,h,2,s,1)
              deallocate(h)
            endif
            recrd(36:37) = space
            if (iprint() >= 30) call awrit0(recrd,' ',-80,stdo)
            goto 999
C     ... Continue to parse header
          else
            k = ip-1
            call tokmat(recrd(ip+1:recl),dir,4,7,' ',i,ip,.false.)
            ip = ip+k
            if (i < 0) then
              if (recrd(ip+1:ip+3) == 'ib=') then
                i = 0
                if (.not. (a2bin(recrd(ip+4:),ib0,2,0,' ',i,recl)))
     .            call rxs('pvgfe7: failed to parse ib in',recrd)
C               matchib = ib == ib0
                if (lio0 == 2 .and. ib == ib0 .and. jb == jb0) then
                  backspace ifi
                  return
                endif
C               File ib0 > ib => already past ib we are looking for ...
                if (lio0 == 0 .and. ib < ib0) then
                  backspace ifi
                  pvgfe7 = -1
                  return
C                 call rx('pvgfe7 not ready')
                endif
              elseif (recrd(ip+1:ip+6) == 'icomp=') then
                i = 0
                if (.not. (a2bin(recrd(ip+7:),icmp0,2,0,' ',i,recl)))
     .            call rxs('pvgfe7: failed to parse icomp in',recrd)
C               icmp0>0 and ib0>=ib => already past this icomp ...
                if (lio0 == 0 .and. icomp < icmp0 .and. ib <= ib0)
     .            then
                  backspace ifi
                  pvgfe7 = -1
                  return
                endif
              elseif (recrd(ip+1:ip+3) == 'jb=') then
                i = 0
                if (.not. (a2bin(recrd(ip+4:),jb0,2,0,' ',i,recl)))
     .            call rxs('pvgfe7: failed to parse jb in',recrd)
C               A match to this (ib,jb) pair
                if (lio0 == 2 .and. ib == ib0 .and. jb == jb0) then
                  backspace ifi
                  return
                endif
C               File jb0 > jb and ib0 >= ib => already past this pair ...
                if (lio0 == 0 .and. jb < jb0 .and. ib <= ib0) then
                  backspace ifi
                  pvgfe7 = -1
                  return
                endif
              elseif (recrd(ip+1:ip+6) == 'jcomp=') then
                i = 0
                if (.not. (a2bin(recrd(ip+7:),jcmp0,2,0,' ',i,recl)))
     .            call rxs('pvgfe7: failed to parse jcomp in',recrd)
C               jcmp0>0 and jb0>=jb and ib0>=ib => already past this jcomp ...
                if (lio0 == 0 .and. jcomp < jcmp0 .and. jb <= jb0
     .              .and. ib <= ib0) then
                  backspace ifi
                  pvgfe7 = -1
                  return
                endif

              endif
              call skp2bl(recrd,recl,ip)
C           Check match in row (or column) dimension
            elseif (i == 0 .or. i == 1) then
              if (.not. (a2bin(recrd,j,2,0,' ',ip,recl))) return
              if (i == 0 .and. j /= n1*n2 .or. i == 1 .and. j /= n3)
     .          call rx('pvgfe7:  file mismatch')
            elseif (i == 2) then
              space = 'rs'
            elseif (i == 3) then
              space = 'qs'
            endif
          endif
          goto 6
        else
          call rx('pvgfe7: file jr not in standard format')
        endif
   99   continue
        pvgfe7 = -2
        return
C       call prm3(recrd,0,s,n1,n2,n3)
      endif
  999 continue
      end

      subroutine pvgfe8(n1,n2,n3,plat,ipc,amom,s_cpasite,nsite,nJq,
     .  JqRR,omega,nomega)
C- Make omega from Jq, multi-band case
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1..n3:mesh of q-points on which Jq is tabulated
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   amom  :magnetic moments, by class
Ci   nlst  :number of sites for which sumRL is calculated
Ci         :nlst=0 -> all sites calculated; list is not used
Ci   list  :list of sites, if nlst>0
Ci   JqRR  :exchange interactions on uniform mesh of q-points
Co Outputs
Co   omega :spin waves
Co  nomega :minimum number of spin waves calculated in the BZ.
Cr Remarks
Cr   Adapt Takao's convention: spin waves satisfy
Cr     det |omega/M_i - Jbar^T(q)_ij| = 0
Cr     S_i J^T_ij(q) S_j = J_ij(q)
Cr     J^T_ij = (2/M_i) J_ij (2/M_j)
Cr     Jbar^T_ij(q) = J^T_ij(q)
Cr                  - [(2/M_i) sum k J^T_ik(0) M_k/2] delta_ij
Cr                  = (2/M_i) J_ij (2/M_j)
Cr                  - [(2/M_i) sum k (2/M_i) J_ik(0) (2/M_k) M_k/2]
Cr                    * delta_ij
Cr                  = 4 * 1/M_i J_ij 1/M_j
Cr                  - 4 (1/M_i)^2 delta_ij  sum_k J_ik(q=0)
Cr   Thus for one atom/cell the ASA formula is (see pvgfej)
Cr      M omega / 4 = J(0) - J(q)
Cu Updates
Cu   16 Jan 13 Reworked site list for CPA compatibility
Cu   11 Aug 08 Spin waves for Ferri case (first attempt)
Cu   20 Feb 07 Spin waves for 2 atom AFM case, multiatom FM case
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer n1,n2,n3,nsite,nJq,ipc(*)
C    .  nbas,ipc(nbas),nlst,list(nbas)
      double precision amom(nJq),omega(n1,n2,n3,nJq),plat(3,3)
      double complex JqRR(n1,n2,n3,nJq,nJq)
C ... For structures
!      include 'structures.h'
      type(str_cpasite):: s_cpasite(nsite)
C ... Local parameters
      character outs*120
      logical lpos,lneg,llshft(3),cmdopt,a2bin
      integer i1,i2,i3,ib,jb,nev,il,jl,nat,nskip,iprint,k,
     .  jj1,jj2,jj3,is(3),ifac(3),nevl,iprm(nJq),nomega
      integer icomp,iib,jcomp,jib
      double precision wk(11*nJq),e(nJq),a,b,c,fac,abnrm,tolimw
      double precision qb(3,3),rb(3,3),qk,q1(3)
      double complex h(nJq,nJq),s(nJq,nJq),z(nJq,nJq),tmp,
     .  evl(nJq)
      double complex zsum
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./
      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr
      tolimw = 1d-4
      i1 = 9
      if (cmdopt('--tolimw=',i1,0,outs)) then
        if (.not. a2bin(outs,tolimw,4,0,' ',i1,-1)) call
     .    rxs2('EXASA: failed to parse "',outs(1:30),'%a"')
      endif

      call info0(20,1,0,' pvgfe8:  Make w(q) from J(q)')
C      call info(20,1,0,' pvgfe8: calculate omega from Jq'//
C     .  '%N  ib    moment',0,0)

C      dc = 'i'
C      if (list(1) == 0) dc = 'I'
C      call arrprt(' Site  amom','%,4i%:-3,4;4d',dc//'d',i,0,4,
C     .  0,' | ',amom,ez,w,w,w,w,w,w)

      lpos = .true.
      lneg = .true.

C ... Check whether all moments are positive or negative
      nat = nJq
      il = 0
      do  iib = 1, nsite
      ib = s_cpasite(iib)%ib
      do  icomp = 0, s_cpasite(iib)%ncomp
      if (icomp == 0 .and. s_cpasite(iib)%ncomp > 0) cycle
      il = il+1

        lpos = lpos .and. amom(ipc(ib)) > 0
        lneg = lneg .and. amom(ipc(ib)) < 0
C        print 333, ib, amom(ipc(ib))
C  333   format(i4,f12.6)

      enddo
      enddo
      if (il /= nat) call rx('pvgfe8: mismatch between nsite,nat')
      nomega = nat

      if (lpos .or. lneg) then
        if (iprint() >= 40) then
          call info0(20,0,0,' pvgfe8:  ferromagnetic case (mRy)'//
     .      '%N  i1  i2  i3%19fq%24fevl ...')
        else
          call info0(20,0,0,' pvgfe8:  ferromagnetic case ...')
        endif
      else if (nat == 2) then
        call info0(20,0,0,' pvgfe8:  antiferromagnetic case ...')
      elseif (iprint() >= 40) then
        call info0(20,0,0,' pvgfe8:  ferri case, real eigenvalues'//
     .    ' of J (mRy)'//
     .    '%N  i1  i2  i3%19fq%24fevl ...')
      else
        call info0(20,0,0,' pvgfe8:  ferri case ...')
      endif

C ... Make spin waves
      call dpzero(h,nJq*nJq*2)
C     If moments are uniformly >0, overlap is positive definite.
C     If moments are uniformly <0, J does not change by mom -> -mom ;
C     To address this case, use -1/mom for overlap instead of 1/mom
C     If neither applies, handle special case nat=2
      nskip = 0
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1

        call dpzero(s,nJq*nJq*2)

        il = 0
        do  iib = 1, nsite
        ib = s_cpasite(iib)%ib
        do  icomp = 0, s_cpasite(iib)%ncomp
        if (icomp == 0 .and. s_cpasite(iib)%ncomp > 0) cycle
        il = il+1

          jl = 0
          do  jib = 1, nsite
          jb = s_cpasite(jib)%ib
          do  jcomp = 0, s_cpasite(jib)%ncomp
          if (jcomp == 0 .and. s_cpasite(jib)%ncomp > 0) cycle
          jl = jl+1

C           Hamiltonian, in units of electron spin (Kotani convention)
            h(il,jl) = -2/amom(ipc(ib))*JqRR(i1,i2,i3,il,jl)
     .                 *2/amom(ipc(jb))

          enddo
          enddo

C         Subtract q=0 term
C          tmp = 0
C          jl = 0
C          do  jib = 1, nsite
C          jb = s_cpasite(jib)%ib
C          do  jcomp = 0, s_cpasite(jib)%ncomp
C          if (jcomp == 0 .and. s_cpasite(jib)%ncomp > 0) cycle
C          jl = jl+1
C            tmp = tmp + JqRR(1,1,1,jl,il)
C          enddo
C          enddo
          tmp = zsum(nat, JqRR(1,1,1,1,il), n1*n2*n3)
          h(il,il) = h(il,il) + (2/amom(ipc(ib)))**2 * tmp

C         Overlap is just 1/amom
          s(il,il) = 1/amom(ipc(ib))
          if (lneg) s(il,il) = -s(il,il)

        enddo
        enddo

C         if (i1 == 2 .and. i2 == 1 .and. i3 == 1) then
C           print *, 'hi',i1,i2,i3
C           call zprm('h',2,h,nJq,nat,nat)
C          print *, 'hi'
C         endif
C       call zprm('h',2,h,nJq,nJq,nJq)

C   ... Debugging check
        do  ib = 1, nat
          do  jb = 1, nat
            if (abs(h(ib,jb)-dconjg(h(jb,ib))) > 1d-8)
     .        call rx('pvgfe8: h is not symmetric')
          enddo
        enddo

C   ... FM system: eigenvalues guaranteed to be real
        if (lpos .or. lneg) then
C          call zprm('s',2,s,nJq,nat,nat)
C          if (nat == 1) then
C            e(1) = h(1,1)/s(1,1)
C          else
            call zhevx(nat,nJq,h,s,1,.true.,0,9d9,nev,wk,.false.,e,
     .        nJq,z)
C          endif
C         call prmx('e',e,nat,nat,1)
          call dcopy(nat,e,1,omega(i1,i2,i3,1),n1*n2*n3)
          nevl = nat
C   ... Simplest AFM system: explicitly check eigenvalues
        else if (nat == 2) then
C         call zprm('s',2,s,4,2,2)
C         call zprm('h',2,h,4,2,2)
          a = s(1,1)*s(2,2)
          b = -(h(2,2)*s(1,1) + h(1,1)*s(2,2))
          c = h(1,1)*h(2,2) - h(1,2)*h(2,1)
          fac = (b/a/2)**2 - c/a
C         Real solution does not exist
          if (fac < 0) then
            omega(i1,i2,i3,1) = 0
            if (abs(fac) > 1d-5) nskip = nskip + 1
          else
            omega(i1,i2,i3,1) = (-b/a/2) + sqrt(fac)
          endif
          nevl = -1
C   ... General case: make (nonhermitian) S^-1 H; find its eigenvalues
        else
C          call zprm('s',2,s,nJq,nat,nat)
C          call zprm('h',2,h,nJq,nat,nat)

          do  ib = 1, nat
            do  jb = 1, nat
              h(ib,jb) = h(ib,jb) / s(ib,ib)
            enddo
          enddo
          if (cmdopt('--wsw',5,0,outs)) then
            q1(1) = qk(1,i1,i2,i3)
            q1(2) = qk(2,i1,i2,i3)
            q1(3) = qk(3,i1,i2,i3)
            print *, 'writing Heisenberg secular matrix, for this q:'
                print 1,i1,i2,i3,q1
            call zprm('h',2,h,nJq,nat,nat)
          endif
          call zgeevs('P','N','V','N',nat,h,nJq,evl,evl,1,s,nJq,abnrm,
     .      il)
C         Look for real solutions, lowest positive and highest negative
          nevl = 0
          call dvset(e,1,nat,999d0)
          do  ib = 1, nat
            if (iprint() >= 50) then
            if (dble(evl(ib)) >= 0) print "(i4,2f12.6)", ib, evl(ib)
            endif
            if (dabs(dimag(evl(ib))/dble(evl(ib))) < tolimw .and.
     .          dble(evl(ib)) > 0d0 .or.
     .          abs(evl(ib)) < tolimw/10) then
              nevl = nevl+1
              e(nevl) = dble(evl(ib))
            endif
          enddo
          if (nevl > 0) then
            call dvheap(1,nevl,e,iprm,0d0,0)
          else
            nskip = nskip + 1
          endif
          nomega = min(nomega,nevl)
          call dcopy(nat,e,1,omega(i1,i2,i3,1),n1*n2*n3)

        endif

        if (iprint() >= 40 .and. nevl >= 0) then
          q1(1) = qk(1,i1,i2,i3)
          q1(2) = qk(2,i1,i2,i3)
          q1(3) = qk(3,i1,i2,i3)
              print 1,i1,i2,i3,q1,(e(ib),ib=1,nevl)
        endif

      enddo
      enddo
      enddo

      if (nskip > 0)
     .  call info2(20,0,0,
     .  ' pvgfe8:  Negative or imaginary values of omega for %i qp'//
     .  ' out of %i',nskip,n1*n2*n3)
      if (nomega /= nat)
     .  call info2(20,0,0,
     .  ' pvgfe8:  Largest number of positive SW frequencies found'//
     .  ' at every omega: %i',nomega,0)
    1 format(3i4,3f12.6,3x,3P,100(f10.3))


C     call prm3('lowest eval omega(q)',0,omega,n1,n2,n3)

      end

      subroutine pvgfe9(plat,jqin,pos,n1,n2,n3,jq)
C- Make J(q) for ib ne jb, shifting by phase
C  xxx this routine is wrong
C  cannot eliminate phase unless inversion symmetry
      implicit none
      integer n1,n2,n3
      double precision jq(n1,n2,n3),plat(3,3),pos(3)
      double complex jqin(n1,n2,n3)
C Local variables
      integer is(3),ifac(3),k,jj1,jj2,jj3,i3,i2,i1
      double precision qb(3,3),rb(3,3),qk,q1(3),twopi,sp
      double complex phase
      logical llshft(3)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./

      twopi = 8d0*datan(1d0)
      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr

      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

C   ... exp(i q pos)
        sp = twopi*(pos(1)*q1(1) + pos(2)*q1(2) + pos(3)*q1(3))
        phase = dcmplx(dcos(sp),-dsin(sp))*jqin(i1,i2,i3)
C   ... debugging
        if (dabs(dimag(phase)) > 1d-10) then
C         print *, 'bug in pvgfe9: phase=',phase
          return
        endif
        jq(i1,i2,i3) = dble(phase)
      enddo
      enddo
      enddo
      end

      subroutine pvgfea(alat,ng,gv,cJ,lambda,cJr2)
C- Given FT cJ, return cJr2 = cJ*g**2*exp(-lambda*g)
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors
Ci   ng    :number of group operations
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cJ    :coefficients to scale
Ci   lambda:factor to scale cJ
Co Outputs
Co   cJr2  :cJ * g**2 * exp(-lambda*g)
Cr Remarks
Cr
Cu Updates
Cu   05 Jan 03 Added extra factor exp(-lambda * g)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng
      double precision gv(ng,3),alat,lambda
      double complex cJ(ng),cJr2(ng)
C ... Local parameters
      integer i
      double precision g2

      if (lambda == 0) then
        do  i = 1, ng
          g2 = alat**2*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
          cJr2(i) = cJ(i)*g2
        enddo
      else
        do  i = 1, ng
          g2 = alat**2*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
          cJr2(i) = cJ(i) * g2 * exp(-lambda*sqrt(g2))
        enddo
      endif
      end

      subroutine pvgfeb(plat,jr,pos,n1,n2,n3,jq)
C- Fourier transform by brute-force
      implicit none
      integer n1,n2,n3
      double precision jr(n1,n2,n3),plat(3,3),pos(3)
      double complex jq(n1,n2,n3)
C Local variables
      integer is(3),ifac(3),k,jj1,jj2,jj3,i3,i2,i1,j1,j2,j3
      double precision qb(3,3),rb(3,3),qk,q1(3),twopi,sp,p(3)
      double complex phase
C     double complex jqi
      logical llshft(3)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./

      twopi = 8d0*datan(1d0)
      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr

      call dpzero(jq,2*n1*n2*n3)

C      print *, qb*plat
C      q1(1) = qk(1,n1,n2,n3)
C      q1(2) = qk(2,n1,n2,n3)
C      q1(3) = qk(3,n1,n2,n3)
C      j1 = 1
C      j2 = 1
C      j3 = 1
C      p(1) = j1*plat(1,1)+j2*plat(1,2)+j3*plat(1,3) - pos(1)*0
C      p(2) = j1*plat(2,1)+j2*plat(2,2)+j3*plat(2,3) - pos(2)*0
C      p(3) = j1*plat(3,1)+j2*plat(3,2)+j3*plat(3,3) - pos(3)*0
C      sp = twopi*(p(1)*q1(1) + p(2)*q1(2) + p(3)*q1(3))
C      print *, sp, qb(1,1:3)*p(1:3),qb(2,1:3)*p(1:3),qb(3,1:3)*p(1:3)

C     jqi = 0
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

        do  j1 = 0, n1-1
        do  j2 = 0, n2-1
        do  j3 = 0, n3-1

C     ... multiply pos*0 for FFT convention
C     ... multiply pos*1 for real jq
          p(1) = j1*plat(1,1)+j2*plat(1,2)+j3*plat(1,3) - pos(1)*1*0
          p(2) = j1*plat(2,1)+j2*plat(2,2)+j3*plat(2,3) - pos(2)*1*0
          p(3) = j1*plat(3,1)+j2*plat(3,2)+j3*plat(3,3) - pos(3)*1*0
          sp = twopi*(p(1)*q1(1) + p(2)*q1(2) + p(3)*q1(3))
          phase = dcmplx(dcos(sp),dsin(sp))
C          if (i1 == 2 .and. i2*i3 == 1) then
C          print 333, i1,i2,i3,j1,j2,j3, p,sp/twopi,jr(j1+1,j2+1,j3+1)
C    333   format(6i4,3f12.6,2x,2f12.6)
C          jqi = jqi + phase*jr(j1+1,j2+1,j3+1)
C          endif

          jq(i1,i2,i3) = jq(i1,i2,i3) + phase*jr(j1+1,j2+1,j3+1)
        enddo
        enddo
        enddo
C      call zfftqp(plat,dcmplx(jr),pos,n1,n2,n3,q1,1,jq(i1,i2,i3))

      enddo
      enddo
      enddo

C     print *, jqi
      end
      subroutine pvgfec(nq,ni,nj,ldi,ldj,nspc,ldg,gfbz,gii)
C- k-integrated g for a subblock
      implicit none
      integer nq,ni,nj,ldi,ldj,nspc,ldg
      double precision gfbz(2,nq,ldi,nspc,ldj,nspc),gii(ldg,nspc,ldg,nspc,2)
      integer i,j,iq,is,js

      do  is = 1, nspc
      do  i = 1, ni
      do  js = 1, nspc
      do  j = 1, nj
        gii(i,is,j,js,1) = 0d0
        gii(i,is,j,js,2) = 0d0
        do  iq = 1, nq
          gii(i,is,j,js,1) = gii(i,is,j,js,1) + gfbz(1,iq,i,is,j,js)
          gii(i,is,j,js,2) = gii(i,is,j,js,2) + gfbz(2,iq,i,is,j,js)
        enddo
      enddo
      enddo
      enddo
      enddo

      end

C      subroutine pvgfez(nq,ni,nj,ldi,ldj,nspc,ldg,gfbz,dmatq)
CC- Energy-integrated g for a subblock
C      implicit none
C      integer nq,ni,nj,ldi,ldj,nspc,ldg
C      double precision gfbz(2,nq,ldi,nspc,ldj,nspc),
C     .  dmatq(ldg,nspc,ldg,nspc,2)
C      integer i,j,iq,is,js
C
C      do  is = 1, nspc
C      do  i = 1, ni
C      do  js = 1, nspc
C      do  j = 1, nj
C        dmatq(i,is,j,js,1) = 0d0
C        dmatq(i,is,j,js,2) = 0d0
C        do  iq = 1, nq
C          dmatq(i,is,j,js,1) = dmatq(i,is,j,js,1) + gfbz(1,iq,i,is,j,js)
C          dmatq(i,is,j,js,2) = dmatq(i,is,j,js,2) + gfbz(2,iq,i,is,j,js)
C        enddo
C      enddo
C      enddo
C      enddo
C      enddo
C      end
C
      subroutine pvgfed(offH,s_cpasitei,nsite,s_site,mxcomp,sumRL,sumR)
C- Printout of J0
C ----------------------------------------------------------------------
Cio Structures
Cio  s_cpasitei :struct for site list with CPA components; see structures.h
Ci     Elts read:  ib,ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   nsite :number of sites in s_cpasitei
Ci   mxcomp:dimensions sumRL
Ci   sumRL :quantities needed for sum rules and on-site exchange
Ci         :interactions.  Real and imaginary parts retained.
Ci         :See pasajj in asajft.f
Co Outputs
Co    On-site exchange interactions are printed out
Cl Local variables
Cr Remarks
Cr   For formulas, see pasajj in asajft.f
Cr   L-resolved on-site exchanges could be further resolved into
Cr   (L,L') blocks.  Here, they are resolved by L but contracted by L'
Cu Updates
Cu   16 Jan 13 Reworked site list for CPA compatibility
Cu   05 Jun 01 first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ib,nsite,mxcomp,nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,*)
      double complex sumRL(11,0:mxcomp,0:mxcomp,*)
      double complex sumR(11,0:mxcomp,0:mxcomp,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_cpasite):: s_cpasitei(nsite)
C ... Local parameters
      integer i,nlma,nglob,offai,stdo,ipr,iib,icomp,jcomp
      double precision fac,facJ0,kboltz,JL(81),JL0(81),J0i(0:mxcomp)

      parameter (kboltz=8.617d-5)
C     logical lerr

      stdo = nglob('stdo')
      call getpr(ipr)
      if (ipr < 20) return
      fac = 1000
      facJ0 = 2d0/3d0 * 13.6d0 / kboltz

C --- Print site-resolved J0 ---
C      write(stdo,1) ' comp  '
C    1 format(/' Sum J_ij for all pairs connected to listed sites (mRy)'/
C     .  ' site',a,'Sum J_ij       J_00',9x,'J_0',4x,'2/3 J_0 (K)')
C
C      do  iib = 1, nsite
C      ib = s_cpasitei(iib)%ib
C      do  icomp = 0, s_cpasitei(iib)%ncomp
CC     if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle
C
CC   ... Printout site-resolved J0
C        i = 20; if (icomp > 0) i = 31
C        call info8(i,0,0,'%,4i%?#n#%-1j%,4i#    #'//
C     .    '%,3;12D%,3;12D%,3;12D%,1;10D',ib,icomp,
C     .      fac*dimag(sumR(1,icomp,0,ib)),
C     .      fac*dimag(sumR(3,icomp,0,ib)),
C     .      fac*dimag(sumR(1,icomp,0,ib)-sumR(3,icomp,0,ib)),
C     .      facJ0*dimag((sumR(1,icomp,0,ib)-sumR(3,icomp,0,ib))),0,0)
C      enddo
C      enddo

      write(stdo,2) ' comp  '
    2 format(/' Sum J_ij for all pairs connected to listed sites (mRy)'/
     .  ' site',a,'Sum J_ij',7x,'J_00',7x,'J_0(ij)',5x,
     .  'J_0(i)',2x,'2/3 J_0 (K)')


C --- Print total mean field parameters for each (site,component) ---
      do  iib = 1, nsite
      ib = s_cpasitei(iib)%ib
      nlma = offH(4,1,ib+1) - offH(4,1,ib)

C ... Extract one-site rotation J0i
C     Non-CPA site: extract J0i from sumR
      if (s_cpasitei(iib)%ncomp == 0) then
        J0i(0) = dimag(-sumR(2,0,0,ib)-sumR(3,0,0,ib))
C     CPA site: extract J0i from s_site(ib)%j0
      else
        J0i(0) = 0
        do  icomp = 1, s_cpasitei(iib)%ncomp
          J0i(icomp) = sum(s_site(ib)%j0(icomp,1:nlma))
          J0i(0) = J0i(0) + J0i(icomp)*s_site(ib)%cpawt(icomp)
        enddo
      endif

      do  icomp = 0, s_cpasitei(iib)%ncomp
C     if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle

C   ... Printout site-resolved J0
        i = 20; if (icomp > 0) i = 31
        call info8(i,0,0,'%,4i%?#n#%-1j%,4i#    #'//
     .    '%,3;12D%,3;12D%,3;12D%,3;12D%,1;10D',ib,icomp,
     .      fac*dimag(sumR(1,icomp,0,ib)),
     .      fac*dimag(sumR(3,icomp,0,ib)),
     .      fac*dimag(sumR(1,icomp,0,ib)-sumR(3,icomp,0,ib)),
     .      fac*J0i(icomp),facJ0*J0i(icomp),0)
      enddo
      enddo

C --- Print L-resolved J0 ---
      write(stdo,4)
    4 format(/' J_0 resolved by L (mRy), from Jij'/
     .  ' site,cmp L  1',6x,'2',6x,'3     ...')

      do  iib = 1, nsite
      ib = s_cpasitei(iib)%ib
      nlma = offH(4,1,ib+1) - offH(4,1,ib)
      offai = offH(4,1,ib)
      call dpzero(Jl0,nlma)
      do  icomp = 0, s_cpasitei(iib)%ncomp
        if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle

C       Contributions from non-CPA sites connected to icomp
        call dpzero(Jl,nlma)
        do  i = 1, nlma
          JL(i) = dimag(sumRL(1,icomp,0,i+offai))
     .          - dimag(sumRL(3,icomp,0,i+offai))
        enddo
C       Contributions from CPA sites connected to icomp
C       sumRL(:,;,jcomp,:) already weighted by cpawt(jcomp); see pasajf
        do  jcomp = 1, mxcomp
        do  i = 1, nlma
          JL(i) = JL(i) +
     .          (dimag(sumRL(1,icomp,jcomp,i+offai))
     .          -dimag(sumRL(3,icomp,jcomp,i+offai)))

        enddo
        enddo
        call dscal(nlma,fac,JL,1)
        if (icomp > 0) then
          call dpadd(JL0,JL,1,nlma,s_site(ib)%cpawt(icomp))
          write(stdo,'(i4,i3,3x,(16f7.3))') ib,icomp,(JL(i), i=1,nlma)
        else
          write(stdo,'(i4,3x,3x,(16f7.3))') ib, (JL(i), i=1,nlma)
        endif
      enddo
      if (s_cpasitei(iib)%ncomp /= 0) then
        write(stdo,'(i4,3x,3x,(16f7.3))') ib, (JL0(i), i=1,nlma)
      endif
      enddo

      if (mxcomp == 0) return
      write(stdo,5)
    5 format(/' J_0 resolved by L (mRy), from one-site rotation'/
     .  ' site,cmp L  1',6x,'2',6x,'3     ...')

      do  iib = 1, nsite
      ib = s_cpasitei(iib)%ib
      if (s_cpasitei(iib)%ncomp == 0) cycle
      nlma = offH(4,1,ib+1) - offH(4,1,ib)
      offai = offH(4,1,ib)
      call dpzero(Jl0,nlma)
      do  icomp = 1, s_cpasitei(iib)%ncomp
        do  i = 1, nlma
          JL(i) = s_site(ib)%j0(icomp,i)
        enddo
        call dscal(nlma,fac,JL,1)
        call dpadd(JL0,JL,1,nlma,s_site(ib)%cpawt(icomp))
        write(stdo,'(i4,i3,3x,(16f7.3))') ib,icomp,(JL(i), i=1,nlma)
      enddo
      write(stdo,'(i4,3x,3x,(16f7.3))') ib, (JL0(i), i=1,nlma)
      enddo

      end

      subroutine pvgfee(nbas,iblst,jblst,n1,n2,n3,lwrite,ipc,amom,Jsum,
     .  JqRR,ifi)
C- Subtract off on-site contribution to JqRR and optionally output
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   iblst :list of sites for which Jqrr defined in 4th (row) dimension
Ci   jblst :list of sites for which Jqrr defined in 5th (col) dimension
Ci   n1..n3:number of mesh points along the three lattice vectors
Ci   lwrite:0 write nothing
Ci         :1 write Jqrr(1,1,1) for mean-field hamiltonian
Ci         :Add 10 not to scale Jqrr by moments
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   amom  :sphere charges by class
Ci   Jsum  :Jsum(1,ib) is q=0 term
Ci         :Jsum(2,ib) is on-site contribution
Ci   ifi   :file handle
Co Outputs
Co  JqRR   :On-site contribution to JqRR is subtracted.
Co         :Parts of JqRR may be written to disk, depending on lwrite
Cu Updates
Cu   23 Oct 02 first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,nbas,iblst(nbas),jblst(nbas),ifi,lwrite,ipc(nbas)
      double precision Jsum(3,nbas),amom(nbas)
      double complex JqRR(n1,n2,n3,nbas,nbas)
C ... Local parameters
      integer ib,jb,i1,i2,i3,ilst,nlst
      double precision kboltz,fac,amoml(nbas)
C     kboltz is in units eV K^-1
      parameter (kboltz=8.617d-5)
      character*(20) fmt
C     double complex h(n1,n2,n3)

C     call fftz3(JqRR,n1,n2,n3,n1,n2,n3,1,0,-1)

      ilst = 0
      do  ib = 1, nbas
        amoml(ib) = 1
        if (lwrite >= 10) amoml(ib) = amom(ib)
        if (iblst(1) > 0) then
          if (iblst(ilst+1) /= ib) cycle
        endif
        ilst = ilst+1

        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          JqRR(i1,i2,i3,ib,ib) = JqRR(i1,i2,i3,ib,ib) - Jsum(2,ib)
        enddo
        enddo
        enddo
      enddo
      nlst = ilst
C     call fftz3(JqRR,n1,n2,n3,n1,n2,n3,1,0,-1)
      fmt = '(9f10.2)'
      fac = 13.6d0/kboltz

      if (lwrite == 1 .or. lwrite == 11) then
        rewind ifi
        if (jblst(1) == 0) then
        call awrit2('%% rows %i cols %i real',' ',80,ifi,nlst,nbas+1)
        else
        call awrit2('%% rows %i cols %i real',' ',80,ifi,nlst,nlst+1)
        endif
        ilst = 0
        do  ib = 1, nbas
          if (iblst(1) > 0) then
            if (iblst(ilst+1) /= ib) cycle
          endif
          ilst = ilst+1

          write(ifi,'(f12.6)') amoml(ipc(ib))
          if (jblst(1) == 0) then
            write(ifi,fmt)
     .      (dble(jqRR(1,1,1,ib,jb))*fac/amoml(ipc(ib))/amoml(ipc(jb)),
     .      jb=1,nbas)
          else
            write(ifi,fmt) (dble(jqRR(1,1,1,ib,jblst(jb)))*fac/
     .        amoml(ipc(ib))/amoml(ipc(jblst(jb))),
     .      jb=1,nlst)
          endif
        enddo
      endif

      end

      subroutine pvgfef(nk1,nk2,nk3,nbas,alat,omega,ng,gv,kv,wk)
C- Spin stiffness
      implicit none
C ... Passed parameters
      integer nk1,nk2,nk3,nbas,ng,kv(ng,3)
      double precision alat,omega(nk1,nk2,nk3,nbas),gv(ng,3)
      double complex wk(nk1,nk2,nk3)
C ... Local parameters
      integer stdo,i,ifi
      integer, parameter :: imax=10
      procedure(integer) :: iprint,nglob,fopna
      double precision facaa,D1,lambda,dsum,facl
      complex(8), pointer :: cJ(:),cJr(:)

      if (iprint() <= 10) return
      call info0(10,1,0,
     .  ' ... Stiffness from '//
     .  'lim(lambda->0) FT 1/6 T**2 omega(T) exp(-lambda T)')
      stdo = nglob('stdo')
      ifi = fopna('sstiff',-1,0)

C --- Make (complex) R.S. omega(T) ---
      call dpzero(wk,2*nk1*nk2*nk3)
      call dcopy(nk1*nk2*nk3,omega,1,wk,2)
C      do  i123 = 1, nk1*nk2*nk3
C        wk(i123,1,1) = omega(i123,1,1,1)
C      enddo
      call fftz3(wk,nk1,nk2,nk3,nk1,nk2,nk3,1,0,-1)
      allocate(cJ(ng))
      call gvgetf(ng,1,kv,nk1,nk2,nk3,wk,cJ)

C --- Sequence T**2 * omega(T) exp(-lambda T) for various lambda ---
      write(stdo,1)
    1 format(4x,'lambda    om(lambda,0)  om(lambda,0)  ln(om)'/
     .       4x,'          Ry-alat**2    meV-A**2     (file sstiff)')
      facl = .10d0
      allocate(cJr(ng))
      do  i = imax, 0, -1
        lambda = alat*facl*i
        call pvgfea(1d0,ng,gv,cJ,lambda,cJr)
C   ... Factor of 6 to convert nabla to coff of q**2
        D1 = -dsum(ng,cJr,2)/6

C       Debugging: h = R**2 * J on a uniform mesh
C       call gvputf(ng,1,kv,nk1,nk2,nk3,cJr,wk)
C       call zprm3('T**2 omega(T)',0,wk,nk1,nk2,nk3)
C       call fftz3(wk,nk1,nk2,nk3,nk1,nk2,nk3,1,0,1)
C       call zprm3('omega''''(q)',0,wk,nk1,nk2,nk3)
C       D1 = -dval(wk,1)/6

C       Scaling to convert from Ry-dimensionless to ev-AA^2
        facAA = 13.6d3 * (alat*.529d0)**2
        if (D1 > 0) then
          write(stdo,'(f11.6,f12.6,f11.3,f14.4)')
     .      lambda,D1*alat**2,D1*facAA,dlog(D1*facAA)
          if (i /= 0)
     .      write(ifi,'(f11.6,f12.6,f20.3)') lambda,dlog(D1*facAA)
        else
          write(stdo,'( '' D1 negative ... skip'')')
        endif
      enddo
      call fclose(ifi)
      deallocate(cJ,cJr)

      end

      subroutine pvgfeg(plat,ipq,n1,n2,n3,jq)
C- Make J(q) for ib ne jb, shifting by phase
C ----------------------------------------------------------------------
Ci Inputs
Ci   plat  :primitive lattice vectors, in units of alat
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   n1
Ci   n2
Ci   n3
Ci   jq
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   02 Aug 05 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,ipq(n1,n2,n3)
      double precision plat(3,3)
      double complex jq(n1,n2,n3)
C ... Local parameters
      integer is(3),ifac(3),k,jj1,jj2,jj3,i3,i2,i1,stdo,nglob
      double precision qb(3,3),rb(3,3),qk,q1(3)
      logical llshft(3)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./

      stdo = nglob('stdo')
C     twopi = 8d0*datan(1d0)
      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr

      write(stdo,356)
  356 format(/' EXASA: Jq-J(q=0)'/' i1..i3',25x,'qp',18x,'iq     J')

      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

        write(stdo,'(3i4,2x,3f12.6,i5,2f12.6)')
     .    i1,i2,i3,q1,ipq(i1,i2,i3),jq(i1,i2,i3)-jq(1,1,1)

      enddo
      enddo
      enddo
      stop
      end

      subroutine pvgfeh(plat,om,qp,n1,n2,n3,wq)
C- Make omega(q) for 1 qp from omega(T)
Cu Updates
Cu   20 Feb 07 First created
      implicit none
      integer n1,n2,n3
      double precision plat(3,3),qp(3)
      double complex om(n1,n2,n3),wq
C Local variables
      integer is(3),ifac(3),j1,j2,j3
      double precision qb(3,3),rb(3,3),twopi,sp,p(3)
      double complex phase
C     double complex jqi
      logical llshft(3)
C     double precision qk
C     integer jj1,jj2,jj3
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
C      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
C     .                    (jj2*ifac(2)-1)*qb(k,2) +
C     .                    (jj3*ifac(3)-1)*qb(k,3)
      data llshft /.false.,.false.,.false./

      twopi = 8d0*datan(1d0)
      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr

      wq = 0
      do  j1 = 0, n1-1
      do  j2 = 0, n2-1
      do  j3 = 0, n3-1
        p(1) = j1*plat(1,1)+j2*plat(1,2)+j3*plat(1,3)
        p(2) = j1*plat(2,1)+j2*plat(2,2)+j3*plat(2,3)
        p(3) = j1*plat(3,1)+j2*plat(3,2)+j3*plat(3,3)
        sp = twopi*(p(1)*qp(1) + p(2)*qp(2) + p(3)*qp(3))
        phase = dcmplx(dcos(sp),dsin(sp))
C        if (i1 == 2 .and. i2*i3 == 1) then
C        print 333, i1,i2,i3,j1,j2,j3, p,sp/twopi,om(j1+1,j2+1,j3+1)
C  333   format(6i4,3f12.6,2x,2f12.6)
C        jqi = jqi + phase*om(j1+1,j2+1,j3+1)
C        endif

        wq = wq + phase*om(j1+1,j2+1,j3+1)
      enddo
      enddo
      enddo

      end
      subroutine pvgfei(Jqij,Jqji,n1,n2,n3)
C- Render susceptibility symmetric
      implicit none
      integer n1,n2,n3
      double precision Jqij(n1,n2,n3),Jqji(n1,n2,n3)
      integer i1,i2,i3
      double precision xx

      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1

        xx = Jqij(i1,i2,i3)/2 + Jqji(i1,i2,i3)/2
        Jqij(i1,i2,i3) = xx
        Jqji(i1,i2,i3) = xx

      enddo
      enddo
      enddo

      end

      subroutine pvgfej(nbas,secmf,ez)
C- Mean-field estimate for Tc, multisite case
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :number of sites participating in MF Tc
Ci   secmf :Secular matrix for MF equations; see Remarks
Co Outputs
Co   ez    :eigenvector of largest eigenvalue (corresponding to Tc)
Co         :relative sizes of moments at Tc.
Cu Updates
Cu   09 Aug 10 First created
Cr Remarks
Cr   See:
Cr     P. W. Anderson, Solid State Physics 14 (1963)
Cr     E. Sasioglu, L. M. Sandratskii, and P. Bruno, Phys. Rev. B 70, 024427 (2004)
Cr     Hortamani, Sandratskii, Kratzer, Mertig, Scheffler, PRB78, 104402 (2008)
Cr     Largest eigenvalue of secmf gives the Curie Temperature
Cr   This routine assumes J's are defined wrt unit conventions for spin:
Cr       H = - Sum_i,j Jij ei.ej
Cr   Jii is by convention, 0
Cr   ei and ej are unit vectors for spin i and spin j.
Cr   Let ezi = <ei . zhat>: average z projection of unit vector
Cr   In MF theory the ezi obey a system of coupled equations:
Cr      ezi = 2/3 1/(kB Tc) sum_j J_ij ezj
Cr   Rearrange to get this eigenvalue problem:
Cr      (kB Tc - S) e = 0 where e is a vector of the ezj and
Cr      S = 2/3 sum_j J_ij
Cr   With one atom/cell this reduces to (Kudrnovsky PRB 64, 174402)
Cr      omega(q) = 4 / M [J(q) - J(0)]
Cr      kB Tc = 2/3 sum_j J_0j = M/6 1/N sum_q omega(q)
Cr   This is because:
Cr      J(q) = sum_T J_0T e^iq.T                             (1)
Cr      J(0) = sum_T J_0T                                    (2)
Cr   and
Cr      N^-1 sum_q e^iq.T = delta_T0                         (3)
Cr   So
Cr      N^-1 sum_q J(q) = sum_T J_0T sum_q e^iq.T = J_00     (4)
Cr   Combine (2) and (4)
Cr      sum_j J_0j = J(0) - N^-1 sum_q J(q) = N^-1 sum_q (J(0) - J(q))
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas
      double precision secmf(nbas,nbas),ez(nbas)
C ... Local parameters
      integer info,lwork
      double precision w(nbas),work(3*nbas),kboltz,facJ0
C     kboltz is in units eV K^-1
      parameter (kboltz=8.617d-5)

      if (nbas == 0) return
C     call prmx('sec mat for MF T',secmf,nbas,nbas,nbas)
      lwork = 3*nbas
      call dsyev('V','L',nbas,secmf,nbas,w,work,lwork,info)
      call rxx(info /= 0,'pvgfej: failed to find evals of secmf')

      facJ0 = 2d0/3d0 * 13.6d0 / kboltz
      call info2(20,0,0,'%N ... Mean-field estimate:  Tc = '//
     .  '%,3;3d meV =  %,1;1d K',w(nbas)*1000*13.6*2d0/3,w(nbas)*facJ0)

      call dcopy(nbas,secmf(1,nbas),1,ez,1)

      end

      subroutine pvgfek(n1,n2,n3,mixmod,s_cpasite,nsite,nJq,Jqij,ez)
C- Tc, Tyablikov formula, multi-atom case
C ----------------------------------------------------------------------
Cio Structures
Cio  s_cpasite :struct for site list with CPA components; see structures.h
Ci     Elts read:  ib,ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   n1..n3:mesh of q-points on which Jq is tabulated
Ci  mixmod :string containing iteration-dependent mixing.  (See parmxp)
Ci   nsite :number of sites in s_cpasite
Ci   nJq   :dimensions JqRR
Ci   Jqij  :exchange interactions on uniform mesh of q-points
Cio Inputs/Outputs
Co Outputs
Cio  ez    :On input, starting (trial) vector of moments for solution
Cio        :of generalized Tyablikov formula (see pvgfe4 for MF estimate)
Cio        :On output, ez corresponds to Tyablikov solution
Cr Remarks
Cr   See:
Cr     Rusz, Turek, Divis, PRB 71, 174408 (2005)
Cr     Hortamani et al, PRB 78, 104402 (2008)
Cr   This routine assumes J's are defined wrt unit conventions for spin:
Cr       H = - Sum_i,j Jij ei.ej
Cr   ei and ej are unit vectors for spin i and spin j.
Cr   Let ezi = <ei . zhat>: average z projection of unit vector
Cr
Cr   Parameters ezi must be determined through a set of coupled
Cr   nonlinear equations.
Cr   From the RTP paper, Eq. 7, define
Cr     Nij(q) = delta_ij [Delta_i + sum_k Jik(0) ezk] - ezi Jij(q)
Cr     Delta_i = 1/2 B field or (Hortamani) anisotropy energy,
Cr               Assumed 0 for now.
Cr   The RTP paper provides a set of equations for the ezi for any
Cr   given temperature.  In the limit T->Tc the equations simplify;
Cr   those are the equations adopted here.
Cr
Cr   As T->Tc the ezi all vanish.  Define a normalized vector
Cr   ez~ = (ez1,ez2,...)/L constructed so that ez~ has unit norm.
Cr   Renormalize Nij(q), substituting ezi -> ez~i; call it N~ij(q).
Cr   The set of simultanous equations defining Tc is
Cr     ezi = 2/3 1/(kB Tc) (1+1/Si) [1/N sum_q [N~ij(q)]^-1]^-1
Cr   This is a family of nJq equations in the ez~1..nJq
Cr   There is in addition a constraint (norm of ez~) and an additional
Cr   unknown Tc, which plays the role of scaling.
Cr
Cr   Interpolation to get converged integral
Cr   RTP use one polynomial interpolation scheme.  Used here is
Cr   a similar ideal but slightly different one; see pvgfe4.
Cr   In the MnPt test case both evaluate to about the same matrix.
Cr   The following mc commands will manually do the interpolations
Cr   given a matrix nij  (assumed to have dimension 2 in the mc test)
Cr   mc -veps=.0003 nij -s1e4 -a f0  f0 f0 -x -1:2 -a o o -s0,eps -+ -i f0 -x
Cr   mc -veps=.0001 nij -s1e4 -a f0 f0 -1:2 -a o o -s0,eps f0 -+ -i -s4 \
Cr      o -s0,eps+eps f0 -+ -i -- o -seps,eps f0 -+ -i -- o -s-eps,eps f0 -+ -i --
Cr   See Bugs for handling the q=0 special point.
Cb Bugs
Cr   These implementation underestimates the contribution from q=0.
Cr   Problem can be resolved with Offset-Gamma method, but not implemented yet.
Cr   Because of it, there is relatively slow convergence wrt k-mesh.
Cu Updates
Cu   16 Jan 13 Reworked site list for CPA compatibility
Cu   09 Aug 10 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer n1,n2,n3,nsite,nJq
      double complex Jqij(n1,n2,n3,nJq,nJq)
      double precision ez(nJq)
      character mixmod*(*)
C ... For structures
!      include 'structures.h'
      type(str_cpasite):: s_cpasite(nsite)
C ... Local parameters
C     logical lherm
      integer i,j,i1,i2,i3,ib,jb,il,jl,iprint,maxit,mmix
      double precision tol,eps,epsl,rms,avg,kboltz,ddot,w
      double complex tmp,delta,zsum,gam(4)
      double complex nij(nJq,nJq),nijbk(nJq,nJq),
     .  wk(nJq,nJq),wk2(nJq,nJq),nii(nJq,4)
C     kboltz is in units eV K^-1
      parameter (tol=1d-6,eps=1d-3,kboltz=8.617d-5,maxit=300,mmix=3)
      integer kpvt(mmix)
      double precision norm(mmix,mmix),wkmx(nJq,2+mmix,2),t(mmix)
      logical parmxp
C     for mixing
      integer broy,nmix,nkill,iter,npmix,amix
      double precision wt(3),rmsdel,rms2,wc,beta,betv,elind,xx
      character fnam*8
      integer iib,icomp,jib,jcomp

      call info2(20,1,0,' ... generalized Tyablikov estimate '//
     .  'for Tc, %i site(s)',nJq,0)

      beta = .3d0
      nmix = mmix

      if (mixmod /= ' ') then
        call dpzero(wt,3)
        rms2 = 0
C        bet2 = .3d0
        elind = 0
        xx = 0
C        locm = 9
        wc = 0
        betv = 0
        broy = 0
        if (.not. parmxp(1,mixmod,len(mixmod),broy,nmix,wt,beta,xx,
     .    elind,xx,xx,fnam,wc,nkill,betv,rms2))
     .    call rx('pvgfek: parse in parmxp failed')
      endif

C      call zprm('Jqij',2,Jqij,n1*n2*n3,n1*n2*n3,nJq**2)
C      call pvgfel(n1,n2,n3,nJq,ipc,amom,list,Jqij,lherm)
C      call zprm('Jqij',2,Jqij,n1*n2*n3,n1*n2*n3,nJq**2)
C      call pvgfel(n1,n2,n3,nJq,ipc,amom,list,Jqij,lherm)
C      call zprm('Jqij',2,Jqij,n1*n2*n3,n1*n2*n3,nJq**2)

C --- Next pass in iterative solution for Tc ---
      delta = 0
      call dpzero(nii,2*nJq*4)

      gam(1) = 0.00040138477872471689*0.5658305d0
      do  j = 2, 4
        gam(j) = j*gam(1)
      enddo

C      wkmx = 0 ; print *, '!!'

      npmix = 0
      call dcopy(nJq,ez,1,wkmx(1,1,2),1)
      do  iter = 1, maxit

C ... BZ integration of each (Nij)^-1 delta_ij
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1

        call dpzero(nij,nJq*nJq*2)
C   ... Build up N for this q
        il = 0
        do  iib = 1, nsite
        ib = s_cpasite(iib)%ib
        do  icomp = 0, s_cpasite(iib)%ncomp
          if (icomp == 0 .and. s_cpasite(iib)%ncomp > 0) cycle
          il = il+1

          jl = 0
          do  jib = 1, nsite
          jb = s_cpasite(jib)%ib
          do  jcomp = 0, s_cpasite(jib)%ncomp
            if (jcomp == 0 .and. s_cpasite(jib)%ncomp > 0) cycle
            jl = jl+1

           nij(il,jl) = - ez(il) * Jqij(i1,i2,i3,il,jl)

          enddo
          enddo

C     ... q=0 term adding to diagonal element
          tmp = 0
          jl = 0
          do  jib = 1, nsite
          jb = s_cpasite(jib)%ib
          do  jcomp = 0, s_cpasite(jib)%ncomp
          if (jcomp == 0 .and. s_cpasite(jib)%ncomp > 0) cycle

          jl = jl+1

            tmp = tmp + Jqij(1,1,1,il,jl) * ez(jl)
          enddo
          enddo

          nij(il,il) = nij(il,il) + (Delta + tmp)


        enddo
        enddo

C        print *, i1,i2,i3
C        call zprm('nij',2,nij,nJq,nJq,nJq)

C   ... Handle q=0 case using Eqn (x) in Remarks
        call zcopy(nJq**2,nij,1,nijbk,1)
        if (nJq > 1 .or. i1*i2*i3 /= 1) then
          do  j = 1, 4
            call zcopy(nJq**2,nijbk,1,nij,1)
            call zgemm('N','N',nJq,nJq,nJq,(1d0,0d0),nij,nJq,nij,
     .        nJq,(0d0,0d0),wk,nJq)
            tmp = zsum(nJq,nij,nJq+1)
            epsl = eps*dble(tmp*dconjg(tmp))/nJq + 1d-16
            epsl = 0.0004013847787247168d0**2
            do  i = 1, nJq
              wk(i,i) = wk(i,i) + gam(j)**2
            enddo
C           call zprm('wk',2,wk,nJq,nJq,nJq)
            call zqinvb('l',wk,nJq,nJq,nJq,wk2,nJq,wk2,nij,nJq,i)
            call rxx(i /= 0,'pvgfek failed to invert nij')
            do  jl = 1, nJq
              nii(jl,j) = nii(jl,j) + nij(jl,jl)
            enddo
          enddo
        endif
      enddo
      enddo
      enddo
      call dscal(2*4*nJq,1/dble(n1*n2*n3),nii,1)

C ... Scale to get Tc's from each independent equation
      rms = 0
      avg = 0
      do  jl = 1, nJq
        nii(jl,1) = dble(2*nii(jl,1) - nii(jl,2))
        nii(jl,1) = 2d0/3d0/nii(jl,1)/ez(jl)
        rms = rms + dble(nii(jl,1))**2
        avg = avg + dble(nii(jl,1))
      enddo
      rms = rms/nJq
      avg = avg/nJq
      rms = dsqrt(max((rms-avg*avg),0d0)/nJq)

C ... Convergence ... hooray!
      if (abs(rms/avg) < tol) then
        call info5(20,0,0,' ... RPA Tc converged in %i iter:'//
     .    '  Tc = '// '%,3;3d meV =  %,1;1d K',iter,avg*13.6*1000,
     .    avg*13.6/kboltz,0,0)
        exit
      endif

      if (iter == maxit) then
        call info2(20,0,0,' ... RPA Tc did NOT converge.  Current '//
     .    'value:  Tc = '// '%,3;3d meV =  %,1;1d K',avg*13.6*1000,
     .    avg*13.6/kboltz)
        exit
      endif

      call info5(35,0,0,'... iter %i : RMS deviation'//
     .  ' in Tc = %,3;3d meV   avg = %,3;3d meV',iter,
     .  rms*13.6*1000,avg*13.6*1000,0,0)
C     Generate new trial eigenvector from output nii
      do  jl = 1, nJq
        ez(jl) = ez(jl) * dble(nii(jl,1))
      enddo
      rms = dsqrt(ddot(nJq,ez,1,ez,1))
      call dscal(nJq,1/rms,ez,1)

C     Anderson mixing of ez
      call dcopy(nJq,ez,1,wkmx,1)
C     call pshpr(45)
      i = amix(nJq,npmix,mmix,0,beta,iprint()-15,10d0,norm,kpvt,
     .  wkmx,t,rmsdel)
C     call poppr
      call dcopy(nJq,wkmx(1,1,2),1,ez,1)
      npmix = min(npmix+1,mmix,nmix)
C     call info0(20,0,0,' ')
      enddo

C     call pvgfel(n1,n2,n3,nbas,ipc,amom,list,Jqij,lherm)
C     call zprm('Jqij',2,Jqij,n1*n2*n3,n1*n2*n3,nbas**2)

      call arrprt(' Site  ez','%,4i%:-3,4;4d','id',nJq,0,4,
     .  0,' | ',s_cpasite(:)%ib,ez,w,w,w,w,w,w)

      end

      subroutine pvgfel(n1,n2,n3,nbas,ipc,amom,list,Jqij,lherm)
C- Scale rows and columns of Jqij by sign of moment
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1..n3:mesh of q-points on which Jq is tabulated
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   amom  :magnetic moments, by class.  Spin is amom/2
Ci   list  :(list(1)>0): list of sites
Ci         :(list(1)=0): all sites are used.
Ci   Jqij  :exchange interactions on uniform mesh of q-points
Co Outputs
Co   Jqij  :scaled by sign(amom(i)) * sign(amom(j))
Co   lherm :T if all amom in list are of the same sign.
Cu Updates
Cu   09 Aug 10 First created
Cr Remarks
Cr   Input Jqij is expected to be hermitian, with real eigenvalues
Cr   If all moments are of the same sign Jqij will remain hermitian;
Cr   but not so if moments have mixed signs.
Cr   Two successive calls to this routine restors Jqij to its original state
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,nbas,ipc(nbas),list(nbas)
      double precision amom(nbas)
      double complex Jqij(n1,n2,n3,nbas,nbas)
      logical lherm
C ... Local parameters
      logical lpos,lneg
      integer i1,i2,i3,ib,jb,il,jl

      lpos = .true.
      lneg = .true.

      il = 0
      do  ib = 1, nbas
        if (list(1) > 0) then
          if (list(il+1) /= ib) cycle
        endif
        il = il+1
        lpos = lpos .and. amom(ipc(ib)) > 0
        lneg = lneg .and. amom(ipc(ib)) < 0
        if (amom(ipc(ib)) < 0) then
          jl = 0
          do  jb = 1, nbas
            if (list(1) > 0) then
              if (list(jl+1) /= jb) cycle
            endif
            jl = jl+1
C            if (ib /= jb) then
              do  i3 = 1, n3
              do  i2 = 1, n2
              do  i1 = 1, n1
C               Jqij(i1,i2,i3,ib,jb) = -Jqij(i1,i2,i3,ib,jb)
                Jqij(i1,i2,i3,jb,ib) = -Jqij(i1,i2,i3,jb,ib)
              enddo
              enddo
              enddo
C            endif
          enddo
        endif
      enddo

      lherm = lpos .or. lneg

      end

      subroutine pvgfem(opt,plat,n1,n2,n3,gausw,sig,Jq,nj)
C- Convolve J(q) with a smoothed delta function (discretized Gaussian)
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0 do nothing
Ci         :1 convolve with Gaussian
Ci   plat  :primitive lattice vectors, in units of alat
Ci  n1..n3:number of k-divisions, 3 lattice vectors
Ci   gausw :discretized Gaussian width
Ci   sig   :Defines direction of FT.
Ci          Use isig=-1 for usual definition of 'forward'
Ci             f~ = sum_x f(x) exp(-i q x) / (n1*n2*n3)
Ci          Use isig=1 for reverse definition:
Ci             f~ = sum_x f(x) exp( i q x)
Ci   Jq    :Jq(:,:,:,1..nj) are nj functions to convolve
Ci   nj    :number of functions Jq to convolve
Co Outputs
Co   Jq    :are convolved with Gaussians
Cr Remarks
Cr   Gaussians are normalized to make discrete sum unity
Cu Updates
Cu   11 Aug 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,opt,sig,nj
      double precision plat(3,3),gausw
      double complex Jq(n1,n2,n3,nj)
C ... Local parameters
      integer is(3),ifac(3),k,jj1,jj2,jj3,i3,i2,i1
      double precision qb(3,3),rb(3,3),qk,q1(3),q,ddot,x,fac
      double precision qlat(9),vol0,d
      complex(8),allocatable:: G(:,:,:)
      double complex zsum
      logical llshft(3)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

      data llshft /.false.,.false.,.false./

      if (opt == 0) return

      call fftz30(n1,n2,n3,i1,i2,i3)
      if (i1 /= n1 .or. i2 /= n2 .or. i3 /= n3)
     .  call rx('pvgfem needs different implementation of FFT')

      call pshpr(1)
      call bzmsh0(plat,llshft,0,n1,n2,n3,is,ifac,rb,qb)
      call poppr
      call mkqlat(plat,qlat,vol0)
      allocate(G(n1,n2,n3))

C ... Make the Gaussian on a discrete mesh
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        d = 1
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)
        call shorbz(q1,q1,qlat,plat)
        q = dsqrt(ddot(3,q1,1,q1,1))
        x = q/gausw
C       d = sm-delta-function; s is integral of d
        d = 0
        if (x < 36) then
          d = dexp(-x*x)
        endif
C       s = .5d0 * derfc(x)
        G(i1,i2,i3) = d

C       G(i1,i2,i3) = fac
C       G(i1,i2,i3) = q

      enddo
      enddo
      enddo

C ... Normalize so integral is unity
      fac = zsum(n1*n2*n3,G,1)/(n1*n2*n3)
      call dscal(n1*n2*n3*2,1/fac,G,1)

C ... Overwrite Jq with convolution of Jq and G
      do  i1 = 1, nj
        call fftz3c(Jq(1,1,1,i1),G,n1,n2,n3,n1,n2,n3,0,sig)
      enddo
      deallocate(G)

C     call zprm3('Jq(q) after convolution',0,Jq,n1,n2,n3)

      end

      subroutine pvgfen(opt,n1,n2,n3,nJq,plat,Jqij,nb1,nb2,nb3,Jqbig)
C- Map Jqij to a doubled mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :not used now.  In future:
Ci         :1 reduce to nlst;
Ci         :2 make option to double, or not.
Ci  n1..n3:number of k-divisions, 3 lattice vectors
Ci   nJq   :number of sites containing Jqij
Ci   plat  :primitive lattice vectors, in units of alat
Ci   Jqij  :exchange interactions on uniform mesh of q-points
Co Outputs
Co nb1..nb3:mesh for Jqbig
Ci   Jqbig :Jqij on a doubled mesh
Cu Updates
Cu   09 Aug 10 First created
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,n1,n2,n3,nJq
      double precision plat(3,3)
      double complex Jqij(n1,n2,n3,nJq,nJq),
     .  Jqbig(2*n1,2*n2,2*n3,nJq,nJq)
C ... Local parameters
      integer ib,jb,nb1,nb2,nb3,ipass
      double precision qlat(3,3),vol0
      if (opt == 0) return

      call info0(30,1,0,' ... generating J on a doubled mesh')

      call mkqlat(plat,qlat,vol0)
      nb1 = 2*n1
      nb2 = 2*n2
      nb3 = 2*n3

C     call zprm('jq',2,jqij,n1*n2,n1*n2,n3)
      ipass = 0
      do  ib = 1, nJq
      do  jb = 1, nJq
        call chgmsh(13,1,qlat,n1,n2,n3,n1,n2,n3,Jqij(1,1,1,ib,jb),
     .    qlat,nb1,nb2,nb3,nb1,nb2,nb3,Jqbig(1,1,1,ib,jb))
        if (ipass == 0) call pshpr(1)
        ipass = 1
      enddo
      enddo
      if (ipass == 1) call poppr

      end

      subroutine pvgfeo(n1,n2,n3,s_cpasite,nsite,tol,nJq,JqRR)
C- Symmetrize JqRR
C ----------------------------------------------------------------------
Cio Structures
Cio  s_cpasite :struct for site list with CPA components; see structures.h
Ci     Elts read:  ib,ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   n1..n3:mesh of q-points on which Jq is tabulated
Ci   nsite :number of sites in s_cpasite
Ci   nJq   :dimensions JqRR
Cio Inputs/Outputs
Cio  JqRR  :exchange interactions on uniform mesh of q-points
Cio        :is symmetrized
Cu Updates
Cu   16 Jan 13 Reworked site list for CPA compatibility
Cu   11 Aug 08 Spin waves for Ferri case (first attempt)
Cu   20 Feb 07 Spin waves for 2 atom AFM case, multiatom FM case
Cr Remarks
Cr   Force JqRR(ib,jb) = JqRR*(jb,ib) for each qp.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer n1,n2,n3,nsite,nJq
      double complex JqRR(n1,n2,n3,nJq,nJq)
      double precision tol
C ... For structures
!      include 'structures.h'
      type(str_cpasite):: s_cpasite(nsite)
C ... Local parameters
      integer i1,i2,i3,ib,jb,il,jl,icmp,iib,icomp,jcmp,jib,jcomp
      double precision errmx
      double complex zz

      errmx = 0
      il = 0

      icmp = 0
      do  iib = 1, nsite
      ib = s_cpasite(iib)%ib
      do  icomp = 0, s_cpasite(iib)%ncomp
      if (icomp == 0 .and. s_cpasite(iib)%ncomp > 0) cycle
      icmp = icmp+1
      il = il+1

        jcmp = 0
        jl = 0
        do  jib = 1, nsite
        jb = s_cpasite(jib)%ib
        do  jcomp = 0, s_cpasite(jib)%ncomp
        if (jcomp == 0 .and. s_cpasite(jib)%ncomp > 0) cycle
        jcmp = jcmp+1
        jl = jl+1

          do  i3 = 1, n3
          do  i2 = 1, n2
          do  i1 = 1, n1

            errmx = max(errmx,
     .        abs(JqRR(i1,i2,i3,icmp,jcmp) -
     .            dconjg(JqRR(i1,i2,i3,jcmp,icmp))))

            zz = (JqRR(i1,i2,i3,icmp,jcmp) +
     .            dconjg(JqRR(i1,i2,i3,jcmp,icmp)))/2
            JqRR(i1,i2,i3,icmp,jcmp) = zz
            JqRR(i1,i2,i3,jcmp,icmp) = dconjg(zz)
          enddo
          enddo
          enddo

        enddo
        enddo

      enddo
      enddo

      call info2(20,1,0,
     .  ' pvgfeo:  Symmetrize J(q) ... max deviation = %;3g',errmx,0)

      if (errmx > tol)
     .  call rx1('Deviation in J(q) outside specified tol = %e',tol)

      end

      subroutine pasajq(s_ctrl,s_pot,s_site,irep,nl,ib,lrel,nbas,iprmb,
     .  offH,nq,zp,wz,pp,vRLshf,vshft,ipc,nrc,nspc,gr,nlma,ldimb,
     .  ghh,ldh,dpf,ddpf,pfdim,qnu)
C- Kernel to integrate qnu from gr for site ib, one energy point
C ----------------------------------------------------------------------
Ci Inputs
Ci   irep  :tells pasajq about the contents of gr
Ci         :0 G is g=(P-S)^-1
Ci         :1 G is proper Green's function matrix
Ci   nl    :(global maximum l) + 1
Ci   ib    :site index
Ci   nbas  :size of basis
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   nq    :number of k-points on a single symmetry line for bands plotting
Ci   zp    :complex energy
Ci   wz    :weights for complex energy integration
Ci   pp    :potential parameters (atomsr.f)
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   vshft :array of site potential shifts
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nrc   :nrclas(i) = number of atoms in class i
Ci   gr    :Green's function connecting site ib to other sites
Ci   nlma  :dimensions gr
Ci   ldimb :dimensions gr
Ci   ghh   :GF for higher bloch
Ci   ldh   :dimension ghh
Ci   dpf   :energy derivative of potential function
Ci   ddpf  :2nd energy derivative of potential function
Co   pfdim :dimension dpf,ddpf
Co Outputs
Cr   qnu   :energy-weighted moments of the sphere charges added
Cu   qnur  :relativistic ASA energy moments
Cu Remarks
Cb Bugs
Cb   Doesn't make the correct qnur yet
Cu Updates
Cu  02 Sep 15 (Vishina) extended to relativistic case, first attempt
Cu  18 Jan 13 Extended to CPA case
Cu  16 Nov 07 Added orbital-dependent potential shifts (vRLshf)
Cu  15 Mar 00  extended to include iwaves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, parameter :: nkap0=4,n0H=5,nsp=2
      integer irep,nl,ib,nbas,nq,nlma,ldimb,nspc,ldh,pfdim,lrel,ispc
      integer iprmb(1),offH(n0H,nkap0,1),ipc(1),nrc(*)
      double precision zp(2),wz(2),vshft(-7:1),vRLshf(*)
      double precision pp(6,nl,*),qnu(3,nl,nsp,*),ghh(ldh)
C     double precision qnur(4,0:nl-1,2*nl,2,2,*),qnur1(4,0:nl-1,2*nl,2,2,*)
      double complex dpf(pfdim,nsp),ddpf(pfdim,nsp),gr(nq,nlma,nspc,ldimb,nsp)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      real(8), allocatable :: dvldu(:)
      integer,allocatable:: ncomp(:)
      double precision,allocatable:: gii(:,:,:,:,:)
C ... Local parameters
      logical lso,bittst
      integer isp,nlml,offl,nlmi,offi,npos,i,mxorb,ic,offpi
      integer indxi(ldimb),lidim,hdim
      real(8), target :: xx(16),xx1(16)
C     double precision gfr(4,4,2,2,nl,2*nl,4),q12
      double precision orbtm(nl,nspc,1)
C     For the relativistic case
      integer nclspd,nclasp,nangl,ldlm,nclass,nclassd,nlspcr,offd,offpd
C     double precision q21
      procedure(integer) :: nglob,iprint

C     double complex gcomp(nlmx,nspc,nlmx,nspc),gkmu(nlmx,nspc,nlmx,nspc)

C     call tcn('pasajq')

      xx = 0
      xx1 = 0

C     Permutation index for this site
C     call dpzero(gfr,4*4*2*2*nl*2*nl*4)

      allocate(gii(ldimb,nspc,ldimb,nspc,2)) ! Real, imag separated in gii
      call dpzero(gii,size(gii))
      call iprmsb(1,ib,ib,nbas,offH,iprmb,indxi,[xx],xx)
      mxorb = nglob('mxorb')

      offl  = offH(1,1,ib)
      offi  = offH(2,1,ib) + offH(1,1,nbas+1)
      nlml  = offH(1,1,ib+1) - offl ! number of orbitals in lower block
      nlmi  = offH(2,1,ib+1) - offH(2,1,ib)
      hdim  = offH(3,1,ib+1) - offH(3,1,ib)
      lidim = nlml + nlmi
      lso = bittst(s_ctrl%lncol,4)
      if (nlmi == 0) offi = 0     ! No intmd block. Avoid violating array bounds.

C ... gfidos needs the entire ncomp array ... actually not (check)
      allocate(ncomp(nbas))
      call sitepack(s_site,1,nbas,'ncomp',1,ncomp,xx)
      allocate(dvldu(lidim*nsp)); call dpzero(dvldu,lidim*nsp)

      do  isp = 1, nsp

C       Make gii: sum over k g_ij, lower block site diagonal, poke into gii(1,1)
        call pvgfec(nq,nlml,nlml,nlma,ldimb,nspc,ldimb,gr(1,1,1,1+offl,isp),gii)
C       Ditto, intermediate block site diagonal; poke into gii(1+nlml,1+nlml)
        if (nlmi /= 0)
     .  call pvgfec(nq,nlmi,nlmi,nlma,ldimb,nspc,ldimb,gr(1,1+nlml,1,1+offi,isp),gii(1+nlml,1,1+nlml,1,1))
C       Non-CPA sites need to be saved in SO case: s_site%gcu is used in place of gii
        if (lso .and. s_site(ib)%ncomp == 1) then
          call dscal(size(gii),1d0/nq,gii,1) ! Scale
          call ztoyy(gii,ldimb*nspc,ldimb*nspc,ldimb*nspc,ldimb*nspc,0,1)
C         Bug if CPA procedure does not order orbitals in gii by l.  CPA has no downfolding now.
          call pokeg0(1,s_site(ib)%norb,[0],ldimb,ldimb,ldimb,nspc,2,s_site(ib)%gcu,gii)
C          print *, 'after pokeg0',ib,sum(s_site(ib)%gcu)
C          call zprm('gcu after pokeg0',2,s_site(ib)%gcu,s_site(ib)%norb*2,s_site(ib)%norb*2,s_site(ib)%norb*2)
        endif

C   ... Convert g to proper Green's function
C       Onsite g taken from gfg2gr in the SO and lrel=2 cases
        if (irep == 0) then
          call pshpr(iprint()+10)
          i = 0
          if (lso) i = 2
          if (lso .and. bittst(s_ctrl%lbxc,4)) i = 12
          if (lrel == 2) i = 1
          if (i /= 0) then
            call gfg2gr(i,s_site(ib),1) ! Returns G in s_site(ib)%gc,1)
          else

C           Spin diagonal blocks only for charge density
            call dscal(2*pfdim*nsp,dble(nq),ddpf,1)
            do  ispc = 1, nspc
              offpi = (ispc-1)*pfdim
              call gfg2g2(0,ldimb*nspc,(ldimb*nspc)**2,nlml,nlml,offpi,offpi,
     .          dpf(1+offl,isp),ddpf(1+offl,isp),gii(1,ispc,1,ispc,1),npos)
              if (nlmi > 0) then
                call gfg2g2(0,ldimb*nspc,(ldimb*nspc)**2,nlmi,nlmi,offpi,offpi,
     .            dpf(1+offi,isp),ddpf(1+offi,isp),gii(1+nlml,ispc,1+nlml,ispc,1),npos)
              endif
            enddo
            call dscal(pfdim*nsp*2,1/dble(nq),ddpf,1)

C            call yprm0('(1p,9e18.10)')
C            call yprm('gii after scaling',2,gii,(ldimb*nspc)**2,(ldimb*nspc),(ldimb*nspc),(ldimb*nspc))

          endif
          call poppr
        endif

C       if (irep == 0 .and. lrel == 1 .and. .not. lso) then

C        if (irep == 0) then
C          call dscal(pfdim*nsp*2,dble(nq),ddpf,1)
C          call gfg2g2(0,ldimb,ldimb**2,nlml,nlml,0,0,
C     .      dpf(1+offl,isp),ddpf(1+offl,isp),gii,npos)
C          call gfg2g2(0,ldimb,ldimb**2,nlmi,nlmi,0,0,
C     .      dpf(1+offi,isp),ddpf(1+offi,isp),gii(1+nlml,1,1+nlml,1,1),npos)
C          call dscal(pfdim*nsp*2,1/dble(nq),ddpf,1)
C        elseif (irep == 0) then
C          i = 0
C          if (lrel == 2) i = 1
C          call gfg2gr(i,s_site,nbas)
C        endif

C       gii = gii/nq

C        call yprm0('(1p,9e18.10)')
C        call yprm('gii',2,gii,(ldimb*nspc)**2,(ldimb*nspc),(ldimb*nspc),(ldimb*nspc))

C       Accumulate qnu for class associated with this site
C       gii is available for only one "site", with tailored permutation index
C       skip the DLM sites, to be processed later (add 10000 to mode)

        if (.not. lso) then
          call gfidos(nl,1,10100,0,1,1,indxi,0,lidim,0,hdim,ldimb,gii,ldh,
     .      ghh,1d0/nq,zp,wz,isp,nsp,nspc,pp,vRLshf,vshft(1),ipc(ib),nrc,
     .      s_site(ib)%ncomp,nl,1,1,1,1,1,1,0,0,xx,xx,xx,qnu,xx,xx,xx,xx,xx)
        endif

C       qnu for DLM sites
        if (s_site(ib)%ncomp > 1 .or. lso) then
          i = 20100; if (lso) i = 100
C         i = i+400; call dpzero(orbtm,size(orbtm))
          if (.not. associated(s_pot%qcorr)) s_pot%qcorr = > xx
          call gfidos2(nl,1,i,0,1,1,indxi,0,lidim,0,hdim,ldimb,ldh,
     .      ghh,1d0/s_site(ib)%ncomp,zp,wz,isp,nsp,nspc,pp,vRLshf,vshft(1),
     .      ipc(ib),nrc,1,1,1,1,1,1,0,0,xx,xx,xx,xx,qnu,xx,xx,orbtm,xx,xx,
     .      s_site(ib)%ncomp,
     .      s_site(ib)%thet,
     .      s_site(ib)%dlmcl,
     .      s_site(ib)%cpawt,
     .      s_site(ib)%gc,
     .      s_site(ib)%gcorr,
     .      s_pot%mxy,
     .      s_pot%qcorr,
     .      s_site(ib)%norb,
     .      s_pot%sop)
        endif

C       qnu for fully relativistic case
        if (lrel == 2) then
!          call rx('check lrel=2 branch')
C          ncomp1 = s_site(ib)%ncomp ; norb = s_site(ib)%norb ! Check the numbers, dim of...
C          n2 = 2*norb
C          gcomp(:,:,:,:) = dcmplx(gii(:,:,:,:,1),gii(:,:,:,:,2)) !Check if works
C          call mstokm(0,nl,ncomp1,norb,gcomp,gkmu,idum)
C          gcomp = gkmu
C          call zcopy(norb*norb*4*ncomp1,gcomp,1,s_site(ib)%gc,1)
          ldlm = s_ctrl%ldlm
          nangl = s_ctrl%nccomp
          nclasp = s_ctrl%nclasp
          nclass = s_ctrl%nclass
          if (ldlm /= 0 .and. nsp /= 2)
     .    call rx('PASAJQ Bug: CPA needs NSPIN=2')
          if (ldlm == 0) nangl = 0
          nclspd = nclasp + nangl
          nclassd = nclass + nangl
          nlspcr = 8*nl*nl*max(nclassd,s_ctrl%nspec)
          ic = s_ctrl%ipc(ib)
          offd = 0 ; offpd = 0
          ic = s_ctrl%ipc(ib)
          call gfidosr(nl,nbas,102,0,ib,iprmb,0,lidim,0,
     .      hdim,lidim,hdim,ghh,1d0,zp,wz,isp,nsp,nspc,s_pot%pprel,
     .      dvldu,vshft(1),s_ctrl%ipc,s_ctrl%nrc,1,1,1,1,1,1,
     .      offd,offpd,xx1,xx1,s_site(ib)%pdos,xx1,s_pot%qnur(1,ic),xx,
     .      orbtm,s_site(ib)%dmat,xx1,ncomp(ib),s_site(ib)%thet,
     .      s_site(ib)%dlmcl,s_site(ib)%cpawt,s_site(ib)%gc,
     .      s_site(ib)%gcorr,s_pot%mxy,s_pot%qcorr,s_site(ib)%norb,
     .      s_pot%GFr(1,ic))

C            call gfidosr(100,s_site,nl,ic,1,s_ctrl%nrc(ic),zp,wz,s_pot%pprel,vshft(1),s_ctrl%ipc,
C     .           s_ctrl%nrc,s_pot%gfr(1,ic),1,1,1,1,xx1,s_pot%qnur(1,ic),q21)
C          call qnu2qnur(21,nl,nsp,qnu(1,1,1,ic),qnur1(1,1,1,1,1,ic))
C          call qnu2qnur(11,nl,nsp,qnu1,qnur)
        endif
        if (nspc == 2) exit
      enddo
      deallocate(gii)
C     call tcx('pasajq')
C      ic = 1
C      call qnu2qnur(21,nl,nsp,qnu(1,1,1,ic),qnur1(1,1,1,1,1,ic))
      end
