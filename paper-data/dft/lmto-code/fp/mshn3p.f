      subroutine mshn3p(nbas,s_site,s_spec,lmet,lrout,lfrce,
     .  qval,ef0,def,sumqv,sumev,n1,n2,n3,k1,k2,k3,smrho,f,lrep)
C- Interpolate the density on the three-point energy mesh for metals
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qkkl qhkl qhhl
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rsma lmxl lmxb lmxa kmxt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   lmet  :0 semiconductor--do nothing
Ci         :>0 interpolate using 3-point sampling scheme; see Remarks
Ci   lrout :0 interpolate ef0,sumqv,sumev only
Ci         :1 also interpolate smrho,oqkkl,f
Ci   lfrce :>0 interpolate forces
Ci   qval  :total valence charge
Ci   def   :Fermi level window.
Ci   k1,k2,k3:dimensions smrho,smpot
Cio Inputs/Outputs
Ci   ef0   :On input: trial  Fermi level; see Remarks
Ci         :On output: interpolated Fermi level
Ci   sumqv :on input, total charge at each of the trial Fermi levels
Ci         :second channel is magnetic moment
Ci         :out output, interpolated charge and magnetic moment
Ci   sumev :on input, eigenvalue sum at each trial Fermi level
Ci         :out output, interpolated eigenvalue sum
Co Outputs
Co   smrho :smooth density on uniform mesh interpolated to true Fermi level
Co   f     :forces interpolated to true Fermi level
Co   lrep  :set to nonzero if Fermi level outside acceptable bounds
Cr Remarks
Cr   Input quantities smrho,qkkl,f were accumulated input for three
Cr   separate fermi levels, ef0, ef0+def, ef0-def.  The true Fermi
Cr   level is calculated, and a quadratic interpolation is made
Cr   for these three quantities.
Cr
Cr   Interpolated data is put back into first slot (iq=1).
Cr   If Ef is crazy and the charges are unusable, set flag lrep and return.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   02 Jan 06 sumqv resolved by spin
Cu   21 Mar 01 Bug fix for case q1,q3 close to qval
Cu   23 Jan 01 Added lrout switch
Cu   17 Jun 00 spin polarized
Cu   22 May 00 Adapted from nfp ip_density.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,lmet,lfrce,lrout,lrep,nbas
      double precision def,ef0,qval
      double precision sumev(2,3),sumqv(3,2),
     .  f(3,nbas,*)
      double complex smrho(k1,k2,k3,3,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer stdo,stdl,ib,ipl,is,kmax,lmxa,lmxh,lmxl,m,nelt1,nelt2,
     .  nelt3,nlma,nlmh,lgunit,ipr,iprint,nsp,nglob,nkaph
      double precision a,b,c,def1,disc,eadd,efold,p(3),px,py,pz,
     .  q1,q2,q3,a1,a2,a3,rsm,sen,sev,sqv,sav,w1,w2,w3,x,x1,x2,xlin
      equivalence (px,p(1)),(py,p(2)),(pz,p(3))

      stdo = lgunit(1)
      stdl = lgunit(2)
      ipr = iprint()
      ipl = ipr
      lrep = 0
      nsp  = nglob('nsp')
      nkaph = nglob('nkaph')

C ... Semiconductors: don't do anything
      if (lmet == 0) then
        sqv = sumqv(1,1) + sumqv(1,2)
C       smm = sumqv(1,2) - sumqv(1,2)
        sev = sumev(1,1)
        if (ipr >= 20) write(stdo,800) sqv,sev
  800   format(/' Semiconductor branch:  sumqv=',f12.6,
     .     '   sumev=',f12.6)
        return
      endif

C --- Metals: make weights w1,w2,w3 for the three trial Ef's ---
      q1 = sumqv(1,1) + sumqv(1,2)
      q2 = sumqv(2,1) + sumqv(2,2)
      q3 = sumqv(3,1) + sumqv(3,2)
      a1 = sumqv(1,1) - sumqv(1,2)
      a2 = sumqv(2,1) - sumqv(2,2)
      a3 = sumqv(3,1) - sumqv(3,2)
      if (ipr >= 30) write(stdo,300) qval,ef0-def,ef0,ef0+def,q1,q2,q3
  300 format(/' Interpolate to Fermi energy for qval=',f9.3
     .   /' trial Fermi energies: ',3f12.5
     .   /' band occupations:     ',3f12.5)
      if (ipr >= 30 .and. nsp == 2) write(stdo,301) a1,a2,a3
  301 format(' magnetic moment:      ',3f12.5)
      if (dabs(q3-q1) < 1d-8) then
        if (dabs(q2-qval) > 1d-8) then
          call awrit1('mshn3p: Input Fermi energy unusable, Q=%d',' ',
     .      80,stdo,q2)
          sumev(1,1) = 0
          sumev(2,1) = 0
          lrep = 1
          return
        endif
        w1 = 0d0
        w2 = 1d0
        w3 = 0d0
        x = 0
        if (ipr >= 10) write(stdo,430) w1,w2,w3
  430   format(' weights (converged):  ',3f12.5)
        goto 90
      endif
      if ((qval > q3 .or. qval < q1) .and. ipr > 1) call awrit0(
     .  ' mshn3p: qval outside window; extrapolating',' ',80,stdo)
      a = (q1+q3-2d0*q2)/(2d0*def*def)
      b = (q3-q1)/(2d0*def)
      c = q2-qval
      disc = b*b - 4d0*a*c
      if (qval <= q2) xlin = -def*(qval-q2)/(q1-q2)
      if (qval > q2) xlin = def*(qval-q2)/(q3-q2)
      if (disc < 0d0) then
        x = xlin
        if (qval >= q2) then
          w1 = 0d0
          w2 = 1d0 - x/def
          w3 = x/def
        else
          w1 = -x/def
          w2 = 1d0 + x/def
          w3 = 0d0
        endif
        if (ipr >= 10) write(stdo,450) w1,w2,w3
  450   format(' weights from lin fit: ',3f12.5)
      else
        x1 = (-b+dsqrt(disc))/(2d0*a)
        x2 = (-b-dsqrt(disc))/(2d0*a)
        x = x1
        if (dabs(x2-xlin) < dabs(x1-xlin)) x=x2
        if (x1*x2 < 0d0) then
          if (qval > q2) x = dmax1(x1,x2)
          if (qval <= q2) x = dmin1(x1,x2)
        endif
        w1 = 0.5d0*(x/def-1d0)*x/def
        w2 = 1d0-(x/def)**2
        w3 = 0.5d0*(x/def+1d0)*x/def
        if (ipr >= 10) write(stdo,460) w1,w2,w3
  460   format(' weights from quad fit:',3f12.5)
      endif

C ... Interpolate sumev,sumqv
  90  sqv = w1*sumqv(1,1) + w2*sumqv(2,1) + w3*sumqv(3,1)
     .    + w1*sumqv(1,2) + w2*sumqv(2,2) + w3*sumqv(3,2)
      sav = w1*sumqv(1,1) + w2*sumqv(2,1) + w3*sumqv(3,1)
     .    - w1*sumqv(1,2) - w2*sumqv(2,2) - w3*sumqv(3,2)
      sev = w1*sumev(1,1) + w2*sumev(1,2) + w3*sumev(1,3)
      sen = w1*sumev(2,1) + w2*sumev(2,2) + w3*sumev(2,3)
      if (ipr >= 20) write(stdo,470)
     .  q1,q2,q3,sqv,
     .  sumev(1,1),sumev(1,2),sumev(1,3),sev,
     .  sumev(2,1),sumev(2,2),sumev(2,3),sen
  470 format(' interpolated charge:  ',3f12.5,'  ->',f11.5
     .   /' eigenvalue sum:       ',3f12.5,'  ->',f11.5
     .   /' entropy term:         ',3f12.5,'  ->',f11.5)
      if (ipr >= 20 .and. nsp == 2) write(stdo,471) a1,a2,a3,sav
  471 format(' interpolated moment:  ',3f12.5,'  ->',f11.5)
      sumqv(1,1) = (sqv + sav)/nsp
      sumqv(1,2) = (sqv - sav)/nsp
      sumev(1,1) = sev
      sumev(2,1) = sen

C --- Interpolate coeffs of local density ---
      if (lrout > 0) then
      do  ib = 1, nbas
        is = s_site(ib)%spec
        p = s_site(ib)%pos
        rsm = s_spec(is)%rsma
        lmxl = s_spec(is)%lmxl
        lmxh = s_spec(is)%lmxb
        lmxa = s_spec(is)%lmxa
        kmax = s_spec(is)%kmxt
C       call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)

        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2

C   ... Case Pkl*Pkl, P*H, H*H
        if (lmxa > -1) then
        nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
        nelt2 = (kmax+1)*nkaph*nlma*nlmh
        nelt3 = nkaph*nkaph*nlmh*nlmh
        call mshn31(w1,w2,w3,nelt1,nsp,s_site(ib)%qkkl)
        call mshn31(w1,w2,w3,nelt2,nsp,s_site(ib)%qhkl)
        call mshn31(w1,w2,w3,nelt3,nsp,s_site(ib)%qhhl)
        endif

      enddo

C --- Interpolate smooth density ---
      call mshn32(w1,w2,w3,n1,n2,n3,k1,k2,k3,nsp,smrho)
C     call zprm3('sm rho-out',0,smrho,k1,k2,k3*nsp)

C --- Interpolate forces ---
      if (lfrce /= 0) then
        do  ib = 1, nbas
          do  m = 1, 3
            f(m,ib,1) = w1*f(m,ib,1) + w2*f(m,ib,2) + w3*f(m,ib,3)
          enddo
        enddo
      endif
      endif

C --- Choose new energy window ---
      efold = ef0
      def1 = def
      eadd = dmin1(x,1.0d0)
      ef0 = ef0+eadd
      def = 2d0*dabs(eadd)
      def = dmin1(def,1.0d0,def1*4d0)
      def = dmax1(def,0.001d0,def1/4d0)
      if (ipr >= 20) write(stdo,455) efold,def1
  455 format(' old energy window:   ef=',f10.5,'   plusminus',f9.5)
      if (ipr >= 10) write(stdo,456) ef0,def
  456 format(' new energy window:   ef=',f10.5,'   plusminus',f9.5)
      if (ipl > 1) write(stdl,711) w1,w2,w3,ef0,def,sev,sen
  711 format('nf w',3f7.3,'  ef',2f8.5,'   sev',f15.9,'  -sgS',f10.6)

      end

      subroutine mshn31(w1,w2,w3,nelt,nsp,qkkl)
C- Interpolate coeffs of local density
      implicit none
C ... Passed parameters
      integer nelt,nsp
      double precision w1,w2,w3,qkkl(nelt,3,nsp)
C ... Local parameters
      integer i,isp

      do  isp = 1, nsp
      do  i = 1, nelt
        qkkl(i,1,isp) = w1*qkkl(i,1,isp) +
     .                  w2*qkkl(i,2,isp) +
     .                  w3*qkkl(i,3,isp)
      enddo
      enddo
      if (nsp == 2) call dcopy(nelt,qkkl(1,1,2),1,qkkl(1,2,1),1)
      end

      subroutine mshn32(w1,w2,w3,n1,n2,n3,k1,k2,k3,nsp,smrho)
C- Interpolate smooth mesh density
      implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,nsp
      double precision w1,w2,w3
      double complex smrho(k1,k2,k3,3,nsp)
C ... Local parameters
      integer i1,i2,i3,isp

      do  isp = 1, nsp
      do  i3 = 1, n3
        do  i2 = 1, n2
          do  i1 = 1, n1
            smrho(i1,i2,i3,1,isp) = w1*smrho(i1,i2,i3,1,isp)
     .                            + w2*smrho(i1,i2,i3,2,isp)
     .                            + w3*smrho(i1,i2,i3,3,isp)
          enddo
        enddo
      enddo
      enddo
      if (nsp == 2)
     .  call dcopy(2*k1*k2*k3,smrho(1,1,1,1,2),1,smrho(1,1,1,2,1),1)
      end
