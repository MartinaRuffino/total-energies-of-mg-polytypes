      subroutine atwf(mode,a,lmxa,nr,nsp,pnu,pnz,rsml,ehl,rmt,z,v0,
     .  nphimx,ncore,konfig,ecore,gcore,gval)
C- Make properties related to core for one sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return ncore, and konfig, and nphimx only;
Ci         :  see description below for contents of nphimx
Ci         :1s digit
Ci         :1 return valence wave functions
Ci         :2 return core wave functions
Ci         :3 combination of 1+2
Ci         :10s digit concerns orthogonalization
Ci         :0 do not orthogonalize
Ci         :1 return orthogonalized to valence orbitals
Ci         :2 return orthogonalized to valence orbitals
Ci         :  using large component only
Ci         :100s digit concerns normalization
Ci         :1000s digit concerns treatment of spin pol core
Ci         :0 normal treatment
Ci         :1 Core w.f. from average of spin-up, spin-down potential
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci         :(not used if mode=0)
Ci   lmxa  :augmentation l-cutoff
Ci   nr    :number of radial mesh points.  Not used if mode=0.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :pnu por local orbitals
Ci   ... The following are not used if mode=0
Ci   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
Ci   ehl   :energy of smoothed Hankel tail for extended local orbital
Ci   rmt   :MT boundary, in a.u.
Ci   z     :nuclear charge
Ci   v0    :spherical potential
Cio Inputs/Outputs
Cio  nphimx:dimensions gval.  Must be at least as large as the
Cio        :number of valence wave functions
Cio        :For mode=0, nphimx is output and is assigned to
Ci         :maximum number radial wave functions for any l channel.
Co Outputs
Co   ncore :number of core levels
Co   konfig:1s digit contains core configuration
Co         :10s digit:
Co         : 0 -> no local orbitals
Co         : 1 -> local orbital with p.q.n. < pnu
Co         : 2 -> local orbital with p.q.n. > pnu
Co   ... The following are not returned if mode=0
Co   ecore :core eigenvalues
Co   gcore :core wave functions
Co   gval  :valence wave functions
Co          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
Co            i=0 : phi
Co            i=1 : phidot
Co            i=2 : local orbital
Cr Remarks
Cu Updates
Cu   14 Nov 14 Calls makrwf with mode=100 for semicore states
Cu             Note: this will affect program execution slightly
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu    4 Sep 04 Adapted to extended local orbitals
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   22 Apr 01 Created by MvS
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nr,nsp,lmxa,ncore,konfig(1+lmxa),n0,nrmx,nphimx
      parameter (n0=10,nrmx=5001)
      double precision rmt,z,a,v0(nr,nsp),pnu(n0,nsp),pnz(n0,nsp),
     .  gval(nr*2,0:lmxa,nphimx,nsp),ecore(*),gcore(nr,2,*),rsml(n0),ehl(n0)
C ... Local parameters
      logical lpz,ltmp,isanrg
      integer l,isp,konf,konfz,k,mode0,mode1,mode3,intopt,j
      double precision sumtc,sumec,e,ez,xx
C     double precision hcrl,val(5),slo(5),pi,tol
C     parameter (tol=1d-12)
      double precision rofi(nrmx),rwgt(nrmx),rhoc(nrmx,2),gp(2*nrmx*4)
      double precision phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp
      procedure(integer) :: nglob

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
C     mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)

C --- Count number of core states ---
      lpz = .false.
      ncore = 0
      do  l = 0, lmxa
        k = l+1
        konfig(k) = pnu(k,1)
        konfz = mod(pnz(k,1),10d0)
        if (konfz == 0) konfz = konfig(k)
C       Sanity check
        ltmp=isanrg(konfz,konfig(k)-1,konfig(k)+5,'atwf:','pnuz',.true.)
C       lpz = konfz /= konfig(k)
        do  konf = l+1, min(konfz,konfig(k))-1
          ncore = ncore+nsp
        enddo
        if (konfz < konfig(k)) then
          konfig(k) = konfz + 10
          lpz = .true.
        elseif (konfz > konfig(k)) then
          konfig(k) = konfig(k) + 20
          lpz = .true.
        endif
      enddo

      if (mode0 == 0) then
        nphimx = 2
        if (lpz) nphimx = 3
        return
      endif

      if (nr > nrmx) call rx('increase nrmx in atwf')
      call radmsh(rmt,a,nr,rofi)
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rmt,a,nr,rwgt)

C --- Valence wave functions ---
      if (mod(mode0,2) == 1) then
        do  l = 0, lmxa
        k = l+1
        do  isp = 1, nsp
          konf = pnu(k,1)

C    ...  Make phi and phidot
C         NB: Write gdot to gp, with extra space for higher derivatives
C         nn  = konf-l-1
C         pi = 4d0*datan(1d0)
C         hcrl = 0
C         val(1) = rofi(nr)
C         slo(1) = 1 + dtan(pi*(0.5d0 - pnu(k,isp)))
C         call phidx(0,z,l,v0(1,isp),hcrl,0d0,rofi,nr,4,tol,e,val,slo,
C    .      nn,gval(1,l,1,isp),gp,xx,xx,xx,xx,pgam,xx,xx,xx,xx)

          call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnu(k,isp),0d0,4,
     .      gval(1,l,1,isp),gp,e,phi,dphi,phip,dphip,p)
C         Copy 1st derivative to passed array
          call dcopy(2*nr,gp,1,gval(1,l,2,isp),1)
C         phi,phidot already orthogonal if mode1=1
          if (mode1 == 2)
     .     call ortrwf(10*(mode1-1)+2,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt,
     .      e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)

C     ... Make local orbital
          if (konf /= konfig(k)) then
            ltmp = isanrg(nphimx,3,3,'atwf:','nphimx',.true.)
            j = 0
            if (int(mod(pnz(l+1,isp),10d0)) < int(pnu(l+1,1))) j = 100
            call makrwf(j,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnz(l+1,isp),0d0,2,
     .        gval(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)

            ltmp = isanrg(mode1,0,2,'atwf:','10s digit mode',.true.)
            if (mode1 == 0) then

C             Extra scaling
C              call ortrwf(0,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt,e,e,ez,
C     .          gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)
C             call dscal(nr*2,1/xx,gval(1,l,3,isp),1)
C             phz = phz/xx
C             dphz = dphz/xx

              call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz,
     .          phzp,dphzp,pnz(1,isp),rsml,ehl,
     .          gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp))
            elseif (mod(pnz(l+1,isp),100d0) < 10) then
              call ortrwf(10*(mode1-1)+1,z,l,v0(1,isp),nr,nr,nr,rofi,
     .          rwgt,e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp),
     .          gval(1,l,3,isp),xx)
            endif
C           call prrmsh('gz',rofi,gval(1,l,3,isp),nr,nr,2)
          endif

        enddo
        enddo

C        call prrmsh('gval (phi)',rofi,gval,nr,nr,2*(1+lmxa))
C        call prrmsh('gval (dot)',rofi,gval(1,0,2,1),nr,nr,2*(1+lmxa))

      endif

C --- Core eigenfunctions and eigenvalues ---
      if (mode0 >= 2) then
        call getcor(10*mode3+1,z,a,pnu,pnz,nr,lmxa,rofi,v0,0,0,0d0,
     .    sumec,sumtc,rhoc,ncore,ecore,gcore)
C       Check normalization:
C       mc -f2f15.10 out.gas -av:nr,1 rmax -p -coll 1 -a r -coll 2:nc -p -xe \
C       r -tog -ccat -int:nord=6 0 rmax \
C       -coll 2:nc -p -coll 1:nc:2 -tog -coll 2:nc:2 -rcat -t -csum
C       call prrmsh('gcore',rofi,gcore,nr,nr,2*ncore)
      endif

      end

      subroutine atwf2l(ifi,jfi,iclass,a,lmxa,nr,nsp,pnz,rmt,
     .  nphimx,konfig,ppnl,ecore,gcore,gval)
C- Translate radial wave functions from shifted to standard log mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit for valence w.f.
Ci   jfi   :file logical unit for core w.f.
Ci   iclass:for printout
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci         :See Remarks
Ci   lmxa  :augmentation l-cutoff
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnz   :p.q.n for local orbitals
Ci   rmt   :MT boundary
Ci   nphimx:dimensions gval.  Should be 2, or 3 if any local orbitals
Ci   konfig:1s digit contains core configuration
Ci         :10s digit:
Ci         : 0 -> no local orbitals
Ci         : 1 -> local orbital with p.q.n. < pnu
Ci         : 2 -> local orbital with p.q.n. > pnu
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci         :ppnl(2) = s00 = <phi | phi>
Ci         :ppnl(7) = s11 = <phidot | phidot>
Ci         :ppnl(8)  = szz = <gz|gz> for local orbital
Ci   ecore :core eigenvalues
Ci   gcore :core wave functions
Ci   gval  :valence wave functions
Ci          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
Ci            i=0 : phi
Ci            i=1 : phidot
Ci            i=2 : local orbital
Cr Remarks
Cr   This routine translates w.f. to a standard log mesh and writes
Cr   them in a spex-readable format.
Cr
Cr   Debugging:
Cr   set sqr = '-p -p -xe -tog -coll 1 -tog -coll 2:nc -ccat'
Cr   mc -qr out.gas $sqr -int 0 2.361911
Cu Updates
Cu   15 Jul 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,jfi,iclass,nr,nsp,lmxa,konfig(0:lmxa),nphimx
      integer n0,nppn,nphmxx
      parameter (n0=10,nppn=12,nphmxx=3)
      double precision rmt,a,pnz(n0),gval(nr,2,0:lmxa,nphimx,nsp),
     .  ecore(*),gcore(nr,2,*),ppnl(nppn,n0,nsp)
C ... Local parameters
      character phitype(nphmxx)*3
      integer ic,icore,iphi,iprint,ir,isp,jcore,jx,konf,l,nphi,
     .  npoly,nrmx,stdo,PRTV,intopt,nrl
      double precision xx1,xx2,xx1n,xx2n,fac,tolg,dot3,sphi(3),ren(2,2)
      double precision x,y,dy,dytop,dytopa  ! For interpolation
      parameter (npoly=6,nrmx=5001,PRTV=60)
      parameter (tolg=1d-8)
      double precision
     .  gwk(nrmx),gwk2(nrmx,2,nphimx),  ! work arrays
     .  rofi(nrmx),rwgt(nrmx),   ! Points and weights for shifted log mesh
     .  rwgtl(nrmx),rofil(nrmx)  ! Points and weights for stnd log mesh
      double precision onorm(2,0:lmxa,nphimx,nsp), ! norm on shifted mesh
     .                  norm(2,0:lmxa,nphimx,nsp)  ! norm on log mesh
C     double precision sum1,sum2
      procedure(integer) :: nglob
      data phitype /'phi','dot','phz'/

C ... Setup
      stdo = nglob('stdo')
      call togpr
      call info0(PRTV,1,0,' ... Transpose radial functions to log mesh')
      if (nr > nrmx) call rx('increase nrmx in atwf2l')

C ... Shifted and standard log mesh
      call radmsh(rmt,a,nr,rofi)
      call mshnrl(1,rmt,a,nr,nrmx,nrl)
C     nrl = nr
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rmt,a,nr,rwgt)  ! weights for shifted mesh
      call radmwt(1000+intopt,rmt,a,nrl,rofil,rwgtl)  ! mesh, weights for std log mesh

C      call prrmsh('standard log mesh',rofi,rofil,nr,nr,1)
C      call prrmsh('standard log weights',rofil,rwgtl,nrl,nrl,1)

C      call info5(0,0,0,'testing nr=%i  a=%d  rmt=%d ... int x dx',
C     .  nr,a,rmt,0,0)
C      sum1 = 0
C      sum2 = 0
C      do  ir  = 1, nr
C        sum1 = sum1 + rofi(ir)*rwgt(ir)
C        sum2 = sum2 + rofil(ir)*rwgtl(ir)
C      enddo
C      print 333, rmt**2/2,sum1-rmt**2/2,sum2-rmt**2/2
C  333 format('exact',f15.10,'  shifted mesh error',1pe10.2,
C     .  '  standard mesh error',1pe10.2)

C --- Valence wave functions ---
C     ic = 2*(lmxa+1)*nphimx*nsp
C     call prrmsh('gval on shifted log mesh',rofi,gval,nr,nr,ic)

      ren(1,1) = 999
      ren(2,1) = 0
      ren(1,2) = 999
      ren(2,2) = 0

      write(ifi) iclass,lmxa,nsp,nrl ! 1st record in class: dimensioning
      call info2(PRTV,1,0,' Valence wave function normalizations for '
     .  //'class %i%N%7fl   spin  mesh'//
     .  '%5flarge%10fsmall%10fsum%12fscale',iclass,isp)
      dytopa = 0
      do  isp = 1, nsp
      do  l = 0, lmxa
        nphi = 2
        sphi(1) = ppnl(2,l+1,isp)
        sphi(2) = ppnl(7,l+1,isp)
        if (pnz(l+1) > 0) then
          nphi = 3
          sphi(3) = ppnl(8,l+1,isp)
        endif
C   ... For each radial function, transpose to standard log mesh
        do  iphi = 1, nphi

C         Alternative way to make xx1=onorm(1,l,iphi,isp)+onorm(2,l,iphi,isp)
C         call addrwf(2,xx1,l,xx1,nr,nr,nr,rofi,rwgt,xx1,xx1,0d0,
C    .      gval(1,1,l,iphi,isp),gval(1,1,l,iphi,isp),xx1)

          do  ic = 1, 2         ! large, small components
            onorm(ic,l,iphi,isp) =
     .        dot3(nr,gval(1,ic,l,iphi,isp),gval(1,ic,l,iphi,isp),rwgt)

C           g(r) = phi(r)*r
C           For small r, g(r)~r^(l+1) ... so fit g(r)/r^(l+1)
            do  ir = 2, nr
              gwk(ir) = gval(ir,ic,l,iphi,isp)/rofi(ir)**(l+1)
            enddo
            jx = 0
            dytop = 0
            do  ir = 1, nrl
              x = rofil(ir)
              call polint(rofi(2),gwk(2),nr-1,npoly,x,tolg,0,jx,y,dy)
              dytop = max(dytop,abs(dy))
              dytopa = max(dytopa,abs(dy))
              gwk2(ir,ic,iphi) = y*rofil(ir)**(l+1)
            enddo
            norm(ic,l,iphi,isp) =
     .        dot3(nrl,gwk2(1,ic,iphi),gwk2(1,ic,iphi),rwgtl)
C            Debugging printout
C            call info8(10,0,0,'interp ic=%i l=%i iphi=%i isp=%i:'//
C       .      ' err ~ %;g  onorm = %;6d  nnorm = %;6d',
C       .      ic,l,iphi,isp,dytop,onorm(ic,l,iphi,isp),
C       .      norm(ic,l,iphi,isp),0)

          enddo
        enddo

C   ... Orthogonalize phidot to phi; get new norm for phi
        iphi = 2
        call ortrwf(23,xx1,l,xx1,nrmx,nrl,nrl,rofil,rwgtl,
     .    xx1,xx1,xx1,gwk2(1,1,1),gwk2(1,1,iphi),gwk2,xx1)
        do  ic = 1, 2           ! large, small components
          norm(ic,l,iphi,isp) =
     .      dot3(nrl,gwk2(1,ic,iphi),gwk2(1,ic,iphi),rwgtl)
        enddo

C   ... Scale radial wave functions to match ppnl; write to disk
        do  iphi = 1, nphi

          if (iprint() >= PRTV) then
            xx1 = onorm(1,l,iphi,isp)
            xx2 = onorm(2,l,iphi,isp)
            write(stdo,344) phitype(iphi), l, isp,
     .        's',xx1,xx2,xx1+xx2
            xx1 = norm(1,l,iphi,isp)
            xx2 = norm(2,l,iphi,isp)
            write(stdo,344) phitype(iphi), l, isp,
     .        'l',xx1,xx2,xx1+xx2,dsqrt(sphi(iphi)/(xx1+xx2))
  344       format(2x,a3,i3,i6,4x,a,4f15.8)
          endif

          xx1 = norm(1,l,iphi,isp)
          xx2 = norm(2,l,iphi,isp)
          fac = dsqrt(sphi(iphi)/(xx1+xx2))
          ren(1,1) = min(ren(1,1),fac)
          ren(2,1) = max(ren(2,1),fac)
          do  ic = 1, 2
            call dscal(nrl,fac,gwk2(1,ic,iphi),1)
          enddo

C        Sanity check
          call addrwf(2,xx1,l,xx1,nrmx,nrl,nrl,rofil,rwgtl,xx1,xx1,0d0,
     .      gwk2(1,1,iphi),gwk2(1,1,iphi),xx1)
C         print *, xx1,xx1-sphi(iphi)

          write(ifi) l,iphi,isp ! indices to wave function being written
          write(ifi) gwk2(1:nrl,1,iphi),gwk2(1:nrl,2,iphi) ! Large component, followed by small component

        enddo
C        call prrmsh('gval (phi)',rofil,gwk2(1,1,1),nrmx,nrl,2)
C        call prrmsh('gval (dot)',rofil,gwk2(1,1,2),nrmx,nrl,2)
      enddo
      enddo

      write(ifi) -1,-1,-1 ! Flags last record for this class
C     Debugging; need to comment overwrite above
C      ic = 2*(lmxa+1)*nphimx*nsp
C      call prrmsh('gval on standard log mesh',rofil,gval,nrl,nrl,ic)

C --- Core eigenfunctions and eigenvalues ---
      write(jfi) iclass,lmxa,nsp,nrl ! 1st record in class: dimensioning
      dytopa = 0
      call info2(PRTV,1,0,' Core wave function normalizations for '
     .  //'class %i%N   n   l   spin  mesh'//
     .  '%5flarge%10fsmall%10fsum%9fecore/scale',iclass,isp)
      icore = 0
      do  l = 0, lmxa
        do  isp = 1, nsp
          jcore = 0
          do  konf = l+1, mod(konfig(l),10)-1
            icore = icore+1
            jcore = jcore+1

            xx1 = dot3(nr,gcore(1,1,icore),gcore(1,1,icore),rwgt)
            xx2 = dot3(nr,gcore(1,2,icore),gcore(1,2,icore),rwgt)

            dytop = 0
            do  ic = 1, 2 ! large, small components

C             g(r) = phi(r)*r
C             For small r, g(r)~r^(l+1) ... so fit g(r)/r^(l+1)
              do  ir = 2, nr
                gwk(ir) = gcore(ir,ic,icore)/rofi(ir)**(l+1)
              enddo
              jx = 0
              do  ir = 1, nrl
                x = rofil(ir)
                call polint(rofi(2),gwk(2),nr-1,npoly,x,tolg,0,jx,y,dy)
                if (ic == 1) then
                  dy = dy*rofil(ir)**(l+1)
                  dytop  = max(dytop,dy)
                  dytopa = max(dytopa,dy)
                endif
                gwk2(ir,ic,1) = y*rofil(ir)**(l+1)
C               gcore(ir,ic,icore) = gwk2(ir,ic,1) ! Overwrite, for debuggging
              enddo
            enddo

            xx1n = dot3(nrl,gwk2(1,1,1),gwk2(1,1,1),rwgtl)
            xx2n = dot3(nrl,gwk2(1,2,1),gwk2(1,2,1),rwgtl)
C           Debugging printout
            call info8(99,0,0,' atwf2l interp l=%i konf=%i isp=%i:'//
     .        ' err ~ %;g  onorm = %,6;6d  nnorm-onorm = %;3g',l,
     .        konf,isp,dytop,xx1,xx1n-xx1,0,0)

C       ... Normalize to unity inside MT sphere
            fac = dsqrt(1d0/(xx1n+xx2n))
            ren(1,2) = min(ren(1,2),fac)
            ren(2,2) = max(ren(2,2),fac)
            do  ic = 1, 2
              call dscal(nrl,fac,gwk2(1,ic,1),1)
            enddo

            if (iprint() >= PRTV) then
              write(stdo,345) konf, l, isp, 's',
     .        xx1, xx2, xx1+xx2, ecore(icore)
              write(stdo,346) konf, l, isp, 'l',
     .          xx1n, xx2n, xx1n+xx2n, fac
  345         format(2i4,i6,4x,a,3f15.8:f15.6)
  346         format(2i4,i6,4x,a,4f15.8)
            endif

C           Sanity check
C           call addrwf(2,xx1,l,xx1,nrmx,nrl,nrl,rofil,rwgtl,xx1,xx1,
C    .        0d0,gwk2(1,1,1),gwk2(1,1,1),xx1)
C           print *, xx1-1

            write(jfi) jcore, l, isp, konf, ecore(icore)
            write(jfi) gwk2(1:nrl,1,1), gwk2(1:nrl,2,1)

          enddo
        enddo
      enddo
      write(jfi) -1,-1,-1,-1,-1d0 ! Flags last record for this class
      call info5(30,0,0,
     .  ' renorm:  min val %,8;8d  max val %,8;8d'//
     .  '  min core %,8;8d  max core %,8;8d',
     .  ren(1,1),ren(2,1),ren(1,2),ren(2,2),0)

      call togpr

      end
