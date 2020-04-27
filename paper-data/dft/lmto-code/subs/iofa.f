      integer function iofa(mode,spid,nxi0,nxi,exi,hfc,hfct,rsm,z,rmt,a,nr,nsp,
     .  qc,ccof,ceh,stc,rho,rhoc,v,ifi)
C- I/O for free-atom data, one species
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :indicates what to read from file (read) or what to write to file (write)
Ci         :bit
Ci         :1  : dimensioning parameters z,a,nsp,lrel,nr,rmt,rsm
Ci         :2  : (read only) do not read dimensioning parameters
Ci         :   : but confirm that they match passed values
Ci         :4  : Interstitial fit nxi; exi,hfc,hfct
Ci         :8  : valence potential
Ci         :16 : valence density
Ci         :32 : core density and core parameters: qc,ccof,ceh,stc
Ci   spid  :species label
Ci   nxi0  :leading dimension of hfc,hfct
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio File I/O
Cio  nxi   :number of energies used to fit tail of valence density
Cio  exi   :energies that fit tail of valence density
Cio  hfc   :coefficients for fit of tail valence density
Cio  hfct  :coefficients for fit of tail valence+core density
Cio  rsm   :smoothing radius for fit of tail valence density
Cio  z     :nuclear charge
Cio  rmt   :muffin-tin radius, in a.u.
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio  nr    :number of radial mesh points
Cio  nsp   :number of spin channels
Cio  qc    :Sphere core charge
Cio  ccof  :coefficient to core density fit by unsmoothed Hankel
Cio  ceh   :hankel function energy to core density fit
Cio  stc   :core kinetic energy
Cio  rho   :valence density
Cio  rhoc  :core density
Cio  v     :spherical potential
Cr Remarks
Cu Updates
Cu   07 Mar 18 iofa and read/write sections of atomic data.  Format changed
Cu   10 Jun 00 spin polarized
Cu   20 May 00 adapted from nfp rw_fa.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ifi,nr,nxi,nxi0,nsp
      character spid*8
      double precision a,ccof,ceh,qc,rmt,rsm,stc,z,hfc(nxi0,2),
     .  hfct(nxi0,2),exi(nxi0)
C     nr isn't known on entry
C     double precision rho(2*nr),v(2*nr),rhoc(2*nr)
C     double precision rho(*),v(*),rhoc(*)
      real(8), target :: rho(*),v(*),rhoc(*)
C ... Local parameters
      real(8), pointer :: prho(:)
      real(8), parameter:: tol=1d-6
      integer i,jfi,lrel,lrel0,nsp0,nr0,mode0,ipr,stdo
      double precision z0,a0,rmt0,rsm0
      procedure(integer) :: nglob,lgunit,iprint
      character msg*23, strn*10

      ipr    = iprint()
      stdo = nglob('stdo')
      lrel = mod(nglob('lrel'),10)
      msg  = '         File mismatch:'
      iofa   = -1

C --- Input ---
      if (ifi > 0) then

        jfi = ifi
        read(jfi,201,end=998,err=998) spid
  201   format(9x,a8)
        read(jfi,*,end=998,err=998) strn, mode0
        if (strn /= 'contents') goto 998

        if (IAND(mode,1) > 0 .and. IAND(mode0,1) == 0) goto 998
        if (IAND(mode,4) > 0 .and. IAND(mode0,4) == 0) goto 998
        if (IAND(mode,8) > 0 .and. IAND(mode0,8) == 0) goto 998
        if (IAND(mode,16) > 0 .and. IAND(mode0,16) == 0) goto 998
        if (IAND(mode,32) > 0 .and. IAND(mode0,32) == 0) goto 998

C   ... Read dimensioning parameters and sanity checks
        if (IAND(mode0,1) > 0) then
          read(jfi,*) strn
          read(jfi,102) z0,a0,nsp0,lrel0,nr0,rmt0,rsm0
          if (IAND(mode,2) == 0 .and. IAND(mode,1) == 1) then
            z=z0; a=a0; nsp=nsp0; lrel=lrel0; nr=nr0; rmt=rmt0; rsm=rsm0
          endif
          if (IAND(mode,1) == 1) then
            call fsanrg(z,z0,z0,tol,msg,'z',.true.)
            call fsanrg(a,a0,a0,tol,msg,'a',.true.)
            call fsanrg(rmt,rmt0,rmt0,tol,msg,'rmt',.true.)
            call fsanrg(rsm,rsm0,rsm0,tol,msg,'rsm',.true.)
            call sanrg(.true.,nr0,nr,nr, msg,'nr')
            call sanrg(.true.,nsp0,nsp,nsp, msg,'nsp')
            call sanrg(.true.,lrel0,lrel,lrel,msg,'lrel')
          endif
        endif

C   ... Read parameters related to interstitial fit
        if (IAND(mode0,4) > 0) then
          read(jfi,*,end=998,err=998) strn, i
          if (strn /= 'nxi') goto 998
          if (IAND(mode,4) == 0) then
            read(jfi,*)
            read(jfi,*)
            read(jfi,*)
            if (nsp == 2) read(jfi,*)
            if (nsp == 2) read(jfi,*)
          else
            nxi = i
            read(jfi,*) (exi(i),i=1,nxi)
            read(jfi,*) (hfc(i,1),i=1,nxi)
            read(jfi,*) (hfct(i,1),i=1,nxi)
            if (nsp == 2) read(jfi,*) (hfc(i,2),i=1,nxi)
            if (nsp == 2) read(jfi,*) (hfct(i,2),i=1,nxi)
          endif
        endif

C   ... Read parameters related to core density
        if (IAND(mode0,32) > 0 .and. IAND(mode,32) > 0) then
          read (jfi,210) qc,ccof,ceh,stc
  210     format(5x,4f16.7)
        elseif (IAND(mode0,32) > 0) then
          read (jfi,210)
        endif

        allocate(prho(nr))
        if (IAND(mode0,16) > 0) then
          read(ifi,"(a)") strn
          if (strn /= 'rho') goto 998
          if (IAND(mode,16) > 0) then
            call dfdump(rho,nr,ifi)
            if (nsp == 2) call dfdump(rho(1+nr),nr,ifi)
          else
            call dfdump(prho,nr,ifi)
            if (nsp == 2) call dfdump(prho,nr,ifi)
          endif
        endif

        if (IAND(mode0,32) > 0) then
          read(ifi,"(a)") strn
          if (strn /= 'rhoc') goto 998
          if (IAND(mode,32) > 0) then
            call dfdump(rhoc,nr,ifi)
            if (nsp == 2) call dfdump(rhoc(1+nr),nr,ifi)
          else
            call dfdump(prho,nr,ifi)
            if (nsp == 2) call dfdump(prho,nr,ifi)
          endif
        endif

        if (IAND(mode0,8) > 0) then
          read(ifi,"(a)") strn
          if (strn /= 'v0') goto 998
          if (IAND(mode,8) > 0) then
            call dfdump(v,nr,ifi)
            if (nsp == 2) call dfdump(v(1+nr),nr,ifi)
          else
            call dfdump(prho,nr,ifi)
            if (nsp == 2) call dfdump(prho,nr,ifi)
          endif
        endif

        deallocate(prho)
      endif

C --- Output ---
      if (ifi < 0)  then
        jfi = -ifi
        write(jfi,101) spid,mode
  101   format('-------- ',a8,' ----------'/'contents ',i4)
        if (IAND(mode,1) > 0) then
          write(jfi,"('     z          a      nsp lrel  nr       rmt          rsm')")
          write(jfi,102) z,a,nsp,lrel,nr,rmt,rsm
  102     format(f8.3,f12.7,2i5,i6,2f13.7,i5)
        endif
        if (IAND(mode,4) > 0) then
          write(jfi,103) nxi
  103     format('nxi',i4)
          write(jfi,105) (exi(i),i=1,nxi)
          write(jfi,105) (hfc(i,1),i=1,nxi)
          write(jfi,105) (hfct(i,1),i=1,nxi)
          if (nsp == 2) write(jfi,105) (hfc(i,2),i=1,nxi)
          if (nsp == 2) write(jfi,105) (hfct(i,2),i=1,nxi)
  105     format(1p,8d16.8)
        endif
        if (IAND(mode,32) > 0) then
          write (jfi,110) qc,ccof,ceh,stc
  110     format('core',4f16.7)
        endif
        if (IAND(mode,16) > 0) then
          write(jfi,"('rho')")
          call dfdump(rho,nr,-jfi)
          if (nsp == 2) call dfdump(rho(1+nr),nr,-jfi)
        endif
        if (IAND(mode,32) > 0) then
          write(jfi,"('rhoc')")
          call dfdump(rhoc,nr,-jfi)
          if (nsp == 2) call dfdump(rhoc(1+nr),nr,-jfi)
        endif
        if (IAND(mode,8) > 0) then
          write(jfi,"('v0')")
          call dfdump(v,nr,-jfi)
          if (nsp == 2) call dfdump(v(1+nr),nr,-jfi)
        endif
      endif

      iofa = 0
      return

C ... Error handling
  998 if (ipr > 1) write(stdo,'('' iofa  : failed to read species .. nothing read'')')

      end
