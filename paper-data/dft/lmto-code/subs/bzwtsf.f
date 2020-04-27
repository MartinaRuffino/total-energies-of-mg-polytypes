      subroutine bzwtsf(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .  fmom,metal,ltet,norder,npts,width,rnge,wtkp,eb,efmax,lswtk,
     .  swtk,efermi,sumev,dosef,wtkb,qval,egap,lwtkb)
C- BZ integration for fermi level, band sum and qp weights, fixed-spin
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbmx  :leading dimension of eb
Ci   nevx  :leading dimension of wtkb and max number of evals calculated
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   n1..n3:number of divisions for the k-point mesh
Ci   nkp   :number of inequivalent k-points (bzmesh.f)
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of tetrahedron;
Ci         :idtet(0,i) number of tetrahedra of the i'th kind
Ci   zval  :valence charge
Ci   metal :0 => nonmetal
Ci         :1 => metal
Ci         :11 => Compute weights using input efermi as Fermi level
Ci   ltet  :0 => sampling integration
Ci         :1 => tetrahedron integration
Ci         :11 => tetrahedron integration for bands, sampling for weights
Ci   norder,npts,width,rnge: parameters for sampling integr. (maknos)
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   eb    :energy bands; alias eband
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   eb    :energy bands
Ci   efmax :eigenvalue limit for eigenvectors calculated in diagno
Ci   lswtk :Flags indicating whether 'spin weights' swtk are available
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Supplies information about spin moment in noncoll case.
Ci         :Used when lswtk is set
Cio Inputs/Outputs
Cio  fmom  :fixed spin moment.  If zero, no constraint is applied.
Cio        :(input) fmom(1) = constrained moment.  Internal global field
Cio                 is applied until moment fmom(1)
Cio                 fmom(1) = (0,NULLR) => no fixed-spin-moment applied
Cio        :(input)  fmom(2) = field applied before constraint
Cio        :(output) fmom(2) = field applied to satisfy constraint.
Cio        :(output) fmom(3) = field * moment (for total energy term)
Cio  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
Cio        :0 weights are not available; no moment calculation
Cio        :if 1, weights were generated with no constraint
Cio        :In this case, print moment, and if fmom ne 0 remake weights
Cio        :with constraint; set to lwtkb=2 on output.
Cio        :if 2, weights were generated with constrained global moment
Cio        :if -1, same as 1
Co Outputs
Co   efermi:Fermi energy (not altered if metal is 11)
Co   sumev :sum of eigenvalues
Co   dosef :DOS at the Fermi level
Co   wtkb  :integration weights (not generated for nonmetal case)
Co   qval  :qval(1) = total charge; qval(2) = magnetic moment
Co   egap  :energy gap (returned if metal=T and insulator found)
Co         :If none sought or metal found, egap returned as NULLI
Cu Updates
Cu   16 Feb 14 dosef now returned to calling program; new argument list
Cu   11 Oct 13 metal now an integer, with allowed values 0,1,11
Cu   09 Aug 13 ltet now an integer, with allowed values 0,1,11
Cu   31 Aug 12 Returns egap
Cu   19 Jul 10 Return effective field, double counting from FSM
Cu             New input/output fmom(1:2)
Cu   12 Jul 08 change arg list in bzwts -- now returns entropy term
Cu   02 Jan 06 return qval (valence charge and moment)
Cu   22 Sep 01 Adapted from bzwts.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,
     .  idtet(5,ntet),lswtk,lwtkb,ltet,metal
      double precision zval,fmom(3),eb(nbmx,nsp,nkp),width,rnge,
     .  wtkp(nkp),wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efmax,efermi,
     .  sumev,qval(2),egap,dosef
C ... Local parameters
      integer stdo,stdl,ipr,itmax,iter,iprint
      double precision amom,dosefl(2),vhold(12),vnow,dv,ef0,ent
      double precision dvcap,NULLR
      parameter (dvcap=.05d0,itmax=50,NULLR=-99999d0)
      procedure(integer) :: nglob
C ... External calls
      external awrit3,awrit5,bzbfde,bzwtsm,dpzero,dvdos,getpr,info0,
     .         info5,poppr,pshpr

      dosefl = 0

C --- Fermi level without spin constraint ---
C ... undo effects of fmom(2)
      if (nsp == 2 .and. fmom(2) /= 0 .and.
     .    fmom(1) /= NULLR .and. fmom(1) /= 0) then
          call bzbfde(nspc == 2,nkp,nsp,nbmx,nevx,-fmom(2),swtk,eb)
      endif
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .  metal,ltet,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .  sumev,wtkb,dosefl,qval,ent,egap)
      fmom(3) = 0
      if (nsp == 1) goto 99

      call getpr(ipr)
      stdo = nglob('stdo')
      stdl = nglob('stdl')

C --- Make and print out magnetic moment ---
      if ((lswtk == 1 .or. nspc == 1) .and. metal /= 0) then
       call bzwtsm(lswtk == 1.and.nspc == 2,nkp,nsp,nevx,wtkb,swtk,amom)
C       if (ipr >= 20) write(stdo,1) amom
C    1  format(9x,'Mag. moment:',f15.6)
       call info2(20,0,0,'%9fMag. moment:%;15,6D',amom,0)
       qval(2) = amom
      else
        call info0(20,0,0,
     .    '%9fspin weights not available ... no spin moment calculated')
        goto 99
      endif

C --- Setup for fixed-spin moment method ---
      if (fmom(1) == 0d0 .or. fmom(1) == NULLR .or. lwtkb == 0) then
C   ... Double counting from fmom(2)
        if (nsp == 2 .and. fmom(2) /= 0) then
          fmom(3) = amom*fmom(2)
          if (ipr >= 20) then
          call awrit3('%9fconst. moment:%;13,6D  Beff=%g  Beff.M=%,6d',
     .      ' ',160,stdo,amom,fmom(2),fmom(3))
          call awrit3('bz fmom %,6d  Beff %g  Beff.M %,6d',
     .      ' ',160,stdl,amom,fmom(2),fmom(3))
          endif
        endif
        goto 99
      endif

      call dpzero(vhold,12)
C     Already have no-field result as first guess
      vnow = 0
      call dvdos(vnow,amom,dosefl,vhold,fmom,dvcap,dv)
C     Start fresh if starting bands have no field
      if (fmom(2) == 0) call dpzero(vhold,12)
      vnow = fmom(2)*2
      ef0 = efermi
      call info0(41,1,0,' Seek potential shift for fixed-spin mom ...')
      iter = 0

C --- Entry point for next iteration (no effect if initial vnow=0) ---
C ... Potential shift
   10 continue
      call bzbfde(nspc == 2,nkp,nsp,nbmx,nevx,vnow/2,swtk,eb)

C ... Fermi level with dv shift
      if (iter >= 0) call pshpr(ipr-50)
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,metal,
     .           ltet,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .           sumev,wtkb,dosefl,qval,ent,egap)
      if (iprint() >= 20 .or. iter == 0 .and. iprint() >= -30) then
        call bzwtsm(lswtk == 1.and.nspc == 2,nkp,nsp,nevx,wtkb,swtk,
     .    amom)
        call awrit3('%9fconst. moment:%;13,6D  Beff=%g  Beff.M=%,6d',
     .              ' ',160,stdo,amom,vnow/2,amom*vnow/2)
        call awrit3('bz fmom %,6d  Beff %g  Beff.M %,6d',' ',160,stdl,
     .              amom,vnow/2,amom*vnow/2)
      endif

      if (iter >= 0) then
        call poppr
        call bzbfde(nspc == 2,nkp,nsp,nbmx,nevx,-vnow/2,swtk,eb)

C   ... Re-entry point for new guess at potential shift
C       iterate to itmax times.  iter<0 => last iteration
        iter = iter + 1

C ...   Magnetic moment
        call bzwtsm(lswtk == 1.and.nspc == 2,nkp,nsp,nevx,wtkb,swtk,
     .    amom)
        if (ipr >= 41) call awrit5(' dv=%;10,6D  yields '//
     .    'ef=%;10,6D  amom=%;10,6D;  seeking %;10,6D',' ',160,
     .    stdo,vnow,efermi,amom,fmom,vnow)
        call dvdos(vnow,amom,dosefl,vhold,fmom,dvcap,dv)
        if (abs(dv) < 1d-6) then
C         A root was found
          if (vhold(12) == -2 .or. vhold(12) == -3 .or.
     .        vhold(12) == 0 .or. vhold(12) == 1) then
            call info5(10,1,0,' BZWTSF: potential shift bracketed.'//
     .        '  Unconstrained efermi=%,6;6d'//
     .        '%N cnst fmom=%,6;6d  actual mom=%,6;6d'//
     .        '  ef=%,6;6d  dv=%,6;6d',
     .        ef0,fmom,amom,efermi,vnow)
            iter = -iter
          endif
        else if (iter == itmax) then
          if (ipr >= 10)
     .      call awrit5('%N BZWTSF: failed to converge potential shift'
     .      //' after %i iterations.'//
     .      '%N constraint fmom=%,6;6d  actual amom=%,6;6d'//
     .      '  ef=%,6;6d  dv=%,6;6d',' ',160,stdo,iter,fmom,amom,efermi,
     .      vnow)
          iter = -iter
        endif
        goto 10
      endif

C     call bzbfde(nspc == 2,nkp,nsp,nbmx,nevx,-vnow/2,swtk,eb)

      fmom(2) = vnow/2
      fmom(3) = amom*vnow/2
      if (lswtk == 1 .and. lwtkb == 1) then
        lwtkb = 2
      elseif (lswtk == 1 .and. lwtkb == 2) then
        lwtkb = 1
      endif

C --- Re-entry point for return: return dosef
   99 continue
      dosef = dosefl(1) + dosefl(2)

      end
      subroutine bzwtsm(lswtk,nkp,nsp,nevx,wtkb,swtk,amom)
C- Determine the magnetic moment, collinear or noncollinear case
C ----------------------------------------------------------------------
Ci Inputs
Ci   lswtk :if true, swtk is used.  Otherwise, collinear case assumed:
Ci         :swtk(*,1,*) = 1  and swtk(*,2,*) = -1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nevx  :Maximum number of bands
Ci   wtkb  :band weights
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Used when lswtk is set
Co Outputs
Co   amom  :magnetic moment
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Jun 07
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lswtk
      integer nkp,nevx,nsp
      double precision wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),amom
C ... Local parameters
      integer ikp,k
      double precision dsum

      if (nsp == 1) return

C      if (lswtk == 1 .and. nspc == 2) then
      if (lswtk) then
        amom = 0
        do  ikp = 1, nkp
          do  k = 1, nevx
            amom = amom + wtkb(k,1,ikp)*swtk(k,1,ikp)
     .                  + wtkb(k,2,ikp)*swtk(k,2,ikp)
          enddo
        enddo
      else
        amom = 0
        do  ikp = 1, nkp
          amom = amom + dsum(nevx,wtkb(1,1,ikp),1) -
     .                  dsum(nevx,wtkb(1,2,ikp),1)
        enddo
      endif
      end
      subroutine bzbfde(lswtk,nkp,nsp,nbmx,nevx,beff,swtk,eb)
C- Shift the energy bands by a given constant field
C ----------------------------------------------------------------------
Ci Inputs
Ci   lswtk :if true, swtk is used.  Otherwise, collinear case assumed:
Ci         :swtk(*,1,*) = 1  and swtk(*,2,*) = -1
Ci   nkp   :number of k-points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbmx  :leading dimension of eb
Ci   nevx  :Maximum number of bands
Ci   beff  :effective field
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Used when lswtk is set
Ci
Cio Outputs
Cio   eb   :energy bands shifted by beff (spin 1), -beff (spin 2)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Jun 07
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lswtk
      integer nkp,nsp,nbmx,nevx
      double precision beff,eb(nbmx,nsp,nkp),swtk(nevx,nsp,nkp)
C ... Local parameters
      integer ikp,ib,isp

C     if (nsp == 1) return

      do  ikp = 1, nkp
      do  ib = 1, nevx
        if (lswtk) then
          eb(ib,1,ikp) = eb(ib,1,ikp) + beff*swtk(ib,1,ikp)
          eb(ib,2,ikp) = eb(ib,2,ikp) + beff*swtk(ib,2,ikp)
        else
          do  isp = 1, nsp
            eb(ib,isp,ikp) = eb(ib,isp,ikp) + (3-2*isp)*beff
          enddo
        endif
      enddo
      enddo

      end
