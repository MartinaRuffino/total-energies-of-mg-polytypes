      subroutine mrqmintb(s_ctrl,s_lat,s_tb,npar,nterm,nlmesp,
     .  nset,nclass,nl,nsp,nspc,nkp,nelts,nfit,ldim,nvar,nvar2,mxnbr,
     .  memode,idim,ivar,ivar2,ip1,cut,cutov,iam,npm,ip2,lmx,ifit,
     .  indxsh,alam,rmaxh,efmax,delta,qp,hso,efit,sig,iwk1,wk2,iwk3,
     .  eband,epar,dpar,beta,tryp,delp,decay,deccf,decov,dcocf,tabme,
     .  tabcf,tabov,tbocf,qnu,vso,pot0,rl,par,chi0,rms0,chi,rmsdel,savc,
     .  alp,cov)
C- Performs one Levenberg-Marquardt iteration, fitting TB hamiltonian to bands
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  ltb lncol lpgf nbas nl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lqp
Cio    Passed to:  bndtbf tbham secmtb getfitptb
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat alat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos
Cio    Passed to:  bndtbf tbham
Cio  s_tb
Ci     Elts read:  alam alsc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   npar  :total number of TB parameters
Ci   nterm :number of parameters for each matrix element;
Ci         :also dimensions parameter tables.
Ci   nlmesp:(number of matrix elts/spin channel) * (# spin channels)
Ci   nset  : number of rules specifying matrix elements
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :Note: this needs to be fixed for nspc=2
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nkp   :number of k-points (bzmesh.f)
Ci   nelts :first dimension of delta
Ci   nfit  :number of eigenvalues to fit out of ldim*nkp total bands
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nvar  :number of parameters to vary out of npar total parameters
Ci   nvar2 :number of sticky parameters
Ci   mxnbr :upper limit to total number of pairs in neighbor table
Ci   memode:Defines the functional form for TB matrix elements (rdtbh.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   ivar  :(1,i) points to the position in the full list of ith variable
Ci         :(2,i) parameter type of ith variable, with types 1 to 6
Ci         :(see getfitptb)
Ci   ivar2 :(1,i) points to position in full list of ith sticky parameter
Ci         :(2,i) points to position in full list of parameter to stick to
Ci         :(see getfitptb)
Ci   ip1   : pointer to locations in full list of variables, types 1-4
Ci         :(see getfitptb)
Ci   ip2   : pointer to locations in full list of variables, types 5,6
Ci         :(see getfitptb)
Ci   ifit  :(1,i)..(2,i): range of bands to fit for ith k-point
Ci   iam   :iam(1,kk) and iam(2,kk) are the two classes for the kkth ME pair
Ci   npm   :table of offsets to iam so that every rule associated
Ci          with a given class can be easily identified (rdtbh.f).
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   rmaxh :cutoff radius for pair neighbor table
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   delta :holds electrostatic increments to the diagonal matrix
Ci   qp    :list of k-points
Ci   hso   :spin-orbit hamiltonian
Ci   efit  :eigenvalues to be fit in a single vector (see routine xcopb)
Ci   sig   :standard deviations for each eigenvalue to be fit
Ci   iwk1,wk2,iwk3: work arrays for call to gaussj, length npar
Ci   tryp  :work array for trial par vector, dimensioned same as par
Ci   delp  :work array for change in parameters, dimensioned same as par
Ci   pot0  :monopole potential at site R (eq. 7.81, Finnis)
Ci   rl    :if true, estat calculated for L>0
Cio Inputs/Outputs
Cio   alp  :curvature (alpha) matrix, see Numerical Recipes
Cio        :NB: alp must be preserved between successive calls to mrqmin
Cio   beta :beta vector (see Numerical Recipes)
Cio        :NB: beta must be preserved between successive calls to mrqmin
Cio  ... All of the following are updated in this routine:
Cio  decay,deccf,decov,dcocf: exponential or power decay parameters
Cio        :matrix element [memode = 2, v_ij d^(-b); 3, v_ij exp(-c d)]
Cio  par   :updated TB parameters in a single array
Co Outputs
Co   alam  :lambda variable, see remarks and Numerical Recipes
Co         :NB: alam must be preserved between successive calls to mrqmin
Co         :except for the final call, when alam=0
Co   eband :work array, holding energy bands for all k-points
Co   epar  :subset of eband, evals to be included in fit (all k-points)
Co   dpar  :derivatives of epar wrt each of 1..nvar TB parameters
Co   chi0  :chi-squared value for fit with initial parameters
Co   rms0  :RMS for fit with initial parameters
Co   chi   :chi-squared value for fit with current parameters,
Co         :(called chisq in Numerical Recipes)
Co   rmsdel: RMS for fit with current parameters
Co   savc  :chi-squared value from previous iteration
Co   cov   :covariance matrix, see Numerical Recipes
Co  ... All of the following TB parameters are updated in this routine:
Co   tabme :table of TB matrix elements; functional form specified by memode
Co   qnu   :diagonal TB parameters
Co   vso   :table of spin-orbit parameters
Co   tabcf :table of crystal field MEs
Co   tabov :table of overlap matrix elements, structured as tabme
Co   tbocf :table of crystal field MEs for overlap, structured as tabme
Cs Command-line switches
Cs   --DC      :
Cs   --RRR     :
Cs   --X       :
Cs   --allvecs :
Cs   --eispack :
Cs   --wvecs   :
Cs   -dumph    :
Cr Remarks
Cr   See Chapter 15 Numerical Recipes, 2nd Edition
Cr   Set alam < 0 for initialization.  Set alam = 0d0 for final call to
Cr   set covariance and curvature (alpha) matrices.
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer npar,nterm,nlmesp,nset,nclass,nl,nsp,nspc,nkp,nelts,nfit,
     .  ldim,nvar,nvar2,mxnbr,memode,idim
      integer ivar(2,npar),ivar2(2,npar),ip1(nterm+1,nlmesp,nset,4),
     .  iam(3,1),npm(2,nclass),ip2(nl,nsp,nclass,5:6),lmx(nclass),
     .  ifit(2,nkp),indxsh(1),iwk1(npar),iwk3(npar)
      double precision alam,rmaxh,efmax,chi0,rms0,chi,rmsdel,savc,
     .  wk2(npar)
      double precision delta(nelts,1),qp(3,nkp),hso(nl**2,nl**2,4,2),
     .  efit(nfit),sig(nfit),eband(ldim*nsp,nkp),epar(nfit),
     .  dpar(nvar,nfit),beta(npar),tryp(npar),delp(npar),
     .  decay(nlmesp,nset),deccf(nlmesp,nset),decov(nlmesp,nset),
     .  dcocf(nlmesp,nset),tabme(nterm,nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),tabov(nterm,nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),qnu(3,nl,nsp,nclass),vso(nl,nclass),
     .  par(npar),alp(npar,npar),cov(npar,npar),pot0(*),
     .  cut(2,nlmesp,1),cutov(2,nlmesp,1)
      logical rl
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_tb)::    s_tb
C ... Local parameters
      integer kk,j,ihit,k,iprint,i1mach,ierr
      double precision alam0,alsc,dif,xx(2)
      character*80 outs
      character*16 ptype(6)
C ... External calls
      external awrit1,awrit2,info5,bndtbf,covsrttb,dcopy,dgemm,dsidi,
     .         dsifa,getfitptb,mrqcof,rx,rxx,upack2

      data ptype /'Hamiltonian','Crystal Field','Overlap',
     .  'Overlap Crys. F.','Ham. diagonal','Ham. spin-orbit'/

      alam0 = s_tb%alam
      alsc = s_tb%alsc

C --- Initialization ---
      if (alam < 0d0) then
        kk = nvar + 1

C --- Does ivar contain a proper permutation of the parameters ---
        do  20  j = 1, npar
          ihit = 0
          do  10  k = 1, nvar
            if (ivar(1,k) == j) ihit = ihit + 1
   10     continue
          if (ihit == 0) then
            ivar(1,kk) = j
            kk = kk + 1
          elseif (ihit > 1) then
            call rx('MRQMIN: duplicate parameter in ivar')
          endif
   20   continue
        call rxx(kk /= (npar+1),'MRQMIN: Improper permutation in ivar')

C   ... Eigenvalues and derivatives for current parameter set
        call bndtbf(s_ctrl,s_lat,nterm,nlmesp,nset,nclass,nl,nsp,
     .    nspc,nkp,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,cut,cutov,iam,
     .    npm,ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,
     .    dcocf,tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,
     .    eband,epar,dpar)

C   ... Make chi, alp, beta
        call mrqcof(nfit,nvar,npar,sig,efit,epar,dpar,chi,alp,beta)

        savc = chi
        chi0 = chi
        rmsdel = dsqrt(chi / nfit)
        rms0 = rmsdel
        call dcopy(npar,par,1,tryp,1)
        alam = alam0
C   ... Printout
        call info5(30,1,0,'%N MRQMIN: Initial values: '
     .    //'chi^2=%g  rms=%g  alam=%g',chi,rmsdel,alam,0,0)
        if (iprint() >= 50) then
          write(*,500)
          do  30  j = 1, nfit
            dif = epar(j) - efit(j)
            write(*,510) j,efit(j),epar(j),dif
   30     continue
  500     format(/4x,'I',5x,'EIGENVALUE',10x,'FIT',12x,'DIFFERENCE')
  510     format(i5,2f16.10,g20.10)
        endif
      endif

C --- Alter linearized fitting matrix by augmenting diagonal elements ---
      do  50  j = 1, nvar
        do  40  k = 1, nvar
          cov(j,k) = alp(j,k)
   40   continue
        cov(j,j) = alp(j,j)*(1d0 + alam)
C       Only needed if calling gaussj below
        delp(j) = beta(j)
   50 continue

C --- Matrix solution ---
C     Solver given in Numerical Recipes
C     call gaussj(nvar,npar,1,1,iwk1,wk2,iwk3,cov,delp)
C     call prmx('mrqmin after inversion: cov',cov,npar,nvar,nvar)
C     Linpack analog
      call dsifa(cov,npar,nvar,iwk1,ierr)
      call rxx(ierr /= 0,'MRQMIN: failed to invert covariance matrix')
      call dsidi(cov,npar,nvar,iwk1,xx,j,wk2,1)
      do  k = 1, nvar
      do  j = 1, k
        cov(k,j) = cov(j,k)
      enddo
      enddo
      call dgemm('N','N',nvar,1,nvar,1d0,cov,npar,beta,npar,0d0,delp,
     .  npar)
C      call prmx('mrqmin after inversion: cov',cov,npar,nvar,nvar)
C      call prmx('mrqmin: beta',beta,npar,nvar,1)
C      call prmx('mrqmin: delp',delp,npar,nvar,1)

C --- Once converged evaluate covariance matrix with alam=0 ---
      if (alam == 0) then
        call covsrttb(nvar,npar,ivar,cov)
        return
      endif

C --- Trial parameters ---
      do  60  j = 1, nvar
        tryp(ivar(1,j)) = par(ivar(1,j)) + delp(j)
   60 continue
      do  70  j = 1, nvar2
        tryp(ivar2(1,j)) = tryp(ivar2(2,j))
   70 continue
C ... Printout
      if (iprint() >= 40) then
        write(*,520)
        do  80  j = 1, nvar
          write(*,530) ivar(1,j),par(ivar(1,j)),tryp(ivar(1,j)),delp(j),
     .      ptype(ivar(2,j))
   80   continue
  520   format(/4x,'I',4x,'OLD PARAM',8x,'TRIAL PARAM',9x,'DIFFERENCE',
     .    9x,'TYPE')
  530   format(i5,3g18.10,3x,a)
      endif
C ... Put trial parameters tryp into tables
      call getfitptb(s_ctrl,1,nterm,nlmesp,nset,nl,nsp,nclass,0,
     .  npar,memode,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,tabme,decay,tabcf,
     .  deccf,tabov,decov,tbocf,dcocf,qnu,vso,0,j,k,0,0,0,0,tryp)
C ... Eigenvalues and derivatives for current parameter set
      call bndtbf(s_ctrl,s_lat,nterm,nlmesp,nset,nclass,nl,nsp,
     .  nspc,nkp,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,cut,cutov,iam,
     .  npm,ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,
     .  dcocf,tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,
     .  eband,epar,dpar)
C ... Make chi, cov, delp
      call mrqcof(nfit,nvar,npar,sig,efit,epar,dpar,chi,cov,delp)

      call awrit2('%x%N MRQMIN: old chi^2=%g  new chi^2=%g',
     .  outs,80,0,savc,chi)

      if (chi < savc*0.9999d0) then

C   --- Success, accept new solution and decrease alam ---
        alam = alam / alsc
        savc = chi
        do  100  j = 1, nvar
          do  90  k = 1, nvar
            alp(j,k) = cov(j,k)
   90     continue
          beta(j) = delp(j)
          par(ivar(1,j)) = tryp(ivar(1,j))
  100   continue
        do  110  j = 1, nvar2
          par(ivar2(1,j)) = tryp(ivar2(2,j))
  110   continue

      else

C   --- Failure, increase alam ---
        alam = alsc*alam
        chi = savc
      endif
      rmsdel = dsqrt(chi / nfit)

C ... Printout
      if (iprint() >= 30) then
        call awrit1('%a  alam=%g',outs,-80,-i1mach(2),alam)
        if (iprint() >= 50) then
          write(*,500)
          do  120  j = 1, nfit
            dif = epar(j) - efit(j)
            write(*,510) j,efit(j),epar(j),dif
  120     continue
        endif
      endif

C --- Restore tables with parameters par ---
      call getfitptb(s_ctrl,1,nterm,nlmesp,nset,nl,nsp,nclass,0,
     .  npar,memode,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,tabme,decay,tabcf,
     .  deccf,tabov,decov,tbocf,dcocf,qnu,vso,0,j,k,0,0,0,0,par)

      end

      subroutine covsrttb(nvar,ncof,ivar,cov)
C- Repack the covariance matrix to the true order
C ----------------------------------------------------------------------
Ci Inputs
Ci   nvar: number of parameters to vary out of ncof total parameters
Ci   ncof: total number of parameters (including fixed parameters)
Ci   ivar: ivar(1,i) = index to ith parm to be varied in the cof vector
Co Outputs
Co   cov: covariance matrix, repacked to true order on output
Cr Remarks
Cr   Elements of cov associated with fixed parameters set to zero.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nvar,ncof
      integer ivar(2,nvar)
      double precision cov(ncof,ncof)
C Local parameters
      integer i,j
      double precision swap

C --- Zero all elements below diagonal ---
      do  20  j = 1, ncof-1
        do  10  i = j+1, ncof
          cov(i,j) = 0d0
   10   continue
   20 continue

C --- Repack off-diagonal elements of fit into correct locations ---
      do  40  i = 1, nvar-1
        do  30  j = i+1, nvar
          if (ivar(1,j) > ivar(1,i)) then
            cov(ivar(1,j),ivar(1,i)) = cov(i,j)
          else
            cov(ivar(1,i),ivar(1,j)) = cov(i,j)
          endif
   30   continue
   40 continue

C --- Temporarily store original diagonal elements in top row ---
      swap = cov(1,1)
      do  50  j = 1, ncof
        cov(1,j) = cov(j,j)
        cov(j,j) = 0d0
   50 continue
      cov(ivar(1,1),ivar(1,1)) = swap

C --- Sort elements into proper order on diagonal ---
      do  60  j = 2, nvar
        cov(ivar(1,j),ivar(1,j)) = cov(1,j)
   60 continue

C --- Fill in above diagonal by symmetry ---
      do  80  j = 2, ncof
        do  70  i = 1, j-1
          cov(i,j) = cov(j,i)
   70   continue
   80 continue

      end
