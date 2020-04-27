      subroutine chkme(iskipp,iovl,nl,nlmesp,nspc,ltb,nterm,tabme,
     .  decay,cut,v0,nset,memodk,ppmode,poly,cutmod,cutpp)
C- Consistency check for data read in ME section
C ----------------------------------------------------------------------
Ci Inputs
Ci   iskipp: 0 do all checks, otherwise skip the PP part
Ci   iovl  : 0 if tabme is a Hamiltonian, 1 if tabme is an overlap matrix
Ci           (ME signs are opposite)
Ci   nl    : global max lmax+1
Ci   nlmesp: (number of ME/spin channel) * (# spin channels)
Ci           for which ME are input
Ci   nspc  : number of spin channels (spin-orbit)
Ci   ltb   : TB switches
Ci   nterm : max number of parameters for each matrix element,
Ci           leading dimension of tabme and a few local auto arrays
Ci   v0    : pair potential parameters
Ci   tabme : a set of parameters that correspond to coefficients of
Ci           Slater-Koster, or Slater-Koster-like matrix elements.
Ci           The meaning of the coefficients, depends on memode
Ci           (see rdtbh)
Ci   decay : exponential or power decay parameters
Ci           matrix element [memode = 2, v_ij d^(-b); 3, v_ij exp(-c d)]
Ci   cut   : For each ME a pair of distances (r1,rc). If cutmod .ne 0,
Ci           the ME is cut-off smoothly to zero between r1 and rc
Ci           (see makvme)
Ci   v0    : parameters for pair potential
Ci           Ordering: a_1 b_1 c_1 a_2 ... c_3
Ci   nset  : number of ME rules specified for each pair of species
Ci   memodk: ME rule for each of the pair of species
Ci           0, fixed MEs
Ci           1, Harrison universal MEs
Ci           2, exponential decay
Ci           3, power decay
Ci           4, ME = \sum_i=1,3 a_i d^b_i exp(-c_i d), the ordering is:
Ci              a_1 b_1 c_1 a_2 ... c_3 for ss-sigma, then sp-sigma, etc
Ci           5, Goodwin-Skinner-Pettifor,
Ci              v_ij (r0/d)^n exp[n (-{d/rc}^nc + {r0/rc}^nc)]
Ci              ordering: v_ij n nc r0 rc for ss-sigma, etc
Ci           6, (nl=1) Slater-Koster + Srini's extra term
Ci              NB: NOT implemented for spin pol.
Ci           7, a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc)
Ci              ordering: a b c d0 for ss-sigma, etc
Ci              nterm=1 for memode=0-3, nterm=9 for memode=4,
Ci              nterm=5 for memode=5, nterm=2 for memode=6,
Ci              nterm=4 for memode=7
Ci              memode >= 10, use canonical TB Hamiltonian
Ci   ppmode: type of repulsive pair potential
Ci           10s digit:
Ci           1, sum of exponents times power law decay
Ci              V0(d) = \sum_i=1,3 a_i d^b_i exp(-c_i d)
Ci           2, quadratic
Ci              V0(d) = a_1 eps + a_2 eps^2 where eps = (d - c_1) / c_1
Ci           3,  Goodwin-Skinner-Pettifor
Ci              V0(d) = A (r0/d)^m exp[m (-{d/rc}^mc + {r0/rc}^mc)]
Ci           1s digit defines particular adjustments to above types (see makvpp.f)
Ci   poly  : degree of polynomial Pn employed to smoothly cut off V0 to zero
Ci           4  cutoff tail matches value, slope, and curvature of V0(r) at r=r1
Ci              whereas only value and slope are zero at r=rc
Ci           5  (default) same as above except that value, slope, and curvature
Ci              of Pn(r1) should all turn to zero
Ci   cutmod: cutoff mode for pair potentials, hopping integrals, and overlap matrix
Ci           (if applicable)
Ci           0  no cutoff
Ci           1  augmentative cutoff: V0(r) is augmented with a polynomial P_n(r)
Ci              at r = r1 up to the second derivatives, whereas P_n(rc),
Ci              P'_n(rc) and P"_n(rc) (if poly = 5) are all zeros
Ci           2  multiplicative cutoff (default mode): V0(r) is multiplied by
Ci              a polynomial P_n(r): P_n(r1)=1 and P_n(rc)=P'_n(rc)=P'_n(r1)=P"_n(r1)=0.
Ci              If poly=5, also P"_n(rc)=0
Ci   cutpp : pair potential cutoffs cutpp(1,*) = r1 and cutpp(2,*) = rc, see above.
Ci           (previously contained within v0)
Ci   memodk, ppmode, poly, cutmod, and cutpp are all arrays defined for each pair of species
Co Outputs
Co   None unless an inconsistency is detected. In which case the program stops and prints
Co   the respective diagnostics.
Cr Remarks
Cr   The program stops after the first inconsistency is encountered. Hence it may take a few
Cr   runs to fully satisfy chkme. Use --nochkme switch if you wish to disable chkme.
Cr
Cr   Signs for Hamiltonian and overlap matrix elements are opposite. Switch iovl takes care
Cr   of this.
Cu Updates
Cu   18 Apr 11 (SL) first created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer, intent(in) :: iskipp,iovl,nl,nlmesp,nspc,nterm,nset,ltb
      integer, intent(in) :: memodk(nset),ppmode(nset),
     .                       poly(nset),cutmod(nset)
      double precision, intent(in) :: tabme(nterm,nlmesp,nset),
     .  decay(nlmesp,nset),cut(2,nlmesp,nset),v0(9,nset),cutpp(2,nset)
C Local parameters
      integer memx,ii,memode,lgunit
      logical cryf,ovl,ocryf
c ... ME modes
      integer, parameter :: memodx = 7      ! max ME mode number
      integer ntm(0:memodx)
c        memode 0 1 2 3 4 5 6 7 - - -
      data ntm /1,1,1,1,9,5,2,4/
      parameter (memx=500)
      logical bittst
      integer icutm,ilme,ierr,mode0,mode1,iexp,iprint
      integer itabsp(2,2)
      data itabsp /0,2,2,1/
      integer nlme(3)
      data nlme /1,4,10/
      integer isp1,isp2,ioff,il,im,chksgn
      double precision vpp(3,3),agsp,m,mc,r0,rc
      integer, parameter :: ilmex=10
      double precision vme(nterm,ilmex),avme(ilmex)
      double precision nme(ilmex),ncme(ilmex)
      double precision r0me(ilmex),rcme(ilmex)
c ... 'acceptance' tolerance
      double precision, parameter :: tol=1d-8
c ... arrays for printout
      character pprtt(0:3)*40,pprtta(0:3)*40,pprt(0:3)*80,pprta(0:3)*80,
     .          mprtt(0:7)*40,psig(10)*3,vmesgn(10)*3,pham(0:1)*11
      data vmesgn /' - ',' + ',' + ',' - ',' - ',
     .             ' - ',' + ',' - ',' + ',' - '/
      data psig   /'sss','sps','pps','ppp','sds',
     .             'pds','pdp','dds','ddp','ddd'/
      data pprtt/'no pair interaction','Exp x Power law','Quadratic',
     .           'GSP'/
      data mprtt /'Fixed','Universal ME','Exp. decay','Power decay',
     .  'Exp. x Power law', 'GSP','Fixed + extension',
     .  'Fermi x Power law'/
      data pprtt/'no pair interaction','Exp x Power law','Quadratic',
     .           'GSP'/
      data pprtta/'',' + Exp','',''/
      data pprt/'',
     .  ', V0(d) = sum [a d^b exp(-c d)]',
     .  ', V0(d) = a1 eps + a2 eps^2, eps = d/d0 - 1',
     .  ', V0(d) = A (r0/d)^m exp{m [(r0/rc)^mc - (d/rc)^mc]}'/
      data pprta/'',' + a exp(-c d)','',''/
      data pham/'Hamiltonian','overlap'/

C --- Initialization ---
      ovl   = bittst(ltb,1)
      cryf  = bittst(ltb,2)
      ocryf = bittst(ltb,4)
      ilme = nlme(nl)
      if (ilme*(2*nspc-1) /= nlmesp) then
        call info5(10,0,0,'chkme: ilme = %i nspc = %i and nlmesp = %i '
     .                     //'do not match',ilme,nspc,nlmesp,0,0)
        call rx('chkme: inconsistent dimensions')
      endif

C --- nset = 0 => nothing to check ---
      if (nset == 0) then
        call info0(30,0,0,
     .    'chkme: Warning! no rules specified for either PP or ME')
        return
      endif

C --- check cutoff switches ---
      do ii = 1, nset
        icutm = cutmod(ii)
        call sanrg(.true.,icutm,0,2,'chkme:','cutmod')
        if (icutm /= 0)
     .    call sanrg(.true.,poly(ii),4,5,'chkme:','poly')
      enddo


      if (iskipp == 0) then
C --- check pair potentials ---
        do ii = 1, nset
          ierr = 0
          mode0 = mod(ppmode(ii),10)
          mode1 = mod(ppmode(ii)/10,10)
          vpp(1:3,1:3) = 0d0
          call dcopy(9,v0(1,ii),1,vpp,1)
          select case (mode1)
            case (0)
c ...       no PP
              call sanrg(.true.,mode0,0,0,'chkme:','ppmode')
            case (1)
c ...       sum of exponentials times power law decay
              call sanrg(.true.,mode0,0,0,'chkme:','ppmode')
              if ((min(vpp(1,1),vpp(1,2),vpp(1,3)) < -tol) .and.
     .            (max(vpp(1,1),vpp(1,2),vpp(1,3)) < tol)) ierr = 10
              do iexp = 1, 3
                if (vpp(2,iexp) > tol .or. vpp(3,iexp) < -tol)
     .          ierr = 10
              enddo
            case (2)
c ...       quadratic PP
              call sanrg(.true.,mode0,0,1,'chkme:','ppmode')
              if (mode0 == 0) then
c ...         cutoff option should be on in Jim Chadi potential
                if (cutmod(ii) == 0) ierr = 20
              else
                if ((min(vpp(1,1),vpp(3,1)) < -tol) .and.
     .              (max(vpp(1,1),vpp(3,1)) < tol)) ierr = 21
                if (vpp(2,1) > tol .or. vpp(1,2) > tol) ierr = 21
                if (abs(vpp(3,2)) >= tol .and. vpp(2,2) <= 0d0)
     .             ierr = 21
              endif
            case (3)
c ...       Goodwin-Skinner-Pettifor PP
              call sanrg(.true.,mode0,0,3,'chkme:','ppmode')
              agsp = vpp(1,1)
              m  = vpp(1,2)
              mc = vpp(2,2)
              r0 = vpp(3,2)
              rc = vpp(1,3)
c ...         agsp, m, mc, r0, rc should all be positive, rc /= 0, and rc >= r0
              if ((min(agsp, m, mc, r0) < -tol)
     .         .or. ((rc-r0) < -tol) .or. (rc <= 0d0)) ierr = 30
c ...         check extra terms in GSP (mode0 /= 0)
              if (mode0 == 1) then
c ...           GSP + v0(2,3)*dexp(-v0(3,3)*d), v0(3,3) >= 0
                if (vpp(3,3) < -tol) ierr = 31
              elseif (mode0 == 2) then
c ...           GSP + v0(2,3)*d**v0(3,3), v0(3,3) <= 0
                if (vpp(3,3) > tol) ierr = 32
              endif
            case default
              call rxi('chkme: PP mode does not exist, ppmode = ',
     .                 ppmode(ii))
          end select
          if (ierr /= 0) then
            call info2(10,1,0,' chkme: Problem with PP detected'//
     .                        ' for set %i, ierr = %i',ii,ierr)
            if (iprint() >= 10) then
              call awrit1('  PP type %i ('//trim(pprtt(mode1))//
     .        trim(pprtta(mode0))//')'//trim(pprt(mode1))//
     .        trim(pprta(mode0)),' ',120,lgunit(1),ppmode(ii))
              if (ierr == 1) then
                call awrit1('  V0 = %9:1d ',' ',120,lgunit(1),vpp)
              else
                call awrit1('  Quadratic PP requires cutoff,'//
     .            ' cutmod = %i',' ',120,lgunit(1),cutmod(ii))
                call awrit1('  V0 = %9:1d ',' ',120,lgunit(1),vpp)
              endif
            endif
            call rx('chkme: bad PP')
          endif
c ...     check cutoff distances
          if (cutmod(ii) /= 0) then
            r0 = cutpp(1,ii)
            if (r0 < -tol .or. cutpp(2,ii) < r0) then
              call info2(10,0,0,
     .          ' Check pair potential cutoff for set %i, cutpp = %2d',
     .        ii,cutpp(1,ii))
              call rx('chkme: bad PP cutoff')
            endif
          endif
        enddo
      endif

C --- check ME (hopping integrals) ---

C ... Loop over spin combinations (spin-orbit)
      do isp1 = 1, nspc
      do isp2 = 1, nspc
        ioff = 1 + itabsp(isp1,isp2)*ilme

        do ii = 1, nset
          ierr = 0
          vme(1:nterm,1:ilme) = 0d0
          memode = memodk(ii)
          if (memode /= 1) then
            call dcopy(nterm*ilme,tabme(1,ioff,ii),1,vme,1)
c ...       check sign of the prefactor
            avme(1:ilme) = vme(1,1:ilme)
            ierr = chksgn(iovl,ilme,avme)
          endif
          select case (memode)
            case (0, 6)
C ...       fixed MEs
              if (memode == 6 .and. nl /= 2)
     .          call rx(' chkme: need nl = 2 for memode 6')
            case (1)
C ...       Harrison universal MEs, check switches
              call rxx(ovl,
     .        'chkme: No overlap matrix for universal Hamiltonian')
              call rxx((cryf .or. ocryf),
     .        'chkme: No crystal field terms for universal Hamiltonian')
            case (2, 3)
C ...       exponential or power law decay
              do il = 1, ilme
                if (decay(il-1+ioff,ii) < -tol) ierr = 2
              enddo
            case (4)
c ...       sum of exponentials times power law decay
              do il = 1, ilme
                do im = 1, 3
                  if  (vme(2+(im-1)*3,il) > tol
     .            .or. vme(3+(im-1)*3,il) < -tol) ierr = 2
                enddo
              enddo
            case (5)
c ...       Goodwin-Skinner-Pettifor
              nme(1:ilme)  = vme(2,1:ilme)
              ncme(1:ilme) = vme(3,1:ilme)
              r0me(1:ilme) = vme(4,1:ilme)
              rcme(1:ilme) = vme(5,1:ilme)
c ...         n, nc, r0, rc should all be positive and rc >= r0
              do il = 1, ilme
                if (min(nme(il),ncme(il),r0me(il),rcme(il)) < -tol
     .           .or. (rcme(il)-r0me(il)) < -tol) ierr = 2
              enddo
            case (7)
c ...       ME = a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc)
              nme(1:ilme)  = vme(2,1:ilme)
              do il = 1, ilme
                if (nme(il) < -tol) ierr = 2
              enddo
            case default
              call rxi('chkme: ME mode does not exist, memode = ',
     .          memode)
          end select
          if (ierr /= 0) then
            call info2(10,1,0,' chkme: Problem with '//trim(pham(iovl))
     .                 //' ME detected for set %i, ierr = %i',ii,ierr)
            if (nspc /= 1)
     .      call info2(10,0,0,'        isp1 = %i, isp2 = %i',isp1,isp2)
            if (iprint() >= 10) then
              print '(a,i4,3a)',' memode = ',memode,
     .          ' (',trim(mprtt(memode)),')'
              print '(14x,10(a,7x))',psig(1:ilme)
              if (ierr. eq. 1)
     .        print '(14x,10(a,7x))',vmesgn(1:ilme)
              print *,' vme = '
              do im = 1, ntm(memode)
                print '(11x,10g10.2)',vme(im,1:ilme)
              enddo
              if (memode == 2 .or. memode == 3) then
                print *,' DECAY = '
                print '(11x,10g10.2)',decay(ioff:ioff+ilme-1,ii)
              endif
            endif
            call rxi('chkme: bad ME for set ',ii)
          endif

c ...   check cutoff distances unless ME are fixed
          if (cutmod(ii) /= 0 .and.
     .        memode /= 0 .and. memode /= 6) then
            do il = 1, ilme
              r0 = cut(1,il-1+ioff,ii)
              rc = cut(2,il-1+ioff,ii)
              if (r0 < -tol .or. rc < r0) then
                call info5(10,0,0,' Check '//trim(pham(iovl))//
     .            ' ME cutoff for set %i, rule = %i, cutme = %2d',
     .            ii,il,cut(1,il-1+ioff,ii),0,0)
                call rx('chkme: bad ME cutoff')
              endif
            enddo
          endif
        enddo

      enddo
      enddo


      end

      integer function chksgn(iovl,n,v)
C- check if the sign of hopping integrals V is correct
C ----------------------------------------------------------------------
Ci Inputs
Ci   iovl  : 0 if v is a Hamiltonian, 1 if v is an overlap matrix
Ci           (ME signs are opposite)
Ci   n     : number of hopping integrals to check
Ci   v     : hopping integrals
Co Outputs
Co  chksgn :0 if all signs are correct, 1 otherwise
C ----------------------------------------------------------------------
      implicit none
C ... Passed variables
      integer, intent(in) :: iovl,n
      double precision, intent(in) :: v(n)
C ... Local variables
      integer, parameter :: nmax=10
      double precision, parameter :: tol=1d-8
      integer vsign(nmax),i,mone
c                 sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd
      data vsign /-1,   1,  1, -1, -1, -1,  1, -1,  1, -1/

      mone = float(1)
      if (iovl /= 0) mone = -mone

      chksgn = 0
      do i = 1, min(n,nmax)
        if (v(i)*vsign(i)*mone < -tol) chksgn = 1
      enddo

      end







