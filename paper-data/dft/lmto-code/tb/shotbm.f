      subroutine shotbm(lscale,alat,tabme,decay,v0,nset,npair,nterm,
     .  nl,nsp,nlmesp,memode,ppmode,poly,cutmod,cutpp,cutme,
     .  nclass,dclabl,iam,ipair,strn)
C- Print table of tight-binding matrix elements and pair potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   lscale : scale (lscale=.true.)/do not scale (lscale=.false.)
Ci            cutoff distances cutme and cutpp with alat
Ci   alat   : lattice constant
Ci   tabme,decay,v0,nterm,memode,iam: see subroutine rdtbh
Ci   nset   : number of matrix element sets
Ci   npair  : number of matrix element pairs
Ci   nl,nsp,nlmesp,nclass,dclabl
Ci   ipair  : 1, print pair potential; 0, do not
Ci   strn   : identifier string for MEs (printed)
Cu Updates
Cu   19 Apr 11 (SL)  species-dependent memode and cutoffs
Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical, intent(in) :: lscale
      integer, intent(in) :: nset,npair,nterm,nl,nsp,nlmesp,nclass,ipair
      integer, intent(in) ::  memode(nset),ppmode(nset),cutmod(nset),
     .                        poly(nset)
      integer, intent(in) ::  iam(3,*)
      double precision, intent(in) ::  alat,tabme(nterm,nlmesp,nset),
     .   v0(3,3,nset),decay(nlmesp,nset),
     .   dclabl(nclass),cutpp(2,nset),cutme(2,nlmesp,nset)
      character*(*), intent(in) ::  strn
C Local variables
      integer ip,np,ic,jc,k,l,m,nme,lgunit
      integer memod0,ppmod0,ppmod1,cutm
      double precision wk(10),ascale
      character clabli*8,clablj*8,strng*40
      character pprtt(0:3)*40,pprtta(0:3)*40,pprt(0:3)*80,pprta(0:3)*80,
     .          cutprt(0:2)*40,mprtt(0:7)*40,mprt(0:7)*80,psig(10)*3
      data psig/'sss','sps','pps','ppp','sds',
     .          'pds','pdp','dds','ddp','ddd'/
      data pprtt/'no pair interaction','Exp x Power law','Quadratic',
     .           'GSP'/
      data mprtt /'Fixed','Universal ME','Exp. decay','Power law decay',
     .  'Exp. x Power law', 'GSP','Fixed + extension',
     .  'Fermi x Power law'/
      data mprt/'',', V(d) = a / d^2',
     .  ', V(d) = a exp (- b d)',', V(d) = a d^-b',
     .  ', V(d) = sum [a d^b exp(-c d)]',
     .  ', V(d) = A (r0/d)^n exp{n [(r0/rc)^nc - (d/rc)^nc]}','',
     .  ', V(d) = a d^-b / {1+exp[c(d-d0)]}'/

      data pprtt/'no pair interaction','Exp x Power law','Quadratic',
     .           'GSP'/
      data pprtta/'',' + Exp','',''/
      data pprt/'',
     .  ', V0(d) = sum [a d^b exp(-c d)]',
     .  ', V0(d) = a1 eps + a2 eps^2, eps = d/d0 - 1',
     .  ', V0(d) = A (r0/d)^m exp{m [(r0/rc)^mc - (d/rc)^mc]}'/
      data pprta/'',' + a exp(-c d)','',''/
      data cutprt/'no cutoff','augmentative','multiplicative'/
      integer nlme(3)
      data nlme /1,4,10/

      if (nset == 0) return
c...  scaling factor
      if (lscale) then
        ascale = alat
      else
        ascale = 1d0
      endif

C --- Print pair potential ---
      if (ipair == 1) then
        call awrit1('%N SHOTBM: %i rules for pair potential',' ',80,
     .    lgunit(1),nset)
        do k = 1, nset
          np = 0
          do ip = 1, npair
C       ... Find the rule connecting this pair
            if (iam(3,ip) == k) then
              np = np + 1
              if (np == 1) then
                ic = iam(1,ip)
                jc = iam(2,ip)
              endif
            endif
          enddo
          if (np == 0) then
            call awrit1(' ... no pairs matching rule %i',' ',80,
     .        lgunit(1),k)
            cycle
          endif
          ppmod0 = mod(ppmode(k),10)
          ppmod1 = mod(ppmode(k)/10,10)
          call r8tos8(dclabl(ic),clabli)
          call r8tos8(dclabl(jc),clablj)
C     ... Print pp mode
          call awrit2(' '//clabli//'%a,'//clablj//
     .      '%a:%?#n>1# (%i more)#%j#',' ',120,lgunit(1),np,np-1)
          if (ppmode(k) == 21) then
            call awrit1('   type %i (Polynomial + Power law),'
     .        ,' ',120,lgunit(1),ppmode(k))
            call awrit0('           V0(d) = '//
     .        'a1 d^b1 + a2 d^b2 + P_3(eps) eps^3 theta(eps < 0)'
     .        ,' ',120,lgunit(1))
            call awrit0('           '//
     .        'where P_3(eps) = p0 (1 + p1 eps + p2 eps^2 + p3 eps^3)'//
     .        ', eps = d/rc - 1',' ',120,lgunit(1))
          else
            call awrit1('   type %i ('//trim(pprtt(ppmod1))//
     .        trim(pprtta(ppmod0))//')'//trim(pprt(ppmod1))//
     .        trim(pprta(ppmod0)),' ',120,lgunit(1),ppmode(k))
          endif

C     ... Print pp coefficients
          if (ppmod1 == 0) then
            cycle
          elseif (ppmod1 == 1) then
c           call awrit3('   coeffs: a1= %;8F b1= %:-1;5#7g c1= %;5g',
c    .        ' ',120,lgunit(1),v0(1,1,k),v0(2,1,k),v0(3,1,k))
          call awrit3('   coeffs: a1= %:-1;5#9g b1= %:-1;5#7g c1= %;5g',
     .        ' ',120,lgunit(1),v0(1,1,k),v0(2,1,k),v0(3,1,k))
          call awrit3('           a2= %:-1;5#9g b2= %:-1;5#7g c2= %;5g',
     .        ' ',120,lgunit(1),v0(1,2,k),v0(2,2,k),v0(3,2,k))
          call awrit3('           a3= %:-1;5#9g b3= %:-1;5#7g c3= %;5g',
     .        ' ',120,lgunit(1),v0(1,3,k),v0(2,3,k),v0(3,3,k))
          elseif (ppmod1 == 2) then
            if (ppmod0 == 0) then
              call awrit3('   coeffs: a1= %;5g a2= %;5g d0= %;5g',
     .          ' ',120,lgunit(1),v0(1,1,k),v0(1,2,k),v0(3,1,k))
            elseif (ppmod0 == 1) then
              call awrit4('   coeffs: a1= %:-1;5#9g b1= %:-1;5#7g '//
     .          'a2= %:-1;5#9g b2= %:-1;5#7g', ' ',120,lgunit(1),
     .          v0(1,1,k),v0(2,1,k),v0(3,1,k),v0(1,2,k))
              call awrit5('           rc= %:-1;4#5g p0= %:-1;5#7g '//
     .          'p1=%:-1;5#7g p2=%:-1;5#7g p3=%:-1;5#7g',
     .          ' ',120,lgunit(1),
     .          v0(2,2,k),v0(3,2,k),v0(1,3,k),v0(2,3,k),v0(3,3,k))
            endif
          elseif (ppmod1 == 3) then
            call awrit8('   coeffs: A= %;5g m= %;5g mc= %;5g'//
     .        ' r0= %;5g rc= %;5g %?@n==1@ a= %;5g c= %;5g@@',
     .        ' ',120,lgunit(1),v0(1,1,k),v0(1,2,k),v0(2,2,k),
     .        v0(3,2,k),v0(1,3,k),ppmod0,v0(2,3,k),v0(3,3,k))
          else
            call rxi(' shotbm: ppmode = %i does not exist',ppmode(k))
          endif

C     ... Print cutoff mode and cutoff distances
          cutm = cutmod(k)
          if (cutm >= 0 .and. cutm <= 2) then
            call awrit5('   cutoff type %i ('//trim(cutprt(cutm))//
     .      ')%?@n==0@@, %ith order polynomial, range [%;5g,'//
     .      ' %;5g]@',' ',120,lgunit(1),cutm,cutm,poly(k),
     .      cutpp(1,k)*ascale,cutpp(2,k)*ascale)
          else
            call rxi(' shotbm: cutmod = %i does not exist',cutm)
          endif
        enddo
      endif

  100 format(
     .  24x,'a1',4x,'b1',4x,'c1',6x,'a2',4x,'b2',4x,'c2',6x,'a3',
     .   4x,'b3',4x,'c3')
C 110 format(1x,a,',',a,':',:,'  (',i4,' more)')
  115 format('   A(r0/d)^m e^...:',1x,'A=',f7.3,2x,'m=',f6.2,
     .  2x,'mc=',f6.2,2x,'r0=',f7.3,2x,'rc=',f7.3,2x,'GSP')
  120 format('   a1 ep + a2 ep^2:',f8.3,12x,f8.3,16x,'d0=',f7.3)
  130 format('   a d^b exp(-c d):',3(f8.3,f6.2,f6.3))

c     if (memod0 >= 10) return

C --- Print MEs ---
      strng = strn
      call awrit1('%N SHOTBM: %i rules for '//strng//'%a',' ',80,
     .  lgunit(1),nset)
c     if (memod0 == 4) then
c     elseif (memod0 == 5) then
c     elseif (memod0 == 6) then
c     elseif (memod0 == 7) then
c     else
c     endif
      nme = nlme(nl)
      do  k = 1, nset
        np = 0
        do  ip = 1, npair
          if (iam(3,ip) == k) then
            np = np + 1
            if (np == 1) then
              ic = iam(1,ip)
              jc = iam(2,ip)
            endif
          endif
        enddo
        if (np == 0) then
          call awrit1(' ... no pairs matching rule %i',' ',80,
     .      lgunit(1),k)
          cycle
        endif
        if (memode(k) >= 10) then
          call awrit1('   ME mode %i >= 10',' ',80,lgunit(1),memode(k))
          cycle
        endif
        memod0 = mod(memode(k),10)
        call r8tos8(dclabl(ic),clabli)
        call r8tos8(dclabl(jc),clablj)
        call awrit2(' '//clabli//'%a,'//clablj//
     .    '%a:%?#n>1# (%i more)',' ',80,lgunit(1),np,np-1)
C     ... Print ME mode
          call awrit1('   type %i ('//trim(mprtt(memod0))//
     .      ')'//trim(mprt(memod0)),' ',120,lgunit(1),memod0)
        if (memod0 <= 3) then
          print '(9x,10(4x,a3))',psig(1:nme)
c         if (nl >= 3) then
c           print 160
c         elseif (nl >= 2) then
c           print 161
c         else
c           print 162
c         endif
        endif
        if (memod0 == 4) then
          print 150
          write (*,170) (tabme(m,1,k),m=1,9)
          if (nsp == 2) then
            write (*,180) (tabme(m,1+1*nme,k),m=1,9),
     .                    (tabme(m,1+2*nme,k),m=1,9)
          endif
          if (nl >= 2) then
            if (nsp == 2) then
              write (*,190) ((tabme(m,l+0*nme,k),m=1,9),
     .                       (tabme(m,l+1*nme,k),m=1,9),
     .                       (tabme(m,l+2*nme,k),m=1,9),l=2,4)
            else
              write (*,200) ((tabme(m,l,k),m=1,9),l=2,4)
            endif
          endif
          if (nl >= 3) then
            if (nsp == 2) then
              write (*,210) ((tabme(m,l+0*nme,k),m=1,9),
     .                       (tabme(m,l+1*nme,k),m=1,9),
     .                       (tabme(m,l+2*nme,k),m=1,9),l=5,10)
            else
              write (*,220) ((tabme(m,l,k),m=1,9),l=5,10)
            endif
          endif
        elseif (memod0 == 5) then
          print 155
          write (*,175) (tabme(m,1,k),m=1,5)
          if (nsp == 2) then
            write (*,185) (tabme(m,1+1*nme,k),m=1,5),
     .                    (tabme(m,1+2*nme,k),m=1,5)
          endif
          if (nl >= 2) then
            if (nsp == 2) then
              write (*,195) ((tabme(m,l+0*nme,k),m=1,5),
     .                       (tabme(m,l+1*nme,k),m=1,5),
     .                       (tabme(m,l+2*nme,k),m=1,5),l=2,4)
            else
              write (*,205) ((tabme(m,l,k),m=1,5),l=2,4)
            endif
          endif
          if (nl >= 3) then
            if (nsp == 2) then
              write (*,215) ((tabme(m,l+0*nme,k),m=1,5),
     .                       (tabme(m,l+1*nme,k),m=1,5),
     .                       (tabme(m,l+2*nme,k),m=1,5),l=5,10)
            else
              write (*,225) ((tabme(m,l,k),m=1,5),l=5,10)
            endif
          endif
        elseif (memod0 == 7) then
          print 157
          write (*,177) (tabme(m,1,k),m=1,4)
          if (nsp == 2) then
            write (*,187) (tabme(m,1+1*nme,k),m=1,4),
     .                    (tabme(m,1+2*nme,k),m=1,4)
          endif
          if (nl >= 2) then
            if (nsp == 2) then
              write (*,197) ((tabme(m,l+0*nme,k),m=1,4),
     .                       (tabme(m,l+1*nme,k),m=1,4),
     .                       (tabme(m,l+2*nme,k),m=1,4),l=2,4)
            else
              write (*,207) ((tabme(m,l,k),m=1,4),l=2,4)
            endif
          endif
          if (nl >= 3) then
            if (nsp == 2) then
              write (*,217) ((tabme(m,l+0*nme,k),m=1,4),
     .                       (tabme(m,l+1*nme,k),m=1,4),
     .                       (tabme(m,l+2*nme,k),m=1,4),l=5,10)
            else
              write (*,227) ((tabme(m,l,k),m=1,4),l=5,10)
            endif
            endif
        elseif (memod0 == 1) then
          call ume(nme,1d0,wk)
          write (*,230) (wk(l),l=1,nme)
          if (nsp == 2) then
            write (*,240) (wk(l),l=1,nme)
            write (*,250) (wk(l),l=1,nme)
          endif
        elseif (memod0 == 6) then
          print 156
          write (*,230) (tabme(1,l,k),l=1,8)
          if (nsp == 2) then
            write (*,240) (tabme(1,l+1*8,k),l=1,8)
            write (*,250) (tabme(1,l+2*8,k),l=1,8)
          endif
        else
          write (*,230) (tabme(1,l,k),l=1,nme)
          if (nsp == 2) then
            write (*,240) (tabme(1,l+1*nme,k),l=1,nme)
            write (*,250) (tabme(1,l+2*nme,k),l=1,nme)
          endif
          if (memod0 == 2 .or. memod0 == 3
     .      .or. memod0 == 7) then
            write (*,260) (decay(l,k),l=1,nme)
            if (nsp == 2) then
              write (*,240) (decay(l+1*nme,k),l=1,nme)
              write (*,250) (decay(l+2*nme,k),l=1,nme)
            endif
          endif
        endif
C     ... Print cutoff mode and cutoff distances
        cutm = cutmod(k)
        if (cutm >= 0 .and. cutm <= 2) then
          call awrit3('   cutoff type %i ('//trim(cutprt(cutm))//
     .    ')%?@n==0@@, %ith order polynomial, range [r1, rc]',
     .    ' ',120,lgunit(1),cutm,cutm,poly(k))
          if (cutm /= 0) then
            print '(9x,10(4x,a3))',psig(1:nme)
            write (*,270) ascale*cutme(1,1:nme,k)
            write (*,280) ascale*cutme(2,1:nme,k)
          endif
        else
          call rxi(' shotbm: cutmod = %i does not exist',cutm)
        endif
      enddo

* 140 format(/' SHOTBM: found ',i4,' sets of ',a,' matrix elements:')
  150 format(13x,'a1',5x,'b1',5x,'c1',7x,'a2',5x,'b2',5x,'c2',7x,'a3',
     .   5x,'b3',5x,'c3')
  155 format(13x,'V',8x,'n',7x,'nc',7x,'r0',7x,'rc')
  156 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp',4x,'usp',4x,'uxy',3x,
     .       'upps',3x,'uppp')
  157 format(13x,'a',8x,'b',8x,'c',7x,'d0')
  160 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp',4x,'sds',4x,'pds',4x,
     .       'pdp',4x,'dds',4x,'ddp',4x,'ddd')
  161 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp')
  162 format(13x,'sss')
C  163 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp',4x,'usp',4x,'uxy',
C     .       4x,'upps',4x,'uppp')
  170 format('   sss:',3(f9.3,f7.2,f7.3))
  175 format('   sss:',5f9.3)
  177 format('   sss:',4f9.3)
  180 format('    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3))
  185 format('    --:',5f9.3/'    +-:',5f9.3)
  187 format('    --:',4f9.3/'    +-:',4f9.3)
  190 format('   sps:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   pps:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3)/
     .       '   ppp:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3))
  195 format('   sps:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   pps:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3/
     .       '   ppp:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3)
  197 format('   sps:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   pps:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3/
     .       '   ppp:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3)
  200 format('   sps:',3(f9.3,f7.2,f7.3)/'   pps:',3(f9.3,f7.2,f7.3)/
     .       '   ppp:',3(f9.3,f7.2,f7.3))
  205 format('   sps:',5f9.3/'   pps:',5f9.3/
     .       '   ppp:',5f9.3)
  207 format('   sps:',4f9.3/'   pps:',4f9.3/
     .       '   ppp:',4f9.3)
  210 format('   sds:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   pds:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3)/
     .       '   pdp:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   dds:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3)/
     .       '   ddp:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   ddd:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3))
  215 format('   sds:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   pds:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3/
     .       '   pdp:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   dds:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3/
     .       '   ddp:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   ddd:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3)
  217 format('   sds:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   pds:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3/
     .       '   pdp:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   dds:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3/
     .       '   ddp:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   ddd:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3)
  220 format('   sds:',3(f9.3,f7.2,f7.3)/'   pds:',3(f9.3,f7.2,f7.3)/
     .       '   pdp:',3(f9.3,f7.2,f7.3)/'   dds:',3(f9.3,f7.2,f7.3)/
     .       '   ddp:',3(f9.3,f7.2,f7.3)/'   ddd:',3(f9.3,f7.2,f7.3))
  225 format('   sds:',5f9.3/'   pds:',5f9.3/
     .       '   pdp:',5f9.3/'   dds:',5f9.3/
     .       '   ddp:',5f9.3/'   ddd:',5f9.3)
  227 format('   sds:',4f9.3/'   pds:',4f9.3/
     .       '   pdp:',4f9.3/'   dds:',4f9.3/
     .       '   ddp:',4f9.3/'   ddd:',4f9.3)
  230 format('   coeff:',10f7.2)
  240 format('      --:',10f7.2)
  250 format('      +-:',10f7.2)
  260 format('   decay:',10f7.2)
  270 format('   r1:   ',10f7.2)
  280 format('   rc:   ',10f7.2)

      end
