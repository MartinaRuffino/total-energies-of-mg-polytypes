      subroutine fpchk(s_spec,s_site)
C- Routines to check various quantities in FP code
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name rmt rsma lmxa kmxt lmxb pz orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  chkxpn uspecb
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  chkxpn
Ci Inputs
Cu Updates
Cu   23 Apr 02 Added option (mode=0) to find MT radii
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)

      call info(5,2,0,
     .   ' ----- Check accuracy of augmentation expansion -----',0,0)
      call chkxpn(s_site,s_spec)
      end

      subroutine chkxpn(s_site,s_spec)
C- Check accuracy of Pkl expansions of augmented functions
Cu   05 Jul 13 Replace f77 pointers with f90 ones
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: ips(:)
C ... Local parameters
      integer k0,nlm0,nrx,n0,nkap0,stdo,lgunit
      parameter (k0=20, nlm0=25, nrx=201)
      parameter (nkap0=4,n0=10)
      integer nglob,nbas,nspec,ib1,is,
     .  iclbsj,nr,j1,nch,lmxa,kmax,lh(nkap0),nkapi,ik,l1,l2,l,
     .  lmax,nlm,k,i,ilm,ipr,intopt
      double precision a,rmt,rsma,e,rsm,aa,r,rl,fith,fitj,fcth,fctj,wgt,
     .  c(0:k0,nlm0),rofi(nrx),rwgt(nrx),pkl(0:k0,0:n0),xi(0:n0),
     .  sumh(0:n0),difh(0:n0),errh(0:n0),errj(0:n0),cj(0:k0,0:n0),
     .  phi(0:n0),psi(0:n0),difj(0:n0),sumj(0:n0),rsmh(n0,nkap0),
     .  eh(n0,nkap0),xx
      character*8 spid
      character lch(0:4)
      data lch /'s','p','d','f','g'/

C ... setup
      stdo = lgunit(1)
      nbas  = nglob('nbas')
      nspec = nglob('nspec')
      call getpr(ipr)
      write (stdo,1)
    1 format(/' RMS errors in expanded head and tail functions.'/
     .        ' Only column ''tail'' matters for program accuracy.')

C     call defi(oips,nbas)
      allocate(ips(nbas))
      call sitepack(s_site,1,nbas,'spec',1,ips,xx)

C --- Loop over species ---
      do  is = 1, nspec
        ib1 = iclbsj(is,ips,-nbas,1)
        if (ib1 < 0) then
          if (ipr >= 20) write(stdo,
     .   '('' chkxpn (warning) no sites corresponding to species'',i3)')
     .      is
          cycle
        endif

        spid = s_spec(is)%name
        rmt = s_spec(is)%rmt
        rsma = s_spec(is)%rsma
        lmxa = s_spec(is)%lmxa
        kmax = s_spec(is)%kmxt
        a = 0.02d0
        nr = 101

        call word(spid,1,j1,nch)
        write (stdo,2) spid(1:nch),rmt,rsma,lmxa,kmax
    2   format(/' Species ',a,':   rmt=',f7.4,'   rsma=',
     .          f7.4/' augment to  lmax=',i2,'   kmax=',i3)

        if (nr > nrx) call rxi('chkxpn: need nrx ge',nr)
        call radmsh(rmt,a,nr,rofi)
        intopt = 10*nglob('lrquad')
        call radwgt(intopt,rmt,a,nr,rwgt)
        call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)

C  ...  Loop over orbitals for this species
        write (stdo,3)
    3   format(/' block    rsm        e      l',6x,'head',8x,'tail')
        do  ik = 1, nkapi
        l2 = -1
        do  l1  = 0, lh(ik)
          e = eh(l1+1,ik)
          rsm = rsmh(l1+1,ik)
          if (rsm /= 0 .and. l1 > l2) then
            l2 = l1-1
    4         continue
              l2 = l2+1
            if (l2 < lh(ik)) then
            if (rsmh(l2+2,ik) == rsm .and. eh(l2+2,ik) == e) goto 4
            endif

          lmax = l2
          nlm = (lmax+1)**2

          if (nlm > nlm0) call rxi('chkxpn: need nlm0 ge',nlm)
          if (kmax > k0)  call rxi('chkxpn: need k0 ge',kmax)

C     ... Compare expanded and exact head,tail functions
          do  l = 0, lmax
            difh(l) = 0d0
            sumh(l) = 0d0
            difj(l) = 0d0
            sumj(l) = 0d0
          enddo

          call hxpos(rsmh,rsma,eh,kmax,nlm,k0,c)
C         call hxps(rsmh,rsma,e,kmax,nlm,k0,c)
          call jxpos(rsma,e,kmax,lmax,k0,cj)
          aa = 1d0/rsm

          do  i = 1, nr
            r = rofi(i)
            call radpkl(r,rsma,kmax,lmax,k0,pkl)
            call hansmr(r,e,aa,xi,lmax)
            call bessl(e*r*r,lmax,phi,psi)
            do  l = 0, lmax
              if (l == 0) then
                rl = 1
              else
                rl = r**l
              endif
              ilm = l*l+1
              fith = 0d0
              fitj = 0d0
              do  k = 0, kmax
                fith = fith + c(k,ilm)*pkl(k,l)*rl
                fitj = fitj + cj(k,l)*pkl(k,l)*rl
              enddo
              fcth = xi(l)*rl
              fctj = phi(l)*rl
              wgt = rwgt(i)*r*r*rl
              difh(l) = difh(l) + wgt* (fith-fcth)**2
              sumh(l) = sumh(l) + wgt* fcth**2
              difj(l) = difj(l) + wgt* (fitj-fctj)**2
              sumj(l) = sumj(l) + wgt* fctj**2
            enddo
          enddo

          do  l = l1, l2
            errh(l) = dsqrt(difh(l)/sumh(l))
            errj(l) = dsqrt(difj(l)/sumj(l))
                write (stdo,5) l,rsm,e,lch(l),errh(l),errj(l)
    5           format(i4,f10.3,f9.3,2x,3x,a1,1p,6d12.2)


          enddo
          endif
        enddo
        enddo
      enddo
      deallocate(ips)
      end
