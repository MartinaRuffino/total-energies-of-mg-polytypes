      subroutine fmain
      implicit none
      integer i,j,n,ndim,ier,nevmx,nevl,nev,itime,i1mach
      logical lov,lnv,lx
      parameter (ndim=4)
      complex*16 h(ndim,ndim),z(ndim,ndim),wk(32),c(ndim,ndim),
     .  s(ndim,ndim),zz(ndim,ndim),scale,
     .  h2(ndim+1,ndim+1),s2(ndim+1,ndim+1),wk2(ndim*ndim),
     .  z2(ndim+1,ndim+1)
      double precision e(ndim),e2(ndim),eevmx
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))

      eevmx=100
      nevmx=100
      lov = .true.
      lx  = .true.
      lnv = .false.
      itime = i1mach(2)
      call sa(h)
      print *, 'test zhev'
      print *, 'h:'
      print 333, h

C --- Check cplx2r ---
C      call zprm('(5f12.6)',h,ndim,ndim)
C      call cplx2r(ndim**2,0,h,z)
C      call zprm2('(5f12.6)',h,ndim,ndim)
C      call cplx2r(ndim**2,1,h,z)
C      call zprm('(5f12.6)',h,ndim,ndim)
C      stop

      call ss(s)
      print *, 's:'
      print 333, s

      call zhev(ndim,h,s,lov,lx,nevmx,eevmx,nev,wk,lnv,itime,e,z)

      print *, 'eigenvalues:'
      print 333, e
      call prtz('eigenvectors:',z,ndim,ndim,nev)
  333 format(4f14.5)

      call sa(h)
      n = ndim
      call zmpy(h,2*n,2,1,z,2*n,2,1,c,2*n,2,1,n,n,n)
      call ss(s)
      call zmpy(s,2*n,2,1,z,2*n,2,1,zz,2*n,2,1,n,n,n)

      print *, 'H Z  -  E O Z'
      do  10  i = 1, ndim
      do  10  j = 1, ndim
   10 print 333, c(i,j) - e(j)*zz(i,j)

      print *
      print *, 'test zhevx with lh=',ndim+1
      call sa(h)
      h2(1:ndim,1:ndim) = h
      call ss(s)
      s2(1:ndim,1:ndim) = s
      call zhevx(ndim,ndim+1,h2,s2,1,lx,nevmx,eevmx,nev,wk2,lnv,
     .  e2,ndim,z)

      print *, 'compare eigenvalues to those of zhev:'
      print 333, e
      print 333, e2

C      print *
C      print *, 'test zhevo with same overlap and lh=',ndim+1
C      call sa(h)
C      h2(1:ndim,1:ndim) = h
C      call ss(s)
C      s2(1:ndim,1:ndim) = s
C      call zhevo(ndim,ndim+1,h2,s2,nevmx,eevmx,1d-5,nevl,nev,e,e2,
C     .  ndim,z)
C      print *, 'eigenvalues of h,s from zhevo:'
C      print 333, e
C      print 333, e2
C      call prtz('eigenvectors:',z,ndim,ndim,nev)

      print *
      print *, 'test zhevo with small overlap and lh=',ndim+1
      call ss2(s)
      s2(1:ndim,1:ndim) = s
      print *, 's:'
      print 333, s
      call zhevx(ndim,ndim,s,s,0,lx,nevmx,eevmx,nev,wk2,lnv,e2,ndim,z)
      print *, 'eigenvalues of s from zhevx:'
      print 333, e2

      call sa(h)
      h2(1:ndim,1:ndim) = h
      call ss2(s)
      s2(1:ndim,1:ndim) = s

      call zhevx(ndim,ndim,h,s,1,lx,nevmx,eevmx,nev,wk2,lnv,
     .  e,ndim,z)
      print *, 'exact eigenvalues of h,s from zhevx:'
      print 333, e
      call prtz('eigenvectors:',z,ndim,ndim,nev)

C      call sa(h)
C      call ss2(s)
C      call zhev_takao(n,n,h,s,n,1d-2,nev,e,z)
C      print *, 'eigenvalues of h from takao:'
C      print 333, e(1:nev)
C      call prtz('eigenvectors:',z,ndim,ndim,nev)

      call sa(h)
      h2(1:ndim,1:ndim) = h
      call ss2(s)
      s2(1:ndim,1:ndim) = s
      call zhevo(ndim,ndim,h,s,nevmx,eevmx,1d-2,0,nevl,nev,e,e2,ndim,z)
C     Should get same result whether this call is commented out or not
      call zhevo(ndim,ndim+1,h2,s2,nevmx,eevmx,1d-2,0,nevl,nev,e,e2,
     .  ndim,z)

      print *, 'eigenvalues of h,s from zhevo:'
      print 333, e(1:nev)
      print 333, e2
      call prtz('eigenvectors:',z,ndim,ndim,nev)

      call sa(h)
      call ss2(s)
      call zgemm('C','N',nev,ndim,ndim,one,z,ndim,h,ndim,zer,s,ndim)
C     call zprm3('z+*h',2,s,ndim,nev,ndim)
      call zgemm('N','N',nev,nev,ndim,one,s,ndim,z,ndim,zer,h,ndim)
C     call zprm3('z+ * h * z',2,h,ndim,nev,nev)
      call prtz('result of z+ * h * z',h,ndim,nev,nev)

      end

      subroutine prtz(strn,z,ldz,ndim,nev)
      character strn*(*)
      integer ldz,ndim,nev
      double complex z(ldz,nev)
      print *, strn
C     print *, ldz,ndim,nev
      do  i = 1, ndim
        print 333, dble(z(i,1:nev))
        print 333, dimag(z(i,1:nev))
      enddo
C  333 format(4f14.5)
  333 format(4f20.10)

      end

      subroutine sa(a)
      complex*16 a(4,4)
      a(1,1) = 3+1
      a(2,2) = 3+1
      a(3,3) = 1+1
      a(4,4) = 1+1
      a(2,1) = 1
      a(3,1) = 0
      a(4,1) = (0d0,-2d0)
      a(3,2) = (0d0,+2d0)
      a(4,2) = 0
      a(4,3) = 1
      a(1,2) = 1
      a(1,3) = 0
      a(1,4) = (0d0,+2d0)
      a(2,3) = (0d0,-2d0)
      a(2,4) = 0
      a(3,4) = 1

      end
      subroutine ss(s)
      complex*16 s(4,4)

      call zinit(s,16)
      s(1,1) = 1
      s(1,3) = (.1d0,.2d0)
      s(3,1) = (.1d0,-.2d0)
      s(2,2) = 2
      s(3,3) = 3
      s(4,4) = 4
      end

      subroutine ss2(s)
      complex*16 s(4,4)

      call zinit(s,16)
      s(1,1) = 1d0 + .27d0 - .003d0 - .0003d0
      s(1,3) = (.1d0,.2d0)
      s(3,1) = (.1d0,-.2d0)
      s(1,4) = (1d0,2d0)
      s(4,1) = (1d0,-2d0)
      s(2,2) = 2
      s(3,3) = 3
      s(4,4) = 4
      end

      subroutine zprm(fmt,s,nr,nc)
      double precision s(2,nr,nc)
C#ifdef APOLLO_BUG
      character*(20) fmt
C#elseC
C      character*(*) fmt
C#endif
      print *, nr, nc
      do  10  i = 1, nr
   10 print fmt, (s(1,i,j), j=1,nc)
      do  20  i = 1, nr
   20 print fmt, (s(2,i,j), j=1,nc)
      end
      subroutine zprm2(fmt,s,nr,nc)
      double precision s(nr,nc,2)
C#ifdef APOLLO_BUG
      character*(20) fmt
C#elseC
C      character*(*) fmt
C#endif
      print *, nr, nc
      do  10  i = 1, nr
   10 print fmt, (s(i,j,1), j=1,nc)
      do  20  i = 1, nr
   20 print fmt, (s(i,j,2), j=1,nc)
      end

      subroutine zprm3(strn,icast,s,ns,nr,nc)
C     implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(10) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      fmt = '(9f15.10)'
      fmt = '(9f18.11)'
C     fmt = '(5f20.15)'
      outs = ' '
      if (icast == 1)  outs = ' real'
      if (icast == 11) outs = ' symm'
      if (icast == 2)  outs = ' complex'
      if (icast == 12) outs = ' herm'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
C#else
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)
C#endif
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(1,i,j), j=1,nc)
      if (mod(icast,10) > 1) then
      write(ifi,*)
      do  20  i = 1, nr
   20 write(ifi,fmt) (s(2,i,j), j=1,nc)
      endif
      call fclose(ifi)
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#else
      outs = ' zprm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in zprm'
C#else
      if (outs == 'q') call rx0('quit in zprm')
C#endif
      end
C
C      subroutine zhev_takao(nlmto,n,h,s,nmx, epsovl,
C     &   nev,e,z)
CC- Eigenvalues and/or some eigenvectors of a Hermitian matrix (throw away poor linear dependency)
CC ----------------------------------------------------------------
CCi Inputs:
CCi   nlmto: not used (dimension of MTO)
CCi   n:    dimension of h
CCi   h,n:  hermitian matrix, dimensioned h(n,n)
CCi   s:    hermitian overlap matrix, (used only if lov is true)
CCi   nmx:  maximum number of eigenvectors to be found
CCi   emx:  eigenvalue limit for eigenvectors to be found
CCi   wk:   work array of length at least 11n
CCi   epsovl: threshold to throw away linear dependence of the overlap matrix s.
CCo Outputs:
CCo   e:    eigenvalues
CCo   nev:  number of eigenvectors found
CCo   z:    eigenvectors (1..nev)  (declared as z(n,*)
CCo   s:    has been decomposed into and LL+ decomposition.
CCo         You can call zhev2 to scale a vector by L
CCr Remarks:
CCr   z must be at least of dimension z(n,n), even though nev<n.
CCr   h and s are destroyed on exit.
CCr   Aborts on exit
CCp Procedures used:
CCp   (lapack)  zhegv
CCr
CCr This routine throw away basis set which are poorly linear dependent on others.
CCr We first find the Hilbert space which is linearly independent.
CCr Then solve <zz|H|zz> of eigenvalue problem, where |zz> spans the Hilbert space.
CCr It eigenfucntion is moved back to the original Hamiltonian's eigenfunction.
CCr So this can be used just as a replacement of zhev.
CCr
CCu Updates
CCu   8  Jul 08 kotani renewed it. Linear-dependency removal
CCu   This is originally from Mark'z zhev.
CC ----------------------------------------------------------------
C      implicit none
CC Passed parameters
C      integer n,nev,nmx,ltime,ngv
C      double precision e(n),wk(1),emx
C      complex(8):: h(n,n),s(n,n),z(n,n)
C      complex(8),allocatable:: omat(:,:),imat(:,:),wk11(:),
C     & zz(:,:),hh(:,:),omat2(:,:),pmat(:,:),omat2i(:,:),oo(:,:)
C      integer(4):: i,j,nlmto,ni,nm,ik,ik2
C      real(8):: omax_mto, omax_g,ex ,  ddd2,ddd,epsovl
C      real(8),allocatable:: scale(:),rwork(:),eo(:)
CC Local parameters
C      integer ier,lwork,ix
C      character jobz
Cc
C
CC      call zprm3('s',2,s,n,n,n)
CC      call zprm3('h',2,h,n,n,n)
C
C      call tcn('zhev_takao')
Cccccccccccccccccccccccccccccc
Cc      do i=1,n
Cc      do j=1,n
Cc        if(i<5) print *,' 444 sss=',j,i,abs(s(j,i))
Cc      enddo
Cc      enddo
Cccccccccccccccccccccccccccc
C
C
C! if mode 2. For mode 1, comment 'goto 100'
Cc      goto 100
C
C
CC --- mode 1. Remove poor linear dependency. All basis are treated on the same footing.
CC ... rescale omat
Cc      allocate(omat2(n,n))
Cc      omat2=s
C      allocate(scale(n))
C      do i=1,n
C        scale(i)= 1d0/sqrt(s(i,i))
C      enddo
C
C      do i=1,n
C      do j=1,n
Cc        omat(i,j)=scale(i)*scale(j)*omat(i,j) !this is based on  1d0/scale*|basis>
C        s(i,j)= scale(i)*scale(j)*s(i,j)
C        h(i,j)= scale(i)*scale(j)*h(i,j) !this is to isolate the degenerated space.
C      enddo
C      enddo
C      allocate(omat(n,n))
C      omat = s !reserved
Ccccccccccccccccccccccccccccc
Cc      do i=1,n
Cc         ddd=0d0
Cc         do j=1,n
Ccc         if(i<5) print *,' sss=',j,i,abs(s(j,i))
Cc            if(j==i) cycle
Cc            if(abs(omat(j,i) )>ddd) then
Cc              ddd2=ddd
Cc              ik2=ik
Cc              ddd=abs(omat(j,i))
Cc              ik=j
Cc           endif
Cc          enddo
Cc         write(6,"(a,2i4,f9.3,i4,f9.3)")'mmm max ', i, ik,ddd,ik2,ddd2
Cc      enddo
Cccccccccccccccccccccccccccc
C
CC ... eigenvalue of s (omat)
C      jobz = 'V'
C      lwork = n*n
C      allocate(wk11(lwork),rwork(max(1,3*n-2)),eo(n))
C      call zheev(jobz,'U',n,omat,n,eo,wk11,lwork,rwork,ier)
C      deallocate(wk11,rwork)
C      write(6,*)'zhev_takao: ovlmat='
C      write(6,"(5(i5,d10.2))") (i,eo(i),i=1,n)
C
C! We solve  (z* h z -  epsilon z* o z) a =0
C      do ix= 1,n
C        if(eo(ix)>epsovl) then
C           ni = ix ! take i = ni...n
C           exit
C        endif
C      enddo
C
CC ... We solve solution in the Hilbert space spanned by zz (dimension = nm)
C      nm = n-ni+1            ! this is the dimension.
C      allocate(zz(n,nm))     ! zz is the projection matrix
C      do ix=ni,n
C        zz(:,ix-ni+1) = omat(:,ix)/sqrt(eo(ix))
C      enddo
C
CC ... Hamiltonian  <zz|H|zz>
C      allocate(hh(nm,nm))
C      hh = matmul(dconjg(transpose(zz)),matmul(h,zz))
CC     call zprm3('hh',2,hh,nm,nm,nm)
C
C      e=99.d0 !initialization
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      do ix=1,nm
C        write(6,"('hh dia',i5,2f13.4)") ix,hh(ix,ix)
C      enddo
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
CC ... diagonalize
C      if (nmx <= 0) then
C        jobz = 'N'
C        lwork = 4*nm
C        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
C        call zheev(jobz,'U',nm,hh,nm,e,wk11,lwork,rwork,ier)
C        nev = 0
C      else
C        jobz = 'V'
C        lwork = nm*nmx
C        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
C        call zheev(jobz,'U',nm,hh,nm,e,wk11,lwork,rwork,ier)
C        z=1d99
C        do i=1,min(nmx,nm)
C        do j=1,n
C           z(j,i) = scale(j) *sum(zz(j,:)*hh(:,i))  !this is eigenfunction for original problem.
C        enddo
C        enddo
C        nev = min(nmx,nm)
C      endif
C      deallocate(wk11,rwork)
C      deallocate(omat,scale,zz,eo)
C
Cc      do i=1,nev
Cc      do j=1,nev
Cc        if(abs(sum( dconjg(z(:,i))*matmul(omat2,z(:,j)) ))>1d-10) then
Cc        write(6,"(' matooo ',2i4,2d13.6)")  i,j,
Cc     &   sum( dconjg(z(:,i))*matmul(omat2,z(:,j)) )
Cc        endif
Cc      enddo
Cc      enddo
C
Ccccccccccccccccccccccccc
Cc      do i=1,nm
Cc        write(6,"(a,i5,12d13.5)") 'eee=',i,e(i)
Cc      enddo
Cc      stop ' --- vvv:zhev --- '
Cccccccccccccccccccccccccc
C      call rxx(ier /= 0,'zhev: zhegv cannot find all evals 111')
C      goto 999 !exit
C
C
C
C
CC --- mode 2. rescpect MTO --- a branch. not good.
C! I tested it; it seems not good. ---------------------------------
C 100  continue
CC ... PW part ! ... project out mto part
C      allocate(omat2(nlmto,nlmto),omat2i(nlmto,nlmto))
C      omat2  = s(1:nlmto,1:nlmto) !omat2=<chi_i|chi_j>
C      omat2i = omat2
C      call matcinv(nlmto,omat2i)  !omat2i is inverse of <chi_i|chi_j>
C      ngv = n - nlmto
C      deallocate(omat2)
C      allocate( omat(ngv,ngv) )
Cc      omat=0d0
Cc      do i=1,ngv
Cc         omat(i,i) = s(nlmto+i,nlmto+i)
Cc        print *,'qqq',i,omat(i,i)
Cc      enddo
C! <PW'|PW'> matrix. PW' = |PW> - |chi_i> O^(-1)_ij <chi_j|PW>
C       omat = s(nlmto+1:n,nlmto+1:n)  ! n=nlmto+ngv
C     .      -  matmul( dconjg(transpose( s(1:nlmto,nlmto+1:n))) ,
C     .                   matmul( omat2i, s(1:nlmto,nlmto+1:n))   )
CC ... eigenvalue of <PW'|PW'>
C      jobz = 'V'
C      lwork = ngv*ngv
C      allocate(wk11(lwork),rwork(max(1,3*ngv-2)),eo(ngv))
C      call zheev(jobz,'U',ngv,omat,ngv,eo,wk11,lwork,rwork,ier)
C      deallocate(wk11,rwork)
C      write(6,*)'zhev_takao: ovlmat='
C      write(6,"(5(i5,d10.2))") (i,eo(i),i=1,ngv)
C! skip low  eigenvalue
C      do ix= 1,ngv
C        if(eo(ix)>epsovl) then
C           ni = ix ! we use i = ni...n. Skip 1...ni-1
C           exit
C        endif
C      enddo
C      nm = ngv-ni+1  +nlmto      ! this is the dimension.
CC ... We solve solution in the Hilbert space spanned by zz (dimension = nm)
C      print *,' nlmto ngv nm=',nlmto,ngv,nm
C      allocate(zz(n,nm))
C      zz=0d0
C      do ix=1,nlmto
C         zz(ix,ix) = 1d0
C      enddo
C      do ix= 1, nm-nlmto
Cc        zz(ix+nlmto,ix+nlmto) = omat(ix,ix) /sqrt(eo(ix))
Cc        print *,' ppp omat',ix,zz(ix,ix)
C        zz(nlmto+1:n,ix+nlmto) = omat(:,ix) /sqrt(eo(ix))
C        zz(1:nlmto,  ix+nlmto) =
C     &    - matmul(omat2i, matmul(s(1:nlmto,nlmto+1:n),omat(:,ix)))
C     &    /sqrt(eo(ix))
C      enddo
CC ... Hamiltonian  <zz|H|zz>  <zz|Ozz>
C      allocate(hh(nm,nm),oo(nm,nm))
C      hh = matmul(dconjg(transpose(zz)),matmul(h,zz))
C      oo = matmul(dconjg(transpose(zz)),matmul(s,zz))
C      e = 99.d0 !initialization
CC ... diagonalize
C      if (nmx <= 0) then
C        jobz = 'N'
C        lwork = 4*nm
C        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
C        call zhegv(1,jobz,'U',nm,hh,nm,oo,nm,e,wk11,lwork,rwork,ier)
C        nev = 0
C      else
C        jobz = 'V'
C        lwork = nm*nmx
C        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
C        call zhegv(1,jobz,'U',nm,hh,nm,oo,nm,e,wk11,lwork,rwork,ier)
C        z=1d99
C        do i=1,min(nmx,nm)
C        do j=1,n
C           z(j,i) = sum(zz(j,:)*hh(:,i))  !this is eigenfunction for original problem.
C        enddo
C        enddo
C        nev = min(nmx,nm)
C      endif
Cccccccccccccccccccccccccc
Cc      do i=1,nm
Cc        write(6,"(a,i5,12d13.5)") 'zhev mode2: eee=',i,e(i)
Cc      enddo
Cc      stop ' --- vvv:zhev --- '
Cccccccccccccccccccccccccc
C      deallocate(wk11,rwork)
C      deallocate(omat,zz,eo,oo,omat2i)
C      call rxx(ier /= 0,'zhev: zhegv cannot find all evals 222')
Cc      stop 'xxxxxxxxxxxxxxxxxxx'
C
C! exit
C 999  continue
C      call tcx('zhev_takao')
C      end
C      subroutine matcinv(n,a)
CC --- a inverse is returned.
C      implicit none
C      integer(4) :: n, info, ipiv(n)
C      complex(8):: a(n,n)
C      complex(8),allocatable:: work(:)
C      call zgetrf(n,n,a,n,ipiv,info)
C      if(info/=0) then
C        print *,' matcinv: zegtrf info=',info
C        stop    ' matcinv: zegtrf '
C      endif
C      allocate(work(n*n))
C      call zgetri(n,a,n,ipiv,work,n*n,info)
C      deallocate(work)
C      if(info/=0) then
C        print *,'matcinv: zegtri info=',info
C        stop    'matcinv: zegtri '
C      endif
C      end
