      subroutine ylmrtg(nlm,g,rmat)
C- Matrix to rotate cubic harmonics for a given rotation matrix
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nlm : (lmax+1)**2 for rmat
Ci   g   : 3x3 rotation matrix
Co Outputs:
Co   rmat: matrix that rotates Y_L(r) to Y_L(g r)
Cr Remarks:
Cr   Y_L(g r)   = sum_M rmat(L,M) Y_M(r)  and also
Cr   Y_L(g-1 r) = sum_M Y_M(r) rmat(M,L)
Cr   rmat is block diagonal in l.
Cr
Cr   ylmrtg makes rmat from (Y_M)^-1 Y_L at a collection of random points.
Cr   The block diagonal character of Y_L is exploited for inversion,
Cr   but not for the product (Y_M)^-1 Y_L, as it is assumed dgemm is fast enough.
Cr
Cr   Y^-1_M is made only once and retained, to speed up calculation of rmat,
Cr   Y^-1_M is remade if when nlm is increased
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer nlm
      double precision rmat(nlm,nlm),g(3,3)
C Local parameters:
      integer lmx,mmx,ndim,lmax,i,l,m,mmax,nm,offs,ierr,nlmi
      parameter ( lmx=12, mmx=2*lmx+1, ndim=(lmx+1)**2 )
      integer iwk(ndim)
      double precision p(3,mmx),rp(mmx,3),rr(mmx),
     .  yp(ndim,ndim),ylm(mmx,ndim),y(ndim,ndim)
      procedure(integer) :: ll
      integer nlmsav
      save nlmsav,y
      data nlmsav /0/
C ... A collection of random vectors, to make rmat = yrot(p) y^-1
      data p/ 0.020,-.025,-.118,
     .        0.419,-.538,0.513,
     .        0.245,-.717,-.600,
     .        -.056,0.224,-.309,
     .        -.034,-.180,0.207,
     .        -.351,-.614,0.950,
     .        -.782,-.134,-.308,
     .        0.568,0.716,-.457,
     .        -.528,-.927,-.562,
     .        -.856,-.443,0.267,
     .        -.111,0.794,0.598,
     .        -.985,-.144,-.617,
     .        0.678,0.400,-.617,
     .        0.730,-.207,-.101,
     .        0.540,-.137,-.773,
     .        -.758,-.992,-.561,
     .        0.321,-.363,-.988,
     .        0.132755,-.363770,0.562634,
     .        0.920691,-.208151,0.901722,
     .        0.020528,0.809339,0.161315,
     .        0.533549,-.520292,0.013359,
     .        -.161926,0.769089,-.874877,
     .        0.064178,-.216956,0.336349,
     .        0.212504,0.203445,0.481381,
     .        -.543106,-.288564,-.353783/

C     call tcn('ylmrtg')

C --- Initialization  ---
      lmax = ll(nlm)
      mmax = 2*lmax+1
      if (lmax > lmx) call rx('increase lmx in ylmrtg')
      call dpzero(rmat,nlm**2)

C     Special case g is a unit matrix
      if (sum(abs(g)) == 3 .and. g(1,1)==1 .and. g(2,2)==1 .and. g(3,3)==1) then
        forall (i = 1:nlm) rmat(i,i) = 1
        return
      endif

C --- Set up and invert matrix ylm(p) for l<=lmax ---
      if (nlm > nlmsav) then
        call dpzero(y,ndim**2)

C   ... Normalize p, make ylm(p)
        do  i = 1, mmax
          call dscal(3,1/dsqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2),p(1,i),1)
          rp(i,1) = p(1,i)
          rp(i,2) = p(2,i)
          rp(i,3) = p(3,i)
        enddo
        call ropyln(mmax,rp,rp(1,2),rp(1,3),lmax,mmx,ylm,rr)

C   ... Generate matrix y(p)
        do  i = 1, mmax
          do  l = 0, lmax
            nm = 2*l+1
            offs = l*l
            if (i <= nm) then
              do  m = 1, nm
                y(m+offs,i+offs) = ylm(i,m+offs)
              enddo
            endif
          enddo
        enddo

C   ... Invert matrix y(p)
        y(1,1) = 1/y(1,1)
        do  l = 1, lmax
          offs = l**2
          nlmi = 2*l+1

C         call prmx('y',y(offs+1,offs+1),ndim,nlmi,nlmi)
          call dgetrf(nlmi,nlmi,y(offs+1,offs+1),ndim,iwk,ierr)
          if (ierr /= 0) call rx('ylmrtg cannot invert YL')
          call dgetri(nlmi,y(offs+1,offs+1),ndim,iwk,yp,ndim,ierr)
C         call prmx('y',y(offs+1,offs+1),ndim,nlmi,nlmi)
        enddo

        nlmsav = nlm
C       print *, 'generated setup for ylmrtg'

      endif

C --- Set up matrix ylm(g*p) ---
C ... Make rp = g*p in with rp dimensioned (3,mmax)
C     call dgemm('N','N',3,mmax,3,1d0,g,3,p,3,0d0,rp,3)
C     call prmx('rp',rp,3,3,mmax)
C ... Make rp = g*p in with rp dimensioned (mmax,3)
      call dgemm('T','T',mmax,3,3,1d0,p,3,g,3,0d0,rp,mmx)
C     call prmx('rp',rp,mmx,mmax,3)
      call ropyln(mmax,rp,rp(1,2),rp(1,3),lmax,mmx,ylm,rr)

C ... Make matrix y(rp)
      do  i = 1, mmax
        do  l = 0, lmax
          nm = 2*l+1
          offs = l*l
          if (i <= nm) then
            do  m = 1, nm
              yp(m+offs,i+offs) = ylm(i,m+offs)
            enddo
          endif
        enddo
      enddo

C --- rmat = yrot * y^-1 ---
      rmat(1,1) = yp(1,1)*y(1,1)
      do  l = 1, lmax
        offs = l**2
        nlmi = 2*l+1

        call dgemm('N','N',nlmi,nlmi,nlmi,1d0,yp(offs+1,offs+1),ndim,
     .    y(offs+1,offs+1),ndim,0d0,rmat(offs+1,offs+1),nlm)

      enddo
C     call prmx('ylmrtg: rmat',rmat,nlm,nlm,nlm)

C     call tcx('ylmrtg')
      end
C#ifdefC TEST
C      subroutine fmain
CC- Test ylmrtg against original rotdlmm
C      implicit none
C      integer ifi,ng,ig,i,l,m1,m2,lmax,nlm,nl
C      parameter (lmax=12, nlm=(lmax+1)**2, nl=lmax+1)
C      integer, parameter :: ngmx=1
C      procedure(integer) :: fopng
C      real(8) rmat(nl*nl,nl*nl)
C      real(8) :: g(3,3,ngmx), dlmm(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ngmx)
C
CC     read group ops
C      ifi = fopng('SYMOPS',-1,1); rewind ifi
C      read(ifi,*) ng
C      ng = min(ng,ngmx)
C      do  ig = 1, ng
C        read(ifi,*)
C        do  i = 1, 3
C          read(ifi,*) g(i,1:3,ig)
C        enddo
C      enddo
C      close(ifi)
C
C      call rotdlmm(g,ng,lmax+1,dlmm) ! Using ORIG branch of GW code
C
C      do  ig = 1, ng
C        call ylmrtg(nl*nl,g(1,1,ig),rmat)  ! lm code
C        print *
C        print *, 'group op',ig
C        do  i = 1, 3
C          call info2(2,0,0,'%3;12,7D',g(i,:,ig),0)
C        enddo
C
CC       Compare
C        print *, 'rotation of real Ylm'
C        do  l = 0, nl-1
C        do  m2 = -l, l
C          write(*,"(25f12.7)") ( dlmm(m2,m1,l,ig), m1 = -l, l)
C        enddo
C        enddo
C        print *, 'subtract rmat'
C        do  l = 0, nl-1
C        do  m2 = -l, l
C        do  m1 = -l, l
C          dlmm(m2,m1,l,ig) = dlmm(m2,m1,l,ig) - rmat(l*l+l+1+m2,l*l+l+1+m1)
C        enddo
C        enddo
C        enddo
C        do  l = 0, nl-1
C        do  m2 = -l, l
C          write(*,"(25f12.7)") ( dlmm(m2,m1,l,ig), m1 = -l, l)
C        enddo
C        enddo
C        call info2(2,0,0,' max diff : %g',maxval(abs(dlmm(:,:,:,ig))),0)
C
C      enddo
C
C      end
C
C#define ORIG
C      subroutine rotdlmm(symops,ng,nl,dlmm)
CC- Generate rotation matrix D^l_{m,m'} for L-representation, given point groups
CC ----------------------------------------------------------------------
CCi Inputs
CCi   symops:point group operations
CCi   ng    :number of group operations
CCi   nl    :lmax+1
CCo Outputs
CCo   dlmm  :Rotation matrix D^l_{m,m'} that rotates Y_L(r) to Y_L(g r),
CCo         :where _L are real harmonics
CCo         :dimensioned dlmm(2*nl-1,2*nl-1,0:nl-1,ng,2)
CCu Updates
CCu   19 Aug 18 Generate D^l_{m,m'} by calling ylmrtg
CCu   04 May 16 Redesign of tolerance check
Cc-----------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer ng,nl
C      double precision symops(9,ng)
C      double precision dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
CC ... Local parameters
C      integer ig,l,m1,m2
C      real(8) rmat(nl*nl,nl*nl)
C
C#ifndefC ORIG
C      do  ig = 1, ng
C        call ylmrtg(nl*nl,symops(1,ig),rmat)  ! lm code
C        do  l = 0, nl-1
C        do  m2 = -l, l
C        do  m1 = -l, l
C          dlmm(m2,m1,l,ig) = rmat(l*l+l+1+m2,l*l+l+1+m1)
C        enddo
C        enddo
C        enddo
C      enddo
C#endifC
C
C
C#ifdefC ORIG
C      double complex   dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
C      double precision det,igann,osq2
C      double complex   msc(0:1,2,2), mcs(0:1,2,2),dum(2)
CC     parameter(Img=(0d0,1d0))
C      integer i,j,md,m,mx,ikap,is
C      double precision am(3,3),fac1,fac2,amv(9),dlength
C      double precision cbeta,beta,sbeta,calpha,salpha,alpha,cgamma,sgamma,co2,so2,gamma,add
C      real(8),parameter:: tolg = 1d-7
C      complex(8),parameter:: Img=(0d0,1d0)
C      equivalence(amv,am)
C
C      do 10 ig =1,ng
C        do 20 i=1,3
C        do 20 j=1,3
C          am(i,j) = symops(i+3*(j-1),ig)
C   20   continue
Cc calculate determinant(signature)
C        det= am(1,1)*am(2,2)*am(3,3)
C     &      -am(1,1)*am(3,2)*am(2,3)
C     &      -am(2,1)*am(1,2)*am(3,3)
C     &      +am(2,1)*am(3,2)*am(1,3)
C     &      +am(3,1)*am(1,2)*am(2,3)
C     &      -am(3,1)*am(2,2)*am(1,3)
C        if(abs(abs(det)-1d0) >= 1d-10) then
C          print *,' rotdlmm: det/=1 ig and det=',ig,det
C          call rx(' ')
C        endif
Cc seek Euler angle   print *,' goto cbeta',ig,det
C        cbeta = am(3,3)/det
Cc added region correction so as to go beyond domain error for functions, dsqrt and acos.
C        if(abs(cbeta-1d0) <= 1d-6) cbeta= 1d0
C        if(abs(cbeta+1d0) <= 1d-6) cbeta=-1d0
C        beta = dacos(cbeta)
C        sbeta= sin(beta)
Cc beta= 0~pi
C        if(sbeta <= 1.0d-6) then
C          calpha= 1d0
C          salpha= 0d0
C          alpha = 0d0
C          cgamma= am(2,2)/det
C          sgamma= am(2,1)/det
C        else
C          salpha =  am(2,3)/sbeta/det
C          calpha =  am(1,3)/sbeta/det
C          sgamma =  am(3,2)/sbeta/det
C          cgamma = -am(3,1)/sbeta/det
C        endif
C        co2 = dcos(beta/2)
C        so2 = dsin(beta/2)
Cc         print *,' calpha=',calpha
C        if(abs(calpha-1.0d0) <= 1.0d-6) calpha= 1.0d0
C        if(abs(calpha+1.0d0) <= 1.0d-6) calpha=-1.0d0
C        if(abs(cgamma-1.0d0) <= 1.0d-6) cgamma= 1.0d0
C        if(abs(cgamma+1.0d0) <= 1.0d-6) cgamma=-1.0d0
C        alpha=dacos(calpha)
C        if(salpha < 0d0) alpha=-alpha
C        gamma=dacos(cgamma)
C        if(sgamma < 0d0) gamma=-gamma
Cc         print *,'alpha beta gamma det=',alpha,beta,gamma,det
C        do 30 l =  0, nl-1
C        do 30 md= -l, l
C        do 30 m = -l, l
Cc  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wiley and Sons.  p.52 (4.13)
C          fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
C          fac2 = 0d0
C          do 40 ikap=0,2*l
C            if(l-md-ikap >= 0 .and. l+m-ikap >= 0 .and. ikap+md-m >= 0) then
C              add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap)
C     &            *igann(ikap+md-m)*igann(ikap) )
C              if(2*l+m-md-2*ikap /= 0) add=add*co2**(2*l+m-md-2*ikap)
C              if(md-m+2*ikap /= 0)     add=add*(-so2)**(md-m+2*ikap)
C              fac2 = fac2+add
C            endif
C   40     continue
Cc l-th rep. is odd or even according to (det)**l
C          dlmmc(md,m,l,ig) = fac1*fac2*det**l*
C     &        cdexp( -Img*(alpha*md+gamma*m) )
C   30   continue
C
C        am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
C        am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
C        am(1,3)= sin(beta)*cos(alpha)
C        am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
C        am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
C        am(2,3)= sin(beta)*sin(alpha)
C        am(3,1)=-sin(beta)*cos(gamma)
C        am(3,2)= sin(beta)*sin(gamma)
C        am(3,3)= cos(beta)
C
C        if (dlength(9,amv(:)*det-symops(1:9,ig),1) > tolg)
C     .    call rx('Rotation by symgrp does not match rotation by Euler angle')
C
CC        if(abs(am(1,1)*det-symops(1,ig)) > 1.0d-8.or.
CC     &     abs(am(2,1)*det-symops(2,ig)) > 1.0d-8.or.
CC     &     abs(am(3,1)*det-symops(3,ig)) > 1.0d-8.or.
CC     &     abs(am(1,2)*det-symops(4,ig)) > 1.0d-8.or.
CC     &     abs(am(2,2)*det-symops(5,ig)) > 1.0d-8.or.
CC     &     abs(am(3,2)*det-symops(6,ig)) > 1.0d-8.or.
CC     &     abs(am(1,3)*det-symops(7,ig)) > 1.0d-8.or.
CC     &     abs(am(2,3)*det-symops(8,ig)) > 1.0d-8.or.
CC     &     abs(am(3,3)*det-symops(9,ig)) > 1.0d-8) then
CC          call rx('Rotation by symgrp does not match rotation by Euler angle')
CC        endif
C
C
Cc        if(iprint() >= 140) then
CC        if(.false.) then
CC          print *;print *;print *,' **** group ops no. ig=', ig
CC          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
CC          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
CC          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
CC          print *,' by Eular angle '
CC          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
CC          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
CC          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
CC        endif
C 1731   format (' ',3f9.4)
C
C   10 continue
Cc conversion to cubic rep. Belows are from csconvs
Cc  msc mcs conversion matrix generation 2->m 1->-m for m>0
C      osq2 = 1d0/sqrt(2d0)
C      do m = 0,1
C        Msc(m,1,1)= osq2*(-1)**m
C        Msc(m,1,2)=-osq2*Img*(-1)**m
C        Msc(m,2,1)= osq2
C        Msc(m,2,2)= osq2*Img
C
C        Mcs(m,1,1)= osq2*(-1)**m
C        Mcs(m,1,2)= osq2
C        Mcs(m,2,1)= osq2*Img*(-1)**m
C        Mcs(m,2,2)=-osq2*Img
C      enddo
Cc
C      do 23 is=1,ng
C        if(.false.) then
Cc        if(iprint() >= 150) then
C          print *; print *,' **** group ops no. ig=', is
C          write(6,1731) symops(1,is),symops(4,is),symops(7,is)
C          write(6,1731) symops(2,is),symops(5,is),symops(8,is)
C          write(6,1731) symops(3,is),symops(6,is),symops(9,is)
C        endif
Cc convert to cubic rep.
C      do 23   l =0,nl-1
C        do 33 m2=-l,l
C        do 33 m1= 1,l
C          dum(1)= dlmmc(m2, m1,l,is)
C          dum(2)= dlmmc(m2,-m1,l,is)
C          mx    = mod(m1,2)
C          dlmmc(m2,  m1,l,is)=
C     &                       dum(1)*msc(mx,1,1)
C     &                      +dum(2)*msc(mx,2,1)
C          dlmmc(m2, -m1,l,is)=
C     &                       dum(1)*msc(mx,1,2)
C     &                      +dum(2)*msc(mx,2,2)
C   33   continue
C        do 43 m2=  1,l
C        do 43 m1= -l,l
C          dum(1)=dlmmc( m2, m1,l,is)
C          dum(2)=dlmmc(-m2, m1,l,is)
C          mx=mod(m2,2)
C          dlmmc( m2, m1,l,is)=
C     &                       mcs(mx,1,1)*dum(1)
C     &                      +mcs(mx,1,2)*dum(2)
C          dlmmc(-m2, m1,l,is)=
C     &                       mcs(mx,2,1)*dum(1)
C     &                      +mcs(mx,2,2)*dum(2)
C   43   continue
C        do 53 m2=-l,l
C        do 53 m1=-l,l
C          dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
C          if( abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12 ) call rx( ''//
C     &     ' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12')
C   53   continue
Cccccccccccccccccccccc
C        if(.false.) then
Cc        if(.true.) then
Cc        if(iprint() >= 41) then
C          print *; print *,'  points ops  ig, l=', is,l,' cubic   '
C          do m2=-l,l
C            write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
Cc    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
C          enddo
C        endif
Ccccccccccccccccccccccc
C   23 continue
C
C#endifC
C
C      end
C#ifdefC ORIG
C      double precision function igann(i)
C      igann  = 1d0
C      do ix =1,i
C        igann=igann*ix
C      enddo
C      end
C#endifC
C
C
C#endif
