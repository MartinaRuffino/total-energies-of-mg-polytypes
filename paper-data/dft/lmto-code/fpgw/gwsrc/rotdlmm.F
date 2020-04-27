      subroutine rotdlmm(symops,ng,nl,dlmm)
C- Generate rotation matrix D^l_{m,m'} for L-representation, given point groups
C ----------------------------------------------------------------------
Ci Inputs
Ci   symops:point group operations
Ci   ng    :number of group operations
Ci   nl    :lmax+1
Co Outputs
Co   dlmm  :Rotation matrix D^l_{m,m'} that rotates Y_L(r) to Y_L(rotp r),
Co         :where _L are real harmonics
Co         :Dimensioned dlmm(-lmax:lmax,-lmax:lmax,0:lmax,ng), lmax = nl-1
Cu Updates
Cu   19 Aug 18 Generate D^l_{m,m'} by calling ylmrtg
Cu   04 May 16 Redesign of tolerance check
c-----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,nl
      double precision symops(9,ng)
      double precision dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
C ... Local parameters
      integer ig,l,m1,m2
      real(8) rmat(nl*nl,nl*nl)

C#ifndef ORIG
      do  ig = 1, ng
        call ylmrtg(nl*nl,symops(1,ig),rmat)  ! From lm/subs
        do  l = 0, nl-1
        do  m2 = -l, l
        do  m1 = -l, l
          dlmm(m2,m1,l,ig) = rmat(l*l+l+1+m2,l*l+l+1+m1)
        enddo
        enddo
        enddo
      enddo
C#endif


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
C#endif

      end
C#ifdefC ORIG
C      double precision function igann(i)
C      igann  = 1d0
C      do ix =1,i
C        igann=igann*ix
C      enddo
C      end
C#endif

C#ifdefC TEST
C      subroutine fmain
CC- Test original rotdlmm against ylmrtg
CC  Symops read from SYMOPS file
C      implicit none
C      integer ifi,ng,ig,i,l,m1,m2
C      integer, parameter :: nl=4, ngmx=48
C      procedure(integer) :: fopng
C      real(8) rmat(nl*nl,nl*nl)
C      real(8) :: symops(3,3,ngmx), dlmm(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ngmx)
C!     real(8),allocatable :: dlmm(:,:,:,:)
C
C
CC     read group ops
C      ifi = fopng('SYMOPS',-1,1); rewind ifi
C      read(ifi,*) ng
C      ng = min(ng,ngmx)
C      do  ig = 1, ng
C        read(ifi,*)
C        do  i = 1, 3
C          read(ifi,*) symops(i,1:3,ig)
C        enddo
C      enddo
C      close(ifi)
C
C      call rotdlmm(symops,ng,nl,dlmm) ! old GW code
C      do  ig = 1, ng
C        call ylmrtg(nl*nl,symops(1,1,ig),rmat)  ! lm code
C        print *
C        print *, 'group op',ig
C        do  i = 1, 3
C          call info2(2,0,0,'%3;12,7D',symops(i,:,ig),0)
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
C
C      enddo
C
C      end
C#endif
