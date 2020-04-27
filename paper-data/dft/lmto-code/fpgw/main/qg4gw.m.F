      program qg4gw
C- Generate required q and G vectors for GW calculations.
C ----------------------------------------------------------------------
Cr input files
Cr   LATTC: contains the the following
Cr    alat       : lattice constant in a.u.
Cr    QpGcut_psi : maxmum of |q+G| in a.u. in the expansion of the eigenfunction.
Cr    QpGcut_Cou : maxmum of |q+G| in a.u. in the expansion of the Coulomb matrix.
Cr    plat(1:3,1): 1st primitive translation vector in the unit of alat
Cr    plat(1:3,2): 2nd primitive translation vector
Cr    plat(1:3,3): 3rd primitive translation vector
Cr   SYMOPS : point group operations
Cr
Cr output files
Cr   QGpsi: q and G vector for the eigenfunction
Cr   QGcou: q and G vector for the Coulomb matrix
Cr   Q0P  : q points near \Gamma points
Cr   QIBZ : q points in the Irreducible BZ.
Cr   EPSwklm : offset Gamma method.
Cr   Qmtet: q vectors for divided-tetrahedron.
Cr --------------------------
Cr The total number of q points are nqnum=(n1q*n2q*n3q)*(1+ nq0i).
Cr QGpsi is written as follows. See mkqg2 in mkqg.F
Cr    open(ifiqg, file='QGpsi')
Cr    write(ifiqg ) nqnum,ngpmx,QpGcut_psi,nqbz,nqi,imx,nqibz
Cr    allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx) ) !inverse mapping table
Cr    ngvecprev=9999
Cr    ngveccrev=9999
Cr    do iq = 1, nqnum
Cr      q = qq(1:3,iq)
Cr      write (ifiqg) q, ngp, irr(iq) ! irr=1 for irreducible points
Cr      do ig = 1,ngp
Cr        nnn3 = ngvecp(1:3, ig) 
Cr        ngvecprev( nnn3(1), nnn3(2),nnn3(3)) = ig
Cr      enddo
Cr      write (ifiqg)  ngvecp,ngvecprev !ngvecprev is added on mar2012takao
Cr      do ig = 1,ngc
Cr        nnn3 = ngvecc(1:3, ig) 
Cr        ngveccrev( nnn3(1), nnn3(2),nnn3(3)) = ig
Cr      enddo
Cr    enddo  
Cr    close(ifiqg)
Cr --------------------------
Cr True q (in a.u. in Cartesian coordinates) is given by
Cr    True_q(1:3)     = 2*pi/alat * q(1:3)
Cr  True G is given by
Cr    qplusG(1:3,igp) = 2*pi/alat * matmul(qlat * ngvec(1:3,igp)) ,igp=1,ngp
C ----------------------------------------------------------------------
      use keyvalue,only: getkeyvalue
      implicit none
      integer(4) ::n1q,n2q,n3q,ifiqg,ifiqgc,ngrp,ifi,i,ig,iq0pin,idummy
      real(8) :: plat(3,3), alat,QpGcut_psi, QpGcut_Cou
      real(8) :: QpGx2,pi,unit
!     real(8) :: QpGx1 ! If QpGcut_psi is given by the code (LAPW case)
      real(8),allocatable :: symops(:,:,:)
C     character(len=150):: recrdxxx
C     character(len=10) :: keyw1='unit_2pioa',keyw2
      logical ::unit2=.false. !  readgwinput,
      integer(4)::nnn(3)
      integer(4):: q0pchoice
      logical:: lnq0iadd=.false. !keepppovl,
      integer(4):: iq0pinxxx ,ifile_handle
      integer:: gammacellctrl=0
      integer,parameter :: NULLI=-99999
      procedure(integer):: verbose

      pi= 4d0*atan(1d0)
      call cputid (0)
      write(6,*)' qg4gw: Generate Q0P->1; Readin Q0P->2; band mode->3; SW(chipm)->4'
      write(6,*)'        Generate Q0P->101(old offset Gamma)'
      write(6,*)'        Generate Q0P and Q0P for xyz ->201 '
      call getjobmode(iq0pin)
      if (iq0pin == NULLI) then
        read(5,*) iq0pin
      endif
      call headver('qg4gw',iq0pin)
      write(6,*) ' mode iq0pin = ',iq0pin
      if(iq0pin==-100.or.iq0pin==1.or.iq0pin==2.or.iq0pin==3.or.iq0pin==101) then
        iq0pinxxx=iq0pin
      elseif(iq0pin==10002) then
         iq0pinxxx=2
         gammacellctrl=1 !Gammacell skip mode
      elseif(iq0pin==20002) then
         iq0pinxxx=2
         gammacellctrl=2 !Gammacell only mode
      elseif(iq0pin==4) then
         iq0pinxxx=2
      elseif(iq0pin==201) then
        iq0pinxxx=1
        lnq0iadd=.true.
      else
        call rxi('qg4gw does not allow iq0pin=',iq0pin)
      endif

c iq0pin=-100: generate GWinput template ... no longer used
c      if(iq0pin==-100) then
c         call conv2gwinput()
c         call rx0( ' OK! qg4gw mode=-100 to generate GWinput')
c      endif
      idummy=q0pchoice()
      write(6,"(' q0pchoice() = ',i4)")  q0pchoice()

      ifi=ifile_handle()
      open (ifi, file='LATTC')
      read(ifi,*) alat
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
!     read(ifi,*) QpGx1 !If QpGcut_psi is always in Gwinput now
      close(ifi)

!! --- readin SYMOPS. point group operations. r'=matmul(symops(:,:),r) for any ig.
      ifi=ifile_handle()
      open (ifi, file='SYMOPS')
      read(ifi,*) ngrp
      write(6,*) ' SYMOPS ngrp=',ngrp
      allocate(symops(3,3,ngrp))
      do ig = 1,ngrp
        read(ifi,*)
        do i=1,3
          read(ifi,*) symops(i,1:3,ig)
        enddo
      enddo
      close(ifi)

C --- Informational printout, just repeat input ---
      write(6,"(' alat (a.u.)',f13.6 )") alat
      write(6,"(' plat a1   =',3f13.6)") plat(1:3,1)
      write(6,"(' plat a2   =',3f13.6)") plat(1:3,2)
      write(6,"(' plat a3   =',3f13.6)") plat(1:3,3)
      if (verbose() > 50) then
      write(6,*) ' --- point group operations ---'
      do ig = 1,ngrp
        print *, ' ig=',ig
        do i=1,3
          write(6,"(3f14.6)") symops(i,1:3,ig)
        enddo
      enddo
      endif

C --- Readin from GWinput n1q n2q n3q k-divisions for GW calculations
      call getkeyvalue("GWinput", "n1n2n3", nnn,3)
      n1q=nnn(1); n2q=nnn(2); n3q = nnn(3)
      call getkeyvalue("GWinput", "QpGcut_psi",QpGx2)
      call getkeyvalue("GWinput", "QpGcut_cou",QpGcut_Cou)
      call getkeyvalue("GWinput", "unit_2pioa",unit2)
      if(unit2) then
        unit = 2d0*pi/alat
        QpGx2     = QpGx2      *unit
        QpGcut_cou= QpGcut_cou *unit
      endif
      QpGcut_psi = QpGx2
      write(6,"(' k mesh from GWinput =',3i5)") nnn
      write(6,"(' |k+G| < QpG(psi) QpG(Cou)=',1p,2d13.6)") QpGcut_psi, QpGcut_Cou

C --- Open output files QGpsi QGcou
      ifiqg  = ifile_handle()
      open(ifiqg ,file='QGpsi',form='unformatted')
      ifiqgc = ifile_handle()
      open(ifiqgc,file='QGcou',form='unformatted')
      if(iq0pin==4) then
         QpGcut_psi=0d0
         QpGcut_Cou=0d0
      endif
!!
      call mkQG2(alat,plat,symops,ngrp,nnn,iq0pinxxx, 
     &     QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc, gammacellctrl,lnq0iadd)
      write(6,*) ' OK! End of qg4gw '
      if(iq0pin ==1) call rx0( ' OK! qg4gw mode=1 normal mode')
      if(iq0pin ==2) call rx0( ' OK! qg4gw mode=2 Readin Q0P mode')
      if(iq0pin ==10002) call rx0( ' OK! qg4gw mode=10002 Readin Q0P. GammaCell skipped.')
      if(iq0pin ==20002) call rx0( ' OK! qg4gw mode=20002 Readin Q0P. GammaCell Only.')
      if(iq0pin ==3) call rx0( ' OK! qg4gw mode=3 band-plot mode')
      if(iq0pin ==4) call rx0( ' OK! qg4gw mode=4 Readin Q0P mode. Set ngp=ngc=0')
      end
