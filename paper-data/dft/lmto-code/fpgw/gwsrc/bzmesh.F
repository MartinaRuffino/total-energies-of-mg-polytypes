      subroutine bzmesh(plat,qb,n1,n2,n3,g,ng,ipq,qp,wgt,nq,nqmx,nadd,qbzshift)
C- Returns irreducible qp and their weights
C--------------------------------------------------------------------
Ci Inputs:
Ci  plat, direct lattice vectors; n1,n2,n3, number of divisions to be
Ci  made along each of the reciprocal lattice vectors to divide the
Ci  reciprocal lattice into microcells; g, symmetry group operations
Ci  output from GRPGEN; ng, number of point operations output from
Ci  nqmx, maximum number of k-points.
Co Outputs:
Co  ipq, wgt, nq (see remarks);
Co  qb, vectors of first microcell for input to BZINTS.
Cr Remarks:
Cr  The reciprocal lattice is divided into n1*n2*n3 microcells which
Cr  are parallelipipeds with 8 corners. The corners are nodes of the
Cr  k-space mesh in the whole reciprocal lattice unit cell. Some of
Cr  these will be symmetry-related leaving nq irreducible k-points.
Cr  These are returned in qp(3,j) j = 1,nq; for each corner defined
Cr  by the triple (i1,i2,i3), ipq(i1,i2,i3) points to the corresponding
Cr  k-vector in qp. wgt(j) then contains the sampling weight associated
Cr  with the j th irreducible k-point in qp : i.e. the number of times
Cr  it occurs in the whole reciprocal unit cell.
Cr
Cr t.kotani. Jan2004, modified 2016
Cr Original version was from  bz.f of ferdi/gw.
C---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,ng,nq,nqmx,nadd,ipq(n1+nadd,n2+nadd,n3+nadd)
      real(8) plat(3,3),qb(3,3),g(3,3,ng),qp(3,nqmx),wgt(nqmx),qbzshift(3)
C ... Local parameters
      real(8) rb(3,3),v(3),v1(3),qlat(3,3)
      real(8):: diff2(3),diff(3),ginv(3,3),xv(3),xvv(3)
      integer(4) :: verbose
      integer:: n123(3),n123a(3)
      integer:: j123(3),i1,i2,i3,ir,k,m,iq,ig,j1,j2,j3,iwgt,ipr
      real(8) :: w0,swgt,xq(3),xqd(3),xqp(3),xqm(3)
      real(8), parameter :: tolq=3e-6
      equivalence (j1,j123(1)),(j2,j123(2)),(j3,j123(3))

      call minv33tp(plat,qlat)
      n123(1) = n1; n123(2) = n2; n123(3) = n3; n123a = n123+nadd
      qb(:,1) = qlat(:,1)/n1; qb(:,2) = qlat(:,2)/n2; qb(:,3) = qlat(:,3)/n3
      rb(:,1) = plat(:,1)*n1; rb(:,2) = plat(:,2)*n2; rb(:,3) = plat(:,3)*n3
      ipr = verbose()
      if (ipr > 40) then
        write(*,"(' bzmesh : ',5x,'plat',31x,'qlat')")
        do k = 1, 3
          write(*,"(3f10.5,5x,3f10.5)") (plat(m,k),m=1,3),(qlat(m,k),m=1,3)
        enddo
      endif
      ipq = 0
      nq  = 0
      swgt = 0

      do  i3 = 1, n123a(3)
      do  i2 = 1, n123a(2)
      do  i1 = 1, n123a(1)

C   --- Add qp to list if not flagged as symmetry-related to a prior ---
        if (ipq(i1,i2,i3) .eq. 0) then

C         Possibly deform q, assign weight for this q
C         xq = (/I1-1,I2-1,I3-1/) / dble(n123)
C         deform xq -> xqd
          xq  = ([I1-1,I2-1,I3-1]) / dble(n123) + qbzshift   ! trial n123*qp, units of qlat
          xqd = xq              ! Possibly deformed qp (future)
          v = matmul(qlat,xqd)  ! Trial qp
          if (i1 <= n1 .and. i2 <= n2 .and. i3 <= n3) then
            xqp = xq-0.5d0/n123
            xqm = xq+0.5d0/n123
            w0 = 2*product(xqm-xqp)
          else
            w0 = 0              ! No weight for nadd points
          endif

C     ... Find each qp in the star of q as equivalent to this qp
          iwgt = 0
          do  ig = 1, ng
            v1 = matmul(g(:,:,ig),v)  ! rotated qp
            xv = matmul(v1,rb)        ! as multiple of qlat
            j123(:) = nint(xv(:)+1-qbzshift(:)) ! should be integer multiple
            xvv(:) = xv(:)+1-qbzshift(:) - j123 ! if so, this will vanish
            if (sum(abs(xvv(:)))>tolq) then
              if (nadd==0) then
                print *,'bzmesh: mapping i1i2i3 to: ',i1,i2,i3,xvv(:)
                call rx('bzmesh not mapped to regular mesh point')
              else
                cycle
              endif
            endif

C           Shift to first quadrant
            if (nadd==0) then    ! use mod(x+|x|*n1, n1) so that x+|x|*n1 >=0.
              j1 = mod(j1-1+n1*abs(j1-1),n1)+1
              j2 = mod(j2-1+n2*abs(j2-1),n2)+1
              j3 = mod(j3-1+n3*abs(j3-1),n3)+1
            else
              if(j1<1 .or. j1>n1+nadd) cycle
              if(j2<1 .or. j2>n2+nadd) cycle
              if(j3<1 .or. j3>n3+nadd) cycle
            endif

            if (ipq(j1,j2,j3) .eq. 0) then
              iwgt = iwgt+1
            endif
            ipq(j1,j2,j3) = nq+1
          enddo ! ig
          nq = nq+1
          qp(:,nq) = v
          wgt(nq) = iwgt*w0
          swgt = swgt + wgt(nq)
        endif ! set ipq(i1,i2,i3) and sym-related ipq(j1,j2,j3)

      enddo
      enddo
      enddo

      if (ipr > 50) then
        print *, 'ipq array ...'
        DO   I3 = 1, N3+nadd
          DO   I2 = 1, N2+nadd
            DO   I1 = 1, N1+nadd
              write(6,*)i1,i2,i3,ipq(i1,i2,i3)
            enddo
          enddo
        enddo
      endif

C ... Rotation check
C      call minv33(qlat,ginv)
C      do  i3 = 1, n3+nadd
C        do  i2 = 1, n2+nadd
C          do  i1 = 1, n1+nadd
C            v = matmul(qb,(/(i1-1),(i2-1),(i3-1)/))
C            do k  = 1,nq
C              do ir = 1,ng
C                diff = matmul(g(:,:,ir),qp(:,k)) - v
C                call rangedq(matmul(ginv,diff), diff2)
C                if(verbose()>60) then
C                  write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5,3d13.5)")
C     &             i1,i2,i3,v,qp(:,k),k,diff2 ! matmul(qp(:,nqx),plat)*8, matmul(v, plat)*8
C                endif
C                if(sum(abs(diff2))< tolq) then
C                  if(verbose()>50)
C     &             write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5)")
C     &             i1,i2,i3,v, qp(:,k),k ! matmul(qp(:,nqx),plat)*8, matmul(v, plat)*8
C                  goto 1022
C                endif
C              enddo
C            enddo
C            call rx( 'bzmesh: not find irotk; it may require accurate symmetry.')
C 1022       continue
C          enddo
C        enddo
C      enddo

C --- Printout ---
      if (dabs(swgt-2d0)> tolq) then
        print *,' swgt=',swgt,swgt/w0*2d0/n123
        call rx( 'BZMESH: error in weights')
      endif
      write(*,744) nq,(n1+nadd)*(n2+nadd)*(n3+nadd),n1+nadd,n2+nadd,n3+nadd
  744 FORMAT(/' BZMESH: ',I4,' irreducible k - points. from ',i5,' (',3i3,' )')
      write(6,663)
  663 format(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')
      do  51  iq = 1, nq
        iwgt = wgt(iq)/w0+ .1
        write(6,661) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,wgt(iq)
   51 continue
  661 format(i5,2x,3f10.4,i10,f16.6)
      end
C      SUBROUTINE BZMESH_old(icase,RBAS,QB,N1,N2,N3,G,NG,IPQ,QP,WGT,NQ,NQMX)
CC-Divides the reciprocal lattice into microcells
CCr new version. icase=2 for regular mesh not including gamma.
CCr t.kotani. Jan2004
CCr Original version was from  bz.f of ferdi/gw.
CC--------------------------------------------------------------------
CCi Inputs:
CCi  rbas, direct lattice vectors; n1,n2,n3, number of divisions to be
CCi  made along each of the reciprocal lattice vectors to divide the
CCi  reciprocal lattice into microcells; g, symmetry group operations
CCi  output from GRPGEN; ng, number of point operations output from
CCi  GRPGEN; nqmx, maximum number of k-points.
CCo Outputs:
CCo  ipq, wgt, nq (see remarks);
CCo  qb, vectors of first microcell for input to BZINTS.
CCr Remarks:
CCr  The reciprocal lattice is divided into n1*n2*n3 microcells which
CCr  are parallelipipeds with 8 corners. The corners are nodes of the
CCr  k-space mesh in the whole reciprocal lattice unit cell. Some of
CCr  these will be symmetry-related leaving nq irreducible k-points.
CCr  These are returned in qp(3,j) j = 1,nq; for each corner defined
CCr  by the triple (i1,i2,i3), ipq(i1,i2,i3) points to the corresponding
CCr  k-vector in qp. wgt(j) then contains the sampling weight associated
CCr  with the j th irreducible k-point in qp : i.e. the number of times
CCr  it occurs in the whole reciprocal unit cell.
CC---------------------------------------------------------------------
Cc icase==2 is for regular-mesh without gamma points.
C      IMPLICIT double precision (A-H,O-Z)
C      DIMENSION RB(3,3),QB(3,3),WGT(NQMX),V(3),V1(3),IPQ(N1,N2,N3),
C     .   RBAS(3,3),QBAS(3,3),QP(3,NQMX),g(3,3,1)
C      real(8):: diff2(3),diff(3),ginv(3,3),xv(3)
C      integer(4) :: verbose
C      hf=0d0
C      if(icase==2) hf=0.5d0
C!     call dinv33(rbas,1,qbas,vol)
C      call minv33(rbas,qbas)
C      WRITE(*,"(' BZMESH : ',5X,'RBAS',31X,'QBAS')")
C      DO K = 1, 3
C        WRITE(*,"(3F10.5,5X,3F10.5)") (RBAS(M,K),M=1,3),(QBAS(M,K),M=1,3)
C      enddo
C      DO  8  M = 1, 3
C        QB(M,1) = QBAS(M,1)/N1
C        QB(M,2) = QBAS(M,2)/N2
C        QB(M,3) = QBAS(M,3)/N3
C        RB(M,1) = RBAS(M,1)*N1
C        RB(M,2) = RBAS(M,2)*N2
C        RB(M,3) = RBAS(M,3)*N3
C    8 continue
C      IPQ = 0
C      W0  = 2D0/(N1*N2*N3)
C      NQ  = 0
C      SWGT = 0D0
C      ixx=0
C
Cc      DO IG = 1, NG
Cc        x1=.5
Cc        x2=1.5
Cc        x3=0
Cc        V = matmul(QBAS(:,:),(/x1,x2,x3/))
Cc        CALL GRPOP(V,V1,G,IG)
Cc        write(6,"(' ig qbas1 rot(qbas1)=', 3f8.3,'  ',3f8.3,'  ',3f8.3)")
Cc     &   v,v1,matmul(v1,rbas)
Cc      enddo
Cc      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
C      DO  20  I3 = 1, N3
C      DO  20  I2 = 1, N2
C      DO  20  I1 = 1, N1
C        IF (IPQ(I1,I2,I3) .EQ. 0) THEN
C          V = matmul(QB,(/(I1-1+hf),(I2-1+hf),(I3-1+hf)/))
C          IWGT = 0
C          DO  22  IG = 1, NG
C            CALL GRPOP(V,V1,G,IG)
C            XV = matmul(V1,rbas)
C            call rangedq_xxx(n1*xv(1)-hf,xv1)
C            call rangedq_xxx(n2*xv(2)-hf,xv2)
C            call rangedq_xxx(n3*xv(3)-hf,xv3)
C
Cc            write(6,"(' i1 i2 i3 ig=', 4i4,2x,3d13.6)") i1,i2,i3,ig
Cc     &      , n1*xv(1), n2*xv(2), n3*xv(3)
Cc            if(abs(xv1)+abs(xv2)+abs(xv3)>4d-6) then
Cc              write(6,"('             ', 18x,3d13.6)") xv1, xv2, xv3
Cc            else
Cc              write(6,"('             ', 18x,3d13.6,' go')") xv1, xv2, xv3
Cc            endif
C
C            if(abs(xv1)+abs(xv2)+abs(xv3)>3d-6) cycle
C            ! this is for the case not on the regular mesh points.
Cc
C            call rangedq_xxx(xv(1)-hf/n1,xv1)
C            call rangedq_xxx(xv(2)-hf/n2,xv2)
C            call rangedq_xxx(xv(3)-hf/n3,xv3)
Cc              write(6,"('             ', 18x,3d13.6,' xx')") xv1, xv2, xv3
C            j1 = xv1*n1+1 +1d-6
C            j2 = xv2*n2+1 +1d-6
C            j3 = xv3*n3+1 +1d-6
C            if(verbose()>=100)
C     &        write(6,"(' i1 i2 i3 j1 j2 j3=', 3i4,2x,3i4)") i1,i2,i3, j1,j2,j3
C            if(abs(j1-xv1*n1-1)/n1+abs(j2-xv2*n2-1)/n2+abs(j3-xv3*n3-1)/n3>3d-6) then
C              print *,'xxx1=',abs(j1-xv1*n1-1)
C              print *,'xxx2=',abs(j2-xv2*n2-1)
C              print *,'xxx3=',abs(j3-xv3*n3-1)
CCstop2rx 2013.08.09 kino              stop 'bzmesh: bug! abs(j1-xv1*n1-1)+... >4d-6'
C              call rx( 'bzmesh: bug! abs(j1-xv1*n1-1)+... >4d-6')
C            endif
C            IF (IPQ(J1,J2,J3) .EQ. 0) then
C              IWGT = IWGT+1
C            endif
C            IPQ(J1,J2,J3) = NQ+1
C            ixx  = ixx+1
C   22     continue
C          NQ = NQ+1
C          QP(:,NQ) = V
C          WGT(NQ) = IWGT*W0
C          SWGT = SWGT + WGT(NQ)
C        ENDIF
C   20 CONTINUE
C      if(verbose()>=100) print *,' Number of Qibz =',nq
C
Ccccccccccccccccccccccccccccc
Cc      stop 'test end of bzmesh 111 ----------------------------'
Ccccccccccccccccccccccccccccc
C
Cc-------------------------------------
Cc rotation check
C      call minv33(qbas,ginv)
C      DO  I3 = 1, N3
C        DO  I2 = 1, N2
C          DO  I1 = 1, N1
C            V = matmul(QB,(/(I1-1+hf),(I2-1+hf),(I3-1+hf)/))
Cc        write(6,"(' ===== qbz=',3f8.3)") matmul(ginv,v)
C            do k  = 1,nq
C              do ir = 1,ng
Cc        if(verbose()>104) print *,' grp=',ir !;      print *, grp(:,ir)
C                diff = matmul(g(:,:,ir),qp(:,k)) - v
C                call rangedq3(matmul(ginv,diff), diff2)
Cc          if(verbose()>104) write(6,"(' matmul(ginv,diff)=',3f8.3,' ',3f8.3)")
Cc     &      diff, matmul(ginv,diff)
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                if(verbose()>60) then
C                  write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5,3d13.5)")
C     &i1,i2,i3,v,qp(:,k),k,diff2 ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
C                endif
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C                if(sum(abs(diff2))< 1d-6) then
C                  if(verbose()>50)
C     &       write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5)")
C     &       i1,i2,i3,v, qp(:,k),k ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
Cc            irotk(k,ir)= kp
Cc            nstar(k)   = nstar(k) + 1
C                  goto 1022
C                endif
C
C              enddo
C            enddo
CCstop2rx 2013.08.09 kino            stop'bzmesh: not find irotk; it may require accurate symmetry.'
C            call rx( 'bzmesh: not find irotk; it may require accurate symmetry.')
C 1022       continue
C          enddo
C        enddo
C      enddo
Cc-----------------------------------------
C      if (dabs(swgt-2) > 1.d-9) then
C        print *,' swgt=',swgt
CCstop2rx 2013.08.09 kino        stop 'BZMESH: error in weights'
C        call rx( 'BZMESH: error in weights')
C      endif
C      write(*,744) nq,n1*n2*n3,n1,n2,n3
C  744 FORMAT(/' BZMESH: ',I4,' irreducible k - points. from ',i5,
C     .          ' (',3i3,' )')
C      write(6,663)
C  663 format(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')
C      do  51  iq = 1, nq
C        iwgt = wgt(iq)/w0+ .1
C        write(6,661) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,wgt(iq)
C   51 continue
C  661 format(i5,2x,3f10.4,i10,f16.6)
Ccccccccccccccccccccccccccccc
Cc      stop 'test end of bzmesh222 ----------------------------'
Ccccccccccccccccccccccccccccc
C      END
C
C
Cc----------------------------------------------------
C
C
C
C
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc This was old version from  bz.f of ferdi/gw.
Cc--------------------------------------------------------------------
C      SUBROUTINE BZMESH_oold(icase,RBAS,QB,N1,N2,N3,G,NG,IPQ,QP,WGT,NQ,NQMX)
CC-Divides the reciprocal lattice into microcells
CC--------------------------------------------------------------------
CCi Inputs:
CCi  rbas, direct lattice vectors; n1,n2,n3, number of divisions to be
CCi  made along each of the reciprocal lattice vectors to divide the
CCi  reciprocal lattice into microcells; g, symmetry group operations
CCi  output from GRPGEN; ng, number of point operations output from
CCi  GRPGEN; nqmx, maximum number of k-points.
CCo Outputs:
CCo  ipq, wgt, nq (see remarks);
CCo  qb, vectors of first microcell for input to BZINTS.
CCr Remarks:
CCr  The reciprocal lattice is divided into n1*n2*n3 microcells which
CCr  are parallelipipeds with 8 corners. The corners are nodes of the
CCr  k-space mesh in the whole reciprocal lattice unit cell. Some of
CCr  these will be symmetry-related leaving nq irreducible k-points.
CCr  These are returned in qp(3,j) j = 1,nq; for each corner defined
CCr  by the triple (i1,i2,i3), ipq(i1,i2,i3) points to the corresponding
CCr  k-vector in qp. wgt(j) then contains the sampling weight associated
CCr  with the j th irreducible k-point in qp : i.e. the number of times
CCr  it occurs in the whole reciprocal unit cell.
CCm Memory:
CCm  No significant internal requirement. No access to heap.
CC---------------------------------------------------------------------
Cc icase==2 is for regular-mesh without gamma points.
Cc
C      IMPLICIT double precision (A-H,O-Z)
C      DIMENSION RB(3,3),QB(3,3),WGT(NQMX),V(3),V1(3),IPQ(N1,N2,N3),
C     .   RBAS(3,3),QBAS(3,3),QP(3,NQMX),g(*)
Cc      print *,'bzmesh: iprint=',iprint(),n1,n2,n3
C      hf=0d0
C      if(icase==2) hf=0.5d0
C!     call dinv33(rbas,1,qbas,vol)
C      call minv33(rbas,qbas)
C      if (iprint() > 50) then
C        WRITE(*,"(' BZMESH : ',5X,'RBAS',31X,'QBAS')")
C        DO K = 1, 3
C          WRITE(*,"(3F10.5,5X,3F10.5)") (RBAS(M,K),M=1,3),(QBAS(M,K),M=1,3)
C        enddo
C      endif
C      DO  8  M = 1, 3
C        QB(M,1) = QBAS(M,1)/N1
C        QB(M,2) = QBAS(M,2)/N2
C        QB(M,3) = QBAS(M,3)/N3
C        RB(M,1) = RBAS(M,1)*N1
C        RB(M,2) = RBAS(M,2)*N2
C        RB(M,3) = RBAS(M,3)*N3
C    8 continue
C      IPQ = 0  !fortran 90
C      W0 = 2D0/(N1*N2*N3)
C      NQ = 0
C      SWGT = 0D0
C      NN1 = 6*N1-1
C      NN2 = 6*N2-1
C      NN3 = 6*N3-1
C      DO  20  I3 = 1, N3
C      DO  20  I2 = 1, N2
C      DO  20  I1 = 1, N1
C        IF (IPQ(I1,I2,I3) .EQ. 0) THEN
C          V(1) = (I1-1+hf)*QB(1,1) + (I2-1+hf)*QB(1,2) + (I3-1+hf)*QB(1,3)
C          V(2) = (I1-1+hf)*QB(2,1) + (I2-1+hf)*QB(2,2) + (I3-1+hf)*QB(2,3)
C          V(3) = (I1-1+hf)*QB(3,1) + (I2-1+hf)*QB(3,2) + (I3-1+hf)*QB(3,3)
C          IWGT = 0
C          DO  22  IG = 1, NG
C            CALL GRPOP(V,V1,G,IG)
C            X1 = V1(1)*RB(1,1) + V1(2)*RB(2,1) + V1(3)*RB(3,1) + .1
C            X2 = V1(1)*RB(1,2) + V1(2)*RB(2,2) + V1(3)*RB(3,2) + .1
C            X3 = V1(1)*RB(1,3) + V1(2)*RB(2,3) + V1(3)*RB(3,3) + .1
C            IF (X1 .GT. 0.D0) X1 = X1+1.D0
C            IF (X2 .GT. 0.D0) X2 = X2+1.D0
C            IF (X3 .GT. 0.D0) X3 = X3+1.D0
C            J1 = X1
C            J2 = X2
C            J3 = X3
C            J1 = MOD(J1+NN1,N1) + 1
C            J2 = MOD(J2+NN2,N2) + 1
C            J3 = MOD(J3+NN3,N3) + 1
C            IF (IPQ(J1,J2,J3) .EQ. 0) IWGT = IWGT+1
CC temp: try to find bug ...
C            if (j1 < 0) j1 = j1+n1
C            if (j2 < 0) j2 = j2+n2
C            if (j3 < 0) j3 = j3+n3
C            IPQ(J1,J2,J3) = NQ+1
C   22     continue
C          NQ = NQ+1
C          QP(1,NQ) = V(1)
C          QP(2,NQ) = V(2)
C          QP(3,NQ) = V(3)
C          WGT(NQ) = IWGT*W0
C          SWGT = SWGT + WGT(NQ)
C        ENDIF
C   20 CONTINUE
CCstop2rx 2013.08.09 kino      if (dabs(swgt-2) > 1.d-9) stop 'BZMESH: error in weights'
C      if (dabs(swgt-2) > 1.d-9) call rx( 'BZMESH: error in weights')
C      if (iprint() >= 30) then
Ccfa
Cc       do  30  j = 1, 2
Cc  30   write(lgunit(j),744) nq,n1*n2*n3,n1,n2,n3
C        write(*,744) nq,n1*n2*n3,n1,n2,n3
C  744   FORMAT(/' BZMESH: ',I4,' irreducible k - points. from ',i5,
C     .          ' (',3i3,' )')
C      endif
C      if (iprint() > 40) then
C        write(6,663)
C  663   format(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')
C        do  51  iq = 1, nq
C          iwgt = wgt(iq)/w0+ .1
C          write(6,661) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,wgt(iq)
C   51   continue
C  661   format(i5,2x,3f10.4,i10,f16.6)
C      endif
C      END
C      SUBROUTINE GRPOP(V,V1,G,I)
C      double precision G(3,3,*),V(3),V1(3)
C      integer i
C      V1(1) = G(1,1,I)*V(1) + G(1,2,I)*V(2) + G(1,3,I)*V(3)
C      V1(2) = G(2,1,I)*V(1) + G(2,2,I)*V(2) + G(2,3,I)*V(3)
C      V1(3) = G(3,1,I)*V(1) + G(3,2,I)*V(2) + G(3,3,I)*V(3)
C      RETURN
C      END
