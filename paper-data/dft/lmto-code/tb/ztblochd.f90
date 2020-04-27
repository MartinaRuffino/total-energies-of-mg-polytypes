
      function islocal(coords,blcksz,procsz,proccr)
         implicit none
         integer, intent(inout) :: coords(2) ! global coordinates of the element
         integer, intent(in), dimension(2) :: blcksz, & ! block size
                                              procsz, & ! process array size
                                              proccr    ! process coordinates
         logical :: islocal

         integer, dimension(2) :: glblcr, & ! global block coordinates of the element
                                  lcblcr

         islocal = .false.
         glblcr = coords/blcksz
         if (any(proccr /= mod(glblcr,procsz))) return
         islocal = .not. islocal
         lcblcr = glblcr/procsz
         coords = lcblcr*blcksz + mod(coords, blcksz)

      end function islocal


      subroutine ztblochd(tbc,lgamma,bk,nl,nspc,nspu,ispu,nbas,plat,lmx,ipc, &
                     & indxsh,nsite,iax,npr,s,vso,hso,lso1,ldim,sk,lds,typesz)
!C- Bloch transform of real-space matrix stored in possibly larger array
!C ----------------------------------------------------------------------
!Ci Inputs
!Ci   lgamma: gamma point only, sk real
!Ci   bk:    k point / 2 pi
!Ci   nl,nbas,plat,lmx,ipc;
!Ci   nspc = 2 for coupled spins (empirical S-O)
!Ci   nspu = 2 if TB+U, ispu is then current spin
!Ci   indxsh: pointers for constructing permutated H and S from makidx
!Ci   nsite: total number of neighbors in all clusters;
!Ci   iax: neighbor lists; npr: see tbham;
!Ci   s:     real-space matrix that is to be summed;
!Ci   vso,hso: table of spin-orbit parameters and the hamiltonian;
!Ci   lso1:  if T add spin-orbit terms
!Ci   ldim: actual size of sk
!Ci   lds: leading dimension of sk
!Ci   typesz: The size of the type in real(8) words.
!Ci           1 for real(8) and 2 for complex(8)
!Co Outputs
!Co   sk: Bloch sum of real-space matrix S
!Cl Locals:
!Cl   indxH:  work array of length nbas
!Cr Remarks
!Cr   s(k;r1,l1,r2,l2) = sum_T s(r1,l1,T+r2,l2) * exp(-i k . T)
!Cr   where r1 and r2 are basis vectors and T = t2-t1 is the difference
!Cr   in primitive lattice translation vectors.
!C ----------------------------------------------------------------------
      use tbprl
      implicit none
!C Passed parameters
      type(tbc_t) :: tbc
      integer, intent(in) :: typesz, lds
      integer nl,nspc,nspu,ispu,nbas,nsite,ldim
      integer, parameter  :: niax=10
      integer lmx(*),ipc(nbas),indxsh(*),iax(0:niax-1,nsite),npr(0:1,nbas),indxH(nbas)
      real(8) bk(3),plat(3,3)
      real(8) s(nl**2,nl**2,nsite*nspc**2,nspu),vso(0:nl-1,1),hso(nl**2,nl**2,4,2),sk(typesz,lds,tbc%loclsz(2))
      logical lgamma,lso1

      interface
         function islocal(coords,blcksz,procsz,proccr)
            implicit none
            integer, intent(inout) :: coords(2) ! global coordinates of the element
            integer, intent(in), dimension(2) :: blcksz, & ! block size
                                                procsz, & ! process array size
                                                proccr    ! process coordinates
            logical :: islocal
         end function islocal
      end interface
!C Local parameters
      integer iat,isite,j,k,ix0,ixR,lm0,lmR,nlm0,nlmR,i,l2,lsp,icl0
      integer nlm,ll,iprint
      real(8) twopi,TdotK,cosT,sinT,v0
      logical lso
      integer, dimension(2) :: ij,blcksz,procsz,proccr
! A function
      nlm(i) = (1+lmx(i))**2

!       l2 = ldim**2
      lsp = ldim / nspc
      sk = 0.0_8
      twopi = 8*atan(1d0)

!C --- Form table of indices that mark offset of Rth block in S(k) ---
      indxH(1) = 0
      do  iat = 2, nbas
        indxH(iat) = indxH(iat-1) + nlm(ipc(iat-1))
      enddo

      blcksz = tbc % blcksz
      procsz = tbc % c2d % szs(0:1)
      proccr = tbc % c2d % crd(0:1)

!C --- For each RR' pair add contribution to Bloch sum ---
      if (nspc == 2) goto 1

      if (lgamma .and. tbc%sl) then
         do  isite = 1, nsite
            ix0 = indxH(iax(0,isite))
            ixR = indxH(iax(1,isite))
            nlm0 = nlm(ipc(iax(0,isite)))
            nlmR = nlm(ipc(iax(1,isite)))
            do  lmR = 1, nlmR
               do  lm0 = 1, nlm0
                  ij = [indxsh(ix0+lm0),indxsh(ixR+lmR)] - 1
                  if (all(ij < lsp)) then
                     if (islocal(ij,blcksz,procsz,proccr)) then
                        ij = ij + 1
                        sk(1,ij(1),ij(2)) = sk(1,ij(1),ij(2)) + s(lm0,lmR,isite,ispu)
                     end if
                  end if
               enddo
            enddo
         enddo
      else if (lgamma .and. .not. tbc%sl) then
         do  isite = 1, nsite
            ix0 = indxH(iax(0,isite))
            ixR = indxH(iax(1,isite))
            nlm0 = nlm(ipc(iax(0,isite)))
            nlmR = nlm(ipc(iax(1,isite)))
            do  lmR = 1, nlmR
               j = indxsh(ixR+lmR)
               do  lm0 = 1, nlm0
                  i = indxsh(ix0+lm0)
                  if (i <= lsp .and. j <= lsp) then
                     sk(1,i,j) = sk(1,i,j) + s(lm0,lmR,isite,ispu)
                  end if
               enddo
            enddo
         enddo
      else if (.not. lgamma .and. tbc % sl) then
         do  isite = 1, nsite
            TdotK = 0.0d0
            do  k = 1,3
               do  j = 1,3
                  TdotK = TdotK + bk(j)*plat(j,k)*iax(1+k,isite)
               enddo
            enddo
            TdotK = TdotK*twopi
            cosT = cos(TdotK)
            sinT = sin(TdotK)
            ix0 = indxH(iax(0,isite))
            ixR = indxH(iax(1,isite))
            nlm0 = nlm(ipc(iax(0,isite)))
            nlmR = nlm(ipc(iax(1,isite)))

            do  lmR = 1, nlmR
               do  lm0 = 1, nlm0
                  ij = [indxsh(ix0+lm0),indxsh(ixR+lmR)] - 1
                  if (all(ij < lsp)) then
                     if (islocal(ij,blcksz,procsz,proccr)) then
                        ij = ij + 1
                        sk(1,ij(1),ij(2)) = sk(1,ij(1),ij(2)) + s(lm0,lmR,isite,ispu)*cosT
                        sk(2,ij(1),ij(2)) = sk(2,ij(1),ij(2)) + s(lm0,lmR,isite,ispu)*sinT
                     end if
                  endif
               enddo
            enddo
         enddo
      else if (.not. lgamma .and. .not. tbc % sl) then
         do  isite = 1, nsite
            TdotK = 0.0d0
            do  k = 1,3
               do  j = 1,3
                  TdotK = TdotK + bk(j)*plat(j,k)*iax(1+k,isite)
               enddo
            enddo
            TdotK = TdotK*twopi
            cosT = cos(TdotK)
            sinT = sin(TdotK)
            ix0 = indxH(iax(0,isite))
            ixR = indxH(iax(1,isite))
            nlm0 = nlm(ipc(iax(0,isite)))
            nlmR = nlm(ipc(iax(1,isite)))

            do  lmR = 1, nlmR
               j = indxsh(ixR+lmR)
               do  lm0 = 1, nlm0
                  i = indxsh(ix0+lm0)
                  if (i <= lsp .and. j <= lsp) then
                     sk(1,i,j) = sk(1,i,j) + s(lm0,lmR,isite,ispu)*cosT
                     sk(2,i,j) = sk(2,i,j) + s(lm0,lmR,isite,ispu)*sinT
                  endif
               enddo
            enddo
         enddo
      end if
      if (lgamma .and. ((iprint() > 30 .and. ldim < 19) .or. iprint() > 60)) then
        print *, 'TBLOCH real hamiltonian :'
        do  j = 1, ldim
          write (*,'(18f10.6)') (sk(1,j,k),k=1,ldim)
        enddo
      endif
      return

!C --- NSPIN=2 and spin-orbit branch ---
    1 continue
      if (lgamma) call rx('tbloch: nspc not set up for gamma only')
      lso = lso1
      do  isite = 1, nsite
        TdotK = 0
        do  j = 1,3
          do  k = 1,3
            TdotK = TdotK + twopi*bk(j)*plat(j,k)*iax(1+k,isite)
          enddo
        enddo
        cosT = cos(TdotK)
        sinT = sin(TdotK)

!C --- Determine whether to add spin-orbit interaction terms ---
        if (lso1) then
          lso = nl >= 2 .and. isite == npr(1,iax(0,isite))+1
          icl0 = ipc(iax(0,isite))
        endif

        ix0 = indxH(iax(0,isite))
        ixR = indxH(iax(1,isite))
        nlm0 = nlm(ipc(iax(0,isite)))
        nlmR = nlm(ipc(iax(1,isite)))

        do  lm0 = 1, nlm0
          do  lmR = 1, nlmR
            j = indxsh(ixR+lmR)
            i = indxsh(ix0+lm0)
            if (i <= lsp .and. j <= lsp) then
              sk(1,i,j) = sk(1,i,j) + s(lm0,lmR,isite,ispu)*cosT
              sk(2,i,j) = sk(2,i,j) + s(lm0,lmR,isite,ispu)*sinT

              sk(1,i,j+lsp) = sk(1,i,j+lsp) + s(lm0,lmR,isite+2*nsite,ispu)*cosT
              sk(2,i,j+lsp) = sk(2,i,j+lsp) + s(lm0,lmR,isite+2*nsite,ispu)*sinT

              sk(1,i+lsp,j) = sk(1,i+lsp,j) + s(lm0,lmR,isite+1*nsite,ispu)*cosT
              sk(2,i+lsp,j) = sk(2,i+lsp,j) + s(lm0,lmR,isite+1*nsite,ispu)*sinT

              sk(1,i+lsp,j+lsp) = sk(1,i+lsp,j+lsp) + s(lm0,lmR,isite+3*nsite,ispu)*cosT
              sk(2,i+lsp,j+lsp) = sk(2,i+lsp,j+lsp) + s(lm0,lmR,isite+3*nsite,ispu)*sinT

              if (lso) then
!C --- Spin-orbit interaction terms ---
                v0 = vso(ll(lm0),icl0)
                sk(1,i,j) = sk(1,i,j) + v0*hso(lm0,lmR,1,1)
                sk(2,i,j) = sk(2,i,j) + v0*hso(lm0,lmR,1,2)

                sk(1,i,j+lsp) = sk(1,i,j+lsp) + v0*hso(lm0,lmR,2,1)
                sk(2,i,j+lsp) = sk(2,i,j+lsp) + v0*hso(lm0,lmR,2,2)

                sk(1,i+lsp,j) = sk(1,i+lsp,j) + v0*hso(lm0,lmR,3,1)
                sk(2,i+lsp,j) = sk(2,i+lsp,j) + v0*hso(lm0,lmR,3,2)

                sk(1,i+lsp,j+lsp) = sk(1,i+lsp,j+lsp) + v0*hso(lm0,lmR,4,1)
                sk(2,i+lsp,j+lsp) = sk(2,i+lsp,j+lsp) + v0*hso(lm0,lmR,4,2)
              endif
            endif
          enddo
        enddo
      enddo

      end subroutine ztblochd
