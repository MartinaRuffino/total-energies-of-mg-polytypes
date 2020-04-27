         subroutine addsigmatog(gij, sigma, strans, csize, sign, nev, nicix, ndsig,nev0,ndham)
          implicit none
          complex(8), intent(inout):: gij(ndham,ndham)
          complex(8), intent(in)   :: sigma(ndsig,nicix)
          complex(8), intent(in)  :: strans(ndsig,nicix,nev,nev)
          integer, intent(in)      :: csize(nicix), sign, nev, nicix, ndsig,nev0,ndham
           !----- locals
          complex(8) :: sigmaij
           integer    :: i, j, icix, it
           do i=1,nev                     ! over bands-1
              do j=1,nev                  ! over bands-2
                 sigmaij=dcmplx(0d0,0d0)
                 do icix=1,nicix
                    do it=1,csize(icix)
                       sigmaij = sigmaij + strans(it,icix,i,j)*sigma(it,icix)
                    enddo
                 enddo
                 sigmaij=dcmplx(0d0,0d0) ! --LS--! TO BE REMOVED
                 gij(i+nev0,j+nev0) = gij(i+nev0,j+nev0) + sign*sigmaij
              enddo
           enddo
         end subroutine addsigmatog
