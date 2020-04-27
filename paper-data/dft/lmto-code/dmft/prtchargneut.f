      subroutine prtchargneut(io,no,iq,nk,ldhb,nit,V,Nnew,Nref,msg)
       implicit none
       integer, intent(in) :: io,no,iq,nk,nit,msg
       real(8), intent(in) :: Nnew,Nref,V
       logical, intent(in) :: ldhb

       select case(msg)

c        first call
         case(1)
           if(io==no.and.iq==nk) write(*,197) nk,no
  197 format(/' Computing V to ensure charge neutrality'/
     .        ' Hamiltonian constructed and diagonalised for',
     .        i5,' k-points and ',i5,' frequencies')

         case(2)
           if(ldhb) write(*,198) Nnew,Nref
  198 format(' Entering rfalsi. Actual valence charge=',f10.5,
     .       ' target = ',f10.5)

         case(3)
           write(*,199)Nnew,V,nit
  199 format(' New valence charge=',f10.5,' with V=',f10.5/
     .       ' rfalsi called ',i3,' times'/)

       end select
      end subroutine prtchargneut
