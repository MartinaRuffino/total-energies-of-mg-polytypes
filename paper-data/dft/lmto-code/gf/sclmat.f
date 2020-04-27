      module scale_matrix
      contains
      subroutine sclmat(mode,norb,nl,l,zm,dl,dr,zsh,zl,zr)
C- Scale matrix on right or left and shift
C-----------------------------------------------------------------------
Ci Inputs
Ci   mode  :job setting
Ci           0 do nothing
Ci           1 scale on left
Ci           2 scale on right
Ci           3 scale on left and right
Ci       add 4 shift by zsh
Ci   nl    :dimension for dl and dr and zsh
Ci    l    :index array norb->l
Ci   dl    :coefficients for left scaling. Coefficients are REAL, dimensioned nl
Ci   dr    :coefficients for right scaling. Coefficients are REAL, dimensioned nl
Ci   zsh   :complex shift, dimensioned nl
Cio Inputs/ Outputs
Cio  zm    :complex matrix to be scaled.  zm has norb dimensions (for lms)
Cu Updates
Cu   02 Apr 14 (Belashchenko) first created
C-----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,norb,nl,l(norb)
      real(8), intent(in), dimension(nl,2),optional :: dl,dr
      complex(8), intent(in), dimension(nl,2),optional :: zl,zr
      complex(8), intent(in), dimension(nl,2), optional :: zsh
      complex(8), intent(inout), dimension(norb,2,norb,2) :: zm
C Local parameters
      integer i1,i2,s1,s2,l1,l2
      logical lleft,lright,lshft,bittst
      complex(8), dimension(nl,2) :: zli,zri

      if (mode == 0) return
      call sanrg(.true.,mode,1,7,'sclmat','mode')
      lleft = bittst(mode,1)
      lright = bittst(mode,2)
      lshft = bittst(mode,4)
      if (lleft) then
        if (present(dl) .and. present(zl))
     .    call rx('sclmat: both dl and zl are present')
        if (present(dl)) then
          zli = dcmplx(dl)
        elseif (present(zl)) then
          zli = zl
        else
          call rx('sclmat: missing dl and zl')
        endif
      endif
      if (lright) then
        if (present(dr) .and. present(zr))
     .    call rx('sclmat: both dr and zr are present')
        if (present(dr)) then
          zri = dcmplx(dr)
        elseif (present(zr)) then
          zri = zr
        else
          call rx('sclmat: missing dr and zr')
        endif
      endif
      if (lshft .and. .not. present(zsh)) call rx('sclmat: missing zsh')

      do i1 = 1, norb
        l1 = l(i1) + 1
        do s1 = 1, 2
          do i2 = 1, norb
            l2 = l(i2) + 1
            do s2 = 1, 2
              if (lleft) zm(i1,s1,i2,s2) = zli(l1,s1) * zm(i1,s1,i2,s2)
              if (lright) zm(i1,s1,i2,s2) = zm(i1,s1,i2,s2) * zri(l2,s2)
            enddo
          enddo
          if (lshft) zm(i1,s1,i1,s1) = zm(i1,s1,i1,s1) + zsh(l1,s1)
        enddo
      enddo

      end subroutine sclmat
      end module scale_matrix

