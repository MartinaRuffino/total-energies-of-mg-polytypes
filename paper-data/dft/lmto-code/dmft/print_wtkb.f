      subroutine print_wtkb(iq,iqfbz,isp,ndham,nsp,nkp,i2,wtkp,chempot,evl,ldawtkb,dmftwtkb)
C- Print out integration weights to accumulate DMFT charge density
C ----------------------------------------------------------------------
Ci Inputs
Ci   iq
Ci   iqfbz
Ci   isp   :current spin channel (1 or 2)
Ci   ndham :dimensioning parameter, largest hamiltonian dimension
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   i2    :max number of weights to print
Ci   dmftwtkb :dmft weights to print
Ci   ldawtkb  :lda/QSGW weights to print
Co Outputs
Co   Table of ldawtkb, dmftwtkb is printed
Cr Remarks
Cr   States need not be ordered
Cu Updates
Cu   09 Aug 16  (L. Sponza) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in)  :: iq,iqfbz,isp,ndham,nsp,i2,nkp
      real(8), intent(in)  :: chempot,evl(ndham,nsp)
      real(8), intent(in)  :: ldawtkb(ndham,nsp,nkp),dmftwtkb(ndham,nsp,nkp),wtkp(nkp)
C ... Local parameters
      integer :: i, nevl
      double precision xx(2),fac
      procedure(integer) :: iprint
      real(8) :: sumwt(2,2) = 0
      save sumwt

      if (iprint() < 40) return

      do  i = i2, 1, -1
        nevl = i
        if (abs(dmftwtkb(i,isp,iq)) > 1d-6) exit
      enddo

      xx(1) = sum(ldawtkb(1:i2,isp,iq))
      xx(2) = sum(ldawtkb(1:i2,isp,iq))
      fac = nsp/wtkp(iq)
      sumwt(isp,1:2) = sumwt(isp,1:2) + xx
      call info5(1,0,0,
     .  ' Weights for iq=%i (iqfbz=%i) spin=%i  sum wt %,7;7d (LDA)  %,7;7d (DMFT)',
     .  iq,iqfbz,isp,sum(ldawtkb(1:i2,isp,iq)),sum(dmftwtkb(1:i2,isp,iq)))
      call arrprt(' state  evl-mu    LDA       DMFT','%,4i%,6;10D%,6;10D%,6;10D',
     .  'Iddd',nevl,0,3,0,' |',
     .  xx,evl(1:nevl,isp)-chempot,ldawtkb(1:nevl,isp,iq)*fac,dmftwtkb(1:nevl,isp,iq)*fac,
     .  xx,xx,xx,xx)

      if (iq == nkp .and. isp == nsp) then
        call info5(1,0,0,' sum of LDA weights %,7;7d  sum of DMFT weights %,7;7d',
     .    sum(sumwt(1:2,1)),sum(sumwt(1:2,2)),3,4,5)
      endif

C      write(*,'("--- iq, iqfbz, ispin = ",3i6)') iq,iqfbz,isp
C      write(*,'("IBAND      LDA wtkb     DMFT wtkb")')
C      do i=1,i2
C        write(*,'(i5,2f14.8)') i, ldawtkb(i,isp,iq), dmftwtkb(i,isp,iq)
C      enddo

      end subroutine print_wtkb
