      subroutine avsoef(lpzi,psi,pzi,nr,nsp,lmxa,
     .  v,dv,enum,ezum,z,rofi,rwgt,facso,sop,grv1,buf)
C- Calculates the SO weighted average of the electric field
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Cr Remarks
Cu Updates
Cu   24 Aug 15 (Scott McKecknie) SO weighted average of electric field
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lpzi(0:lmxa),nr,nsp,lmxa
      character*(*) buf
      double precision psi(nr,0:lmxa,nsp,3),
     .  pzi(nr, 0:lmxa,nsp),v(nr,nsp),dv(nr),enum(0:8,nsp),
     .  ezum(0:8,nsp),z,rofi(nr),rwgt(nr),facso,
     .  sop(0:lmxa,nsp,nsp,3),grv1(nr,3)

C ... Dynamically allocated local arrays
      real(8),allocatable :: efwi(:,:),
     .  sopn(:,:,:,:,:), efavg(:,:), sopzn(:,:,:,:,:)
C ... Local parameters
      integer i,j,l,stdo,nglob

C ... Setup
      stdo = nglob('stdo')

C ... Modify weights
      allocate(efwi(nr,3))
      do j=1,3
        do i=2,nr
          efwi(i,j)=grv1(i,j)*rwgt(i)
        enddo
      enddo

!      do i=1,nr
!        write(stdo,224) v(i,1), v(i,2)
!        write(stdo,224) grv1(i,1),grv1(i,2),grv1(i,3)
!      enddo

C ... Call soprm for each component of electric field
      allocate(sopn(0:lmxa,nsp,nsp,6,3), sopzn(0:lmxa,nsp,nsp,6,3))
      call pshpr(20) ! suppress printed output
      call soprm(1,lpzi,psi(1,0,1,1),psi(1,0,1,2),pzi,nr,nsp,lmxa,
     .  lmxa,v,dv,enum,ezum,z,rofi,efwi(1,1),facso,sopn(0,1,1,1,1),
     .  sopzn(0,1,1,1,1),buf)
      call soprm(1,lpzi,psi(1,0,1,1),psi(1,0,1,2),pzi,nr,nsp,lmxa,
     .  lmxa,v,dv,enum,ezum,z,rofi,efwi(1,2),facso,sopn(0,1,1,1,2),
     .  sopzn(0,1,1,1,2),buf)
      call soprm(1,lpzi,psi(1,0,1,1),psi(1,0,1,2),pzi,nr,nsp,lmxa,
     .  lmxa,v,dv,enum,ezum,z,rofi,efwi(1,3),facso,sopn(0,1,1,1,3),
     .  sopzn(0,1,1,1,3),buf)
      call poppr

C ... Divide by SO normalisation term (from main program call)
      allocate(efavg(lmxa,3))
      do j=1,3
        do l=1,lmxa
          efavg(l,j) = 0
          if (sop(l,1,1,1) /= 0) then
            efavg(l,j)=sopn(l,1,1,1,j)/sop(l,1,1,1)
          endif
        enddo
      enddo

C ... Print results
      write(stdo,222) ' avsoef:  Average electric field. '
      do l=1,lmxa
        write(stdo,223) '          ',
     .    l,efavg(l,1),efavg(l,2),efavg(l,3)
      enddo

  222 format(a)
  223 format(1x,a,i3,3(x,f20.8))
  224 format(3(1x,f20.8))

      end subroutine avsoef
