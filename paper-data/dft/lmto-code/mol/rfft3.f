      subroutine rfft3(r,n1,n2,n3,iord,isig)
C- Three-dimensional real FFT; output Hermitian
Ci   r, array of points to transform
Ci   c must be dimensioned  (n1+2),n2,n3: size of array
Ci   isig:-1 for forward transform, 1 for reverse.
Ci   iord: 1, returns normal transform; 0: reverses order of mesh.
Ci   Note: iord=1 executes faster.
Co   r, FFT of input r
      implicit none
      integer id,n1,n2,n3,iord,isig
      double precision r(n1+2,n2,n3)
      integer ow1,ow2,ow3,oiwk,ierr,iopt,i,i1,i2,i3,n23
      double precision scale
      real w(1)
      common /w/ w

C --- If forward transform, zero out r(n1+1,*) r(n1+2,*) ---
      n23 = n2*n3
      if (isig == -1) then
        do  20  i = 1, n23
          r(n1+1,i,1) = 0d0
          r(n1+2,i,1) = 0d0
   20   continue
      endif
C --- Work arrays for transform ---
      id = n1+2
      iopt = 0
      if (n3 < 32) iopt = 1
      call defrr(ow1, 6*n1+14)
      call defrr(ow2, 4*n2*((1-iopt)+iopt*(id+1))+14)
      call defrr(ow3, 4*n3+14)
      call defrr(oiwk,max(n1,n2,n3))
C --- Do the FFT ---
      call r3fft(r,id,n1,n2,n3,w(ow1),w(ow2),w(ow3),iopt,0,
     .  iord,w(oiwk),ierr)
      if (ierr == 2) call rx('r3fft needs prime factors 2,3,5')
      if (ierr /= 0) call rx('r3fft called improperly')
      call r3fft(r,id,n1,n2,n3,w(ow1),w(ow2),w(ow3),iopt,isig,
     .  iord,w(oiwk),ierr)
      if (ierr /= 0) call rx('r3fft called improperly')
      call rlse(ow1)

C --- Renormalize forward transform ---
      if (isig > 0) return
      scale = 1/dble(n1*n2*n3)
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1+2
   10 r(i1,i2,i3) = scale*r(i1,i2,i3)

      end
