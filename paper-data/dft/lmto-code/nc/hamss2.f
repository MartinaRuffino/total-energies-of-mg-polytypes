      subroutine hamss2(mode,ldim,theta,pph,sll,
     .  ccor,vmtz,elin,diags,wk,slld,h)
C- Spin-spiral, two-center hamiltonian
Ci mode: 1, copy s(k-q) to (11) block of h, 2: make 2-center, SS h.
Ci ldim: dimension of strux; theta, SS rotation angle
Co h: SS hamiltonian
      implicit none
      logical ccor
      integer ldim,mode
      double precision theta,sll(ldim,ldim,2),slld(ldim,ldim,2,2),
     .  h(ldim,2,ldim,2,2),
     .  pph(5,ldim,2),vmtz,elin,diags(ldim,0:2),wk(ldim)
      integer i,j,i1,j1
      double precision s,cc,ss,h1,h2,u1(2,2),u2(2,2),usu(2,2,2),x1,x2

C      call prmx('sll in hamss2',sll,ldim,ldim,ldim)
C      call prmx('slld1 in hamss2',slld,ldim,ldim,ldim)
C      call prmx('slld2 in hamss2',slld(1,1,1,2),ldim,ldim,ldim)

C --- Copy s(k-q) to (11) block of h ---
      if (mode == 1) then
        do i = 1, ldim
          do j = 1, ldim
          h(i,1,j,1,1) = sll(i,j,1)
          h(i,1,j,1,2) = sll(i,j,2)
          enddo
        enddo
        return
      endif

C --- Copy s(k+q) to (22) block of h ---
      do i = 1, ldim
        do j = 1, ldim
        h(i,2,j,2,1) = sll(i,j,1)
        h(i,2,j,2,2) = sll(i,j,2)
        enddo
      enddo

C --- Elements of rotation matrix for s ---
      cc = dcos(theta/2)**2
      ss = dsin(theta/2)**2
      s  = dsin(theta)/2
      u1(1,1) = cc
      u1(1,2) = -s
      u1(2,1) = -s
      u1(2,2) = ss
      u2(1,1) = ss
      u2(1,2) =  s
      u2(2,1) =  s
      u2(2,2) = cc

      x1 = vmtz-elin
      if (ccor) then
        do j = 1, ldim
        do i = 1, ldim
        x2 = (1 + (vmtz-elin)*(diags(i,1) + diags(j,1)))
        do i1 = 1, 2
        do j1 = 1, 2
          usu(i1,j1,1) = (h(i,1,j,1,1)*x2+slld(i,j,1,1)*x1)*u1(i1,j1) +
     .                   (h(i,2,j,2,1)*x2+slld(i,j,1,2)*x1)*u2(i1,j1)
          usu(i1,j1,2) = (h(i,1,j,1,2)*x2+slld(i,j,2,1)*x1)*u1(i1,j1) +
     .                   (h(i,2,j,2,2)*x2+slld(i,j,2,2)*x1)*u2(i1,j1)
        enddo
        enddo
        do i1 = 1, 2
          do j1 = 1, 2
            h(i,i1,j,j1,1) = pph(3,i,i1)*usu(i1,j1,1)*pph(3,j,j1)
            h(i,i1,j,j1,2) = pph(3,i,i1)*usu(i1,j1,2)*pph(3,j,j1)
          enddo
        enddo
        enddo
        enddo
      endif

      if (.not. ccor) then
C --- Make rotated (12,21) blocks of s ---
        do j = 1, ldim
          do i = 1, ldim
            h(i,2,j,1,1) = s*(h(i,2,j,2,1) - h(i,1,j,1,1))
            h(i,2,j,1,2) = s*(h(i,2,j,2,2) - h(i,1,j,1,2))
            h(i,1,j,2,1) = h(i,2,j,1,1)
            h(i,1,j,2,2) = h(i,2,j,1,2)
          enddo
        enddo

C --- Make rotated (11,22) blocks of s ---
        do j = 1, ldim
          do i = 1, ldim
            h1 = cc*h(i,2,j,2,1) + ss*h(i,1,j,1,1)
            h2 = cc*h(i,2,j,2,2) + ss*h(i,1,j,1,2)
            h(i,1,j,1,1) = cc*h(i,1,j,1,1) + ss*h(i,2,j,2,1)
            h(i,1,j,1,2) = cc*h(i,1,j,1,2) + ss*h(i,2,j,2,2)
            h(i,2,j,2,1) = h1
            h(i,2,j,2,2) = h2
          enddo
        enddo

C     call prmx('after 40',h,ldim*2,ldim*2,ldim*2)

C --- Add srdel S srdel into h ---
        do j = 1, ldim
        do i1 = 1, 2
        do j1 = 1, 2
        do i = 1, ldim
          h(i,i1,j,j1,1) = pph(3,i,i1)*h(i,i1,j,j1,1)*pph(3,j,j1)
          h(i,i1,j,j1,2) = pph(3,i,i1)*h(i,i1,j,j1,2)*pph(3,j,j1)
        enddo
        enddo
        enddo
        enddo

      endif

C --- H += C + (vmtz-eln)*<k|k>_constant ---
      call daxpy(ldim,1d0,pph(2,1,1),5,h,2*ldim+1)
      call daxpy(ldim,1d0,pph(2,1,2),5,h(1,2,1,2,1),2*ldim+1)
      if (ccor) then
        do  i = 1, ldim
          h(i,1,i,1,1) = h(i,1,i,1,1) +
     .      (vmtz-elin)*diags(i,0)*pph(3,i,1)**2
          h(i,2,i,2,1) = h(i,2,i,2,1) +
     .      (vmtz-elin)*diags(i,0)*pph(3,i,2)**2
      enddo
      endif

C     call prmx('after hamss',h,ldim*2,ldim*2,ldim*2)

      end

