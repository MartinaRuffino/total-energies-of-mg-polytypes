      subroutine gfbz2r(n1,n2,n3,nsp,idim,jdim,ldiag,gij,gji)
C- Subblock of real-space Green's function by inverse FT of G(k)
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1,n2,n3  QP mesh
Ci   idim,jdim dimensions of G(i,j) subblock
Ci   ldiag     T, gji and gij are identical
Cio Inputs/Outputs
Cio  gij,gji  Green's function connecting orbitals ij and ji
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      logical ldiag
      integer idim,jdim,n1,n2,n3,nsp
      double complex gij(n1,n2,n3,idim,jdim,nsp)
      double complex gji(n1,n2,n3,jdim,idim,nsp)
C Local variables
      integer id,jd,isp,iset
C --- Loop over orbital pairs ---
      iset= 0
      call fftz3s(n1,n2,n3,n1,n2,n3,iset)
      do  10  isp = 1, nsp
      do  10  id = 1, idim
      do  10  jd = 1, jdim

        print *, id,jd
        call zprm3(' gij %i',id+100+jd,gij(1,1,1,id,jd,isp),n1,n2,n3)
c        call zprm3(' gji %i',id+100+jd,gji(1,1,1,id,jd,isp),n1,n2,n3)

        call fftz3(gij(1,1,1,id,jd,isp),n1,n2,n3,n1,n2,n3,1,iset,-1)
        if (.not. ldiag)
     .  call fftz3(gji(1,1,1,id,jd,isp),n1,n2,n3,n1,n2,n3,1,iset,-1)

        call zprm3(' gij %i',id+100+jd,gij(1,1,1,id,jd,isp),n1,n2,n3)
C       call zprm3(' gji %i',id+100+jd,gji(1,1,1,id,jd,isp),n1,n2,n3)

   10 continue


      end
