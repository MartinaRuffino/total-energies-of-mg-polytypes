      subroutine spvab0(lmxa,ib,nbas,nel,el,nrad,rad,
     .  vl,h,y,bi,bj,vab0)
C- Makes <J V J> for each pair of energies in one sphere
C  for the spherical part of V. spvab0 complements routine spvls.
      implicit none
      integer lmxa,ib,nbas,nrad,nel
      double precision vab0(9,nbas,1),bi(nrad,0:1),bj(nrad,0:1),
     .  el(1),rad(1),vl(1),h(1),y(1)
      double precision sam,y0
      integer ie,ir,iv,je,l,m,lsp,nsp,ila

      y0 = 1/dsqrt(16*datan(1d0))
      nsp=lsp()+1
      if (nsp /= 1) call rx('spvab0: nsp not spin pol')

      ila = 1
      iv = 0
      do  10  ie = 1, nel
      call ropj(rad,el(ie),lmxa,y,h,bi,nrad)
      do  10  je = ie, nel
        iv = iv+1
        call ropj(rad,el(je),lmxa,y,h,bj,nrad)
        do  20  l = 0, lmxa
          sam = 0d0
          do  22  ir = 1, nrad
   22     sam = sam + bi(ir,l)*bj(ir,l)*vl(ir)*y0
          vab0(l+1,ib,iv) = sam
   20   continue
   10 continue

      end
