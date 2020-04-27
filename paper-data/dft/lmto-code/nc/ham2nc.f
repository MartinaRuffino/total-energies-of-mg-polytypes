      subroutine ham2nc(nbas,nsp,nl,ldim,ipc,lmx,indxsh,alp,bet,gam,
     .  u,pph,sll,ccor,vmtz,elin,diags,slld,h)
C- Non-collinear two-center hamiltonian (ASA only)
      implicit none
      logical ccor
      integer ldim,nbas,nsp,nl,lmx(*),ipc(nbas),indxsh(*)
      double precision alp(nbas),bet(nbas),gam(nbas),diags(ldim,0:2),
     .  pph(5,ldim,2),sll(ldim,ldim,2),slld(ldim,ldim,2),vmtz,elin,
     .  h(ldim,2,ldim,2,2)
      double complex u(2,2,nbas)
      integer i,j,ib,ic,il,im,lmi,i1,j1,jb,jc,jl,jm,lmj
      double precision xxc,xxs
      double complex xx,uij(2,2)

C --- Spinor rotation matrices for all sites ---
      do  ib = 1, nbas
        xxc = dcos(bet(ib)/2)
        xxs = dsin(bet(ib)/2)
        u(1,1,ib) =  xxc*cdexp(dcmplx(0d0,(alp(ib)+gam(ib))/2))
        u(1,2,ib) =  xxs*cdexp(dcmplx(0d0,(-alp(ib)+gam(ib))/2))
        u(2,1,ib) = -xxs*cdexp(dcmplx(0d0,(alp(ib)-gam(ib))/2))
        u(2,2,ib) =  xxc*cdexp(dcmplx(0d0,(-alp(ib)-gam(ib))/2))
      enddo

C --- Add  srdel U+ Sll U srdel  into hamiltonian ---
      lmi = 0
      do  ib = 1, nbas
        ic = ipc(ib)
        do  il = 0, nl-1
        do  im = -il, il
        lmi = lmi+1
        if (il <= lmx(ic)) then
        i = indxsh(lmi)
        lmj = 0
        do  jb = 1, nbas
          jc = ipc(jb)
          do  i1 = 1, 2
            do  j1 = 1, 2
              uij(i1, j1) = dconjg(u(1,i1,ib))*u(1, j1, jb)
     .                     +dconjg(u(2,i1,ib))*u(2, j1, jb)
            enddo
          enddo
          do  jl = 0, nl-1
            do  jm = -jl, jl
              lmj = lmj+1
              if (jl <= lmx(jc)) then
                j = indxsh(lmj)
C           ... Add srdel * s-dot * (vmtz - elin) srdel into h
                if (ccor) then
                  do  i1 = 1, 2
                    do  j1 = 1, 2
                      xx = dcmplx(sll(i,j,1), sll(i,j,2))
     .                  *(1+(vmtz-elin)*(diags(i,1)+diags(j,1)))
     .                  +dcmplx(slld(i,j,1), slld(i,j,2))*(vmtz-elin)
                      xx = xx*pph(3,i,i1)*pph(3,j,j1)*uij(i1,j1)
                      h(i,i1,j,j1,1) = dble(xx)
                      h(i,i1,j,j1,2) = dimag(xx)
                    enddo
                  enddo
                else
                  do  i1 = 1, 2
                    do  j1 = 1, 2
                      xx = dcmplx(sll(i,j,1),sll(i,j,2))*uij(i1,j1)*
     .                      pph(3,i,i1)*pph(3,j,j1)
                      h(i,i1,j,j1,1) = h(i,i1,j,j1,1) + dble(xx)
                      h(i,i1,j,j1,2) = h(i,i1,j,j1,2) + dimag(xx)
                    enddo
                  enddo
                endif
              endif
            enddo
          enddo
        enddo
        endif
        enddo
        enddo
      enddo

C --- H += C + (vmtz-eln)*<k|k>_constant ---
      call daxpy(ldim,1d0,pph(2,1,1),5,h,2*ldim+1)
      call daxpy(ldim,1d0,pph(2,1,2),5,h(1,2,1,2,1),2*ldim+1)
      if (ccor) then
        do  i = 1, ldim
          h(i,1,i,1,1) = h(i,1,i,1,1)+(vmtz-elin)*diags(i,0)
     .                       *pph(3,i,1)**2
          h(i,2,i,2,1) = h(i,2,i,2,1)+(vmtz-elin)*diags(i,0)
     .                       *pph(3,i,2)**2
        enddo
      endif
C     call prmx('h at end of ham2nc',h,ldim*2,ldim*2,ldim*2)
      end

      subroutine prmz(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
      implicit none
      integer nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(10) fmt, strn*(*)
      integer i,j,fopna
      fmt = '(9f22.17)'
      ifi = fopna('out',29,0)
      write(ifi,*) nr, nc
      do  i = 1, nr
        write (ifi, fmt) (s(1,i,j), j=1, nc)
      enddo
      write(ifi,*)
      do  i = 1, nr
        write (ifi, fmt) (s(2,i,j), j=1, nc)
      enddo
      call fclose(ifi)
      print *, 'prm: pausing after writing data ',strn
C     pause
      end
