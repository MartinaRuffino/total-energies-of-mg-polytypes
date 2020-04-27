      subroutine gfefbz(nl,nbas,isp,nsp,pos,offH,indxsh,istab,nk1,nk2,
     .  nk3,ipq,offi,offj,ni,nj,ldim,nbmx,eval,evec,g,ag,igstar,ifac,qb,
     .  h1,ze,ldiag,Gij,Gji)
C- GF subblock in full BZ from eigenvectors in the irreducible part
Ci   ldiag:  T, do not calculate Gji
      implicit none
      logical ldiag
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nbas,nl,isp,nsp,ldim,offH(n0H,nkap0,1),indxsh(*),nk1,nk2,
     .  nk3,ni,nj,istab(nbas,1),nbmx,ipq(1),igstar(0:1),offi,offj,
     .  ifac(3)
      double precision pos(3,*),g(3,3,*),ag(3,*),h1(ldim,ldim,2),qb(3,3)
      double precision evec(ldim,ldim,2,nsp,*),eval(nbmx,nsp,*)
      double complex Gij(nk1,nk2,nk3,ni,nj,isp),ze,
     .               Gji(nk1,nk2,nk3,nj,ni,isp)
C Local variables
      integer iq1,i1,i2,i3,k,jj1,jj2,jj3,ipq1,igrp1,nlx,i,j,nu
      double complex wk
      parameter (nlx=4)
      double precision q1(3),qk,rmat(nlx**4)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C     call yprm('1,1',2,evec(1,1,1,1,1),ldim*ldim,ldim,ldim,ldim)

      iq1 = 0
      do  10  i3 = 1, nk3
      do  10  i2 = 1, nk2
      do  10  i1 = 1, nk1

        iq1 = iq1+1
        ipq1 = ipq(iq1)
        igrp1 = igstar(iq1)
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

        if (igrp1 /= 1) then
          call rotwf(0,nl,nbas,1,pos,offH,indxsh,istab(1,igrp1),
     .      g(1,1,igrp1),ag(1,igrp1),q1,rmat,0,ldim,ldim,ldim,
     .      evec(1,1,1,isp,ipq1),h1)
        else
          call dcopy(ldim*ldim*2,evec(1,1,1,isp,ipq1),1,h1,1)
        endif

        do  20  i = offi+1, offi+ni
        do  20  j = offj+1, offj+nj
          wk = 0
          do  22  nu = 1, ldim
            wk = wk + (dcmplx(h1(i,nu,1),h1(i,nu,2)) *
     .                 dcmplx(h1(j,nu,1),-h1(j,nu,2))) /
     .                (Ze-eval(nu,isp,ipq1))
   22     continue
C         if (i+j == 2) print *, i1,i2,i3,wk
          Gij(i1,i2,i3,i-offi,j-offj,isp) = wk
          if (.not. ldiag) then
            wk = 0
            do  24  nu = 1, ldim
              wk = wk + (dcmplx(h1(j,nu,1),h1(j,nu,2)) *
     .                   dcmplx(h1(i,nu,1),-h1(i,nu,2))) /
     .                  (Ze-eval(nu,isp,ipq1))
   24       continue
            Gji(i1,i2,i3,j-offj,i-offi,isp) = wk
          endif
   20   continue

C      i = nk1*nk2*nk3
C      call zcopy(ni*nj,gij(i1,i2,i3,1,1,isp),i,h1,1)
C      call zprm('gij',2,h1,ni,ni,nj)

   10 continue

C     call zprm3('gij(1,1)',0,gij,nk1,nk2,nk3)
C     call zprm3('gij(1,10)',0,gij(1,1,1,1,10,isp),nk1,nk2,nk3)

      end
