      subroutine makdia(nl,nbas,ndim,indxsh,lidim,ipc,
     .                  wsr,avw,alpha,adot,ccd)
C- make constant, linear and bilinear diagonal matrices for combined corr.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ndim  :leading dimension of ccd
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   lidim :size lower+downfolded block
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   wsr   :Wigner-Seitz radius, in a.u. (input; alias rmax)
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   alpha :tight-binding screening parameters
Ci   adot  :(kappa*avw)^2-derivative of tight-binding screening const.
Co Outputs
Co   ccd   :ccd(*,n) are the nth-order terms in S^alpha that appear
Co          in eq. 3.87 Kanpur notes, scaled by 1/avw**2.
Cr Remarks
Cu Updates
Cu   14 Sep 99 discard uneeded argument lmx; use lidim
C-----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nl,nbas,ndim,indxsh(*),lidim,ipc(nbas)
      double precision alpha(*),adot(*),ccd(ndim,0:2),wsr(*),avw
C Local parameters
      integer ibas,l,l2p3,l2m1,l2p1,ml2m1,m,lmr
      double precision a,sw,dl,t1,t2,t3

      lmr = 0
      do  ibas = 1, nbas
        do  l = 0, nl-1
          do  m = -l, l
            lmr = lmr + 1
            if (indxsh(lmr) <= lidim) then
            ccd(indxsh(lmr),2) = adot(lmr)
            endif
          enddo
        enddo
      enddo

      lmr = 0
      do  ibas = 1, nbas
        sw = wsr(ipc(ibas))/avw
        do  l = 0, nl-1
          l2p3 = 2*l + 3
          l2p1 = 2*l + 1
          l2m1 = 2*l - 1
          ml2m1 = -l2m1
          dl = sw**ml2m1 / ( .5d0 * l2m1 )
          do  m = -l, l
            lmr = lmr + 1
            if (indxsh(lmr) <= lidim) then
            a = alpha(lmr)
            t1 = -sw**l2p3 / (2d0 * l2p1*l2p1 * l2p3)
            t2 = a * sw*sw / l2p1
            t3 = a*a * 2d0 * sw**ml2m1 / l2m1
            ccd(indxsh(lmr),2) = ccd(indxsh(lmr),2) + t1 + t2 + t3
            ccd(indxsh(lmr),1) = sw*sw / ( 2d0 * l2p1 ) + a * dl
            ccd(indxsh(lmr),0) = dl
            endif
          enddo
        enddo
      enddo
      end
