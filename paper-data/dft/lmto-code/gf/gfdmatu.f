      subroutine gfdmatu(nbas,nsp,nspc,lmaxu,nl,nlibu,idu,dmat,dmatu)
C- Extract LDA+U output density from dmat (ldmat=1) (gfasa.f)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   lmaxu :dimensioning parameter for U matrix
Ci   nl    :(global maximum l) + 1
Ci   nlibu :number of LDA+U blocks
Ci   dmat  :output density matrix for all (RL,RL')
Ci   dmatu :output density matrix for LDA+U
Co Outputs
Co   dmatu :output density matrix for LDA+U blocks
Cu Updates
Cu   07 Sep 13 (Belashchenko) dmat is now complex
Cu   08 Nov 07 (J. Xu) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nsp,lmaxu,nl,nlibu
      integer nspc,idu(4,nbas)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex dmat(nl**2,nl**2,nbas,nsp,nspc)
C ... Local parameters
      integer isp,jsp,ksp,ib,l,m1,m2,lcount,iblu,id

C     Collinear case:  spin diagonal part of dmat is dmat(..,isp,1)
C     Noncollinear case: "    "      "    "  " is  dmat(..,isp,isp)
C     ksp = first spin index for spin diagonal part of dmat
C     jsp = second spin index for spin diagonal part of dmat
      do  isp = 1, nsp
      do  jsp = 1, nspc
      ksp = max(jsp,isp)

        iblu = 0
        do  ib = 1, nbas
          id = 0
C         m1 and m2 zero => no LDA+U block in ib
          call imxmn(4,idu(1,ib),1,m1,m2)
          if (m1 /= 0 .or. m2 /= 0) then
            lcount = 0
            do  l = 0, nl-1
              if (idu(l+1,ib) /= 0) then
                iblu = iblu + 1
                do  m1 = -l, l
                  do  m2 = -l, l
                    dmatu(m1,m2,ksp,iblu)=
     .              real(dmat(id+m1+l+1,id+m2+l+1,ib,ksp,jsp))
                  end do
                end do
                lcount = lcount + 1
              end if
              id = id + 2*l + 1
            end do
          end if
        end do
      end do
      end do

      end
