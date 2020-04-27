      subroutine paugb(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2,
     .  lbf,sodb,nbf,bfield,nlx1,nlx2,ppiz)
C- Add to ppi constribution from external magnetic field
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of 'bra' function types for each l
Ci   nf1s  :number of 'bra' function types for each l, for which there
Ci         :is a smooth part to be subtracted (which also corresponds
Ci         :to the functions which connect to envelope functions
Ci         :corresponds to the functions which connect to envelope
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   v1    :values of f1 at rofi(nr) (not multiplied by r)
Ci   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
Ci   nf2   :number of 'ket' function types for each l
Ci   nf2s  :number of 'bra' function types for each l, for which there
Ci         :is a smooth part to be subtracted (which also corresponds
Ci         :to the functions which connect to envelope functions
Ci         :corresponds to the functions which connect to envelope
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   v2    :values of f2 at rofi(nr) (not multiplied by r)
Ci   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
Ci   sodb  :matrix elements of radial functions in unit potential;
Ci         :see potpus.f
Ci         :sodb(:,:,:,1) = matrix elements for spin diagonal parts
Ci         :sodb(:,:,:,2) = matrix elements for spin off diagonal parts
Ci         :sodb(1,l+1,isp,:)=<u|SO|u>
Ci         :sodb(2,l+1,isp,:)=<u|SO|s>
Ci         :sodb(3,l+1,isp,:)=<s|SO|u>
Ci         :sodb(4,l+1,isp,:)=<s|SO|s>
Ci         :sodb(5,l+1,isp,:)=<u|SO|z>
Ci         :sodb(6,l+1,isp,:)=<s|SO|z>
Ci         :sodb(7,l+1,isp,:)=<z|SO|z>
Ci         :sodb(8,l+1,isp,:)=<z|SO|u>
Ci         :sodb(9,l+1,isp,:)=<z|SO|s>
Ci   nbf   :dimensions bfield
Ci   bfield :external magnetic field for this site
Ci   nlx1  :dimensions ppiz
Ci   nlx2  :dimensions ppiz
Co Outputs
Co   ppiz  :Contribution from bfield is added into ppiz
Co         :This routine does not initialize ppiz!
Cr Remarks
Cr   Adds B.sigma*sodb to ppiz
Cr   for B constant in space, but orbital dependent.
Cr   This is the potential acting on electrons (spin 1/2), since:
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   Electron spin operator is S = 1/2 sigma or (hbar/2 sigma)
Cr   Coupling is g0 S B => E = (g0.ms.muB) B
Cr   Splitting between ms=1/2 and ms=-1/2 is  g0 muB B = 2 muB B.
Cu Updates
Cu   17 Jul 13 Bug fix for local orbitals case
Cu   01 Apr 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmx1,lmx2,nf1,nf1s,nf2,nf2s,nlx1,nlx2,nbf,n0,nab,lbf
      integer lx1(nf1),lx2(nf2)
      parameter (n0=10,nab=9)
      double precision sodb(nab,n0,2,2),bfield(nbf,3),
     .  v1(0:lmx1,nf1),d1(0:lmx1,nf1),
     .  v2(0:lmx2,nf2),d2(0:lmx2,nf2)
      double complex ppiz(nf1,nf2,nlx1,nlx2,2,2)
C ... Local parameters
      integer i1,i2,lmax1,lmax2,lmax,l,l1,m1,ibf,isp,k
      double precision facB,tmp1,tmp2
      double complex bdots(2,2)
C     facB covers convention for B:
C     facB= 1 => +B induces positive -M
C     facB=-1 => +B induces positive +M
      parameter (facB=1d0)

      if (lbf == 0) return

C ... Part involving envelope functions
      do  i1 = 1, nf1s
        do  i2 = 1, nf2s
          lmax1 = lx1(i1)
          lmax2 = lx2(i2)
          lmax = min0(lmax1,lmax2)
          l1 = 0
          do  l = 0, lmax
          ibf = 1
          if (nbf > 1) ibf = l
          k = l+1
          bdots(1,1) =  bfield(ibf,3)                         *facB
          bdots(1,2) =  dcmplx(bfield(ibf,1),-bfield(ibf,2))  *facB
          bdots(2,1) =  dcmplx(bfield(ibf,1), bfield(ibf,2))  *facB
          bdots(2,2) = -bfield(ibf,3)                         *facB
          do  m1 = -l, l
            l1 = l1 + 1
            do  isp = 1, 2
              tmp1 =
     .          (v1(l,i1)*sodb(1,k,isp,1)*v2(l,i2)
     .         + v1(l,i1)*sodb(2,k,isp,1)*d2(l,i2)
     .         + d1(l,i1)*sodb(3,k,isp,1)*v2(l,i2)
     .         + d1(l,i1)*sodb(4,k,isp,1)*d2(l,i2))
              tmp2 =
     .          (v1(l,i1)*sodb(1,k,isp,2)*v2(l,i2)
     .         + v1(l,i1)*sodb(2,k,isp,2)*d2(l,i2)
     .         + d1(l,i1)*sodb(3,k,isp,2)*v2(l,i2)
     .         + d1(l,i1)*sodb(4,k,isp,2)*d2(l,i2))
              ppiz(i1,i2,l1,l1,isp,1) = ppiz(i1,i2,l1,l1,isp,1) +
     .          bdots(isp,isp) * tmp1
              if (lbf == 1) then
              ppiz(i1,i2,l1,l1,isp,2) = ppiz(i1,i2,l1,l1,isp,2) +
     .          bdots(isp,3-isp) * tmp2
              endif
            enddo
          enddo
          enddo
        enddo
      enddo

C ... Local orbitals part
      do  i1 = 1, nf1
        do  i2 = 1, nf2
          if (i1 > nf1s .or. i2 > nf2s) then
          lmax1 = lx1(i1)
          lmax2 = lx2(i2)
          lmax = min0(lmax1,lmax2)
          l1 = 0
          do  l = 0, lmax
          ibf = 1
          if (nbf > 1) ibf = l
          k = l+1
          bdots(1,1) =  bfield(ibf,3)                         *facB
          bdots(1,2) =  dcmplx(bfield(ibf,1),-bfield(ibf,2))  *facB
          bdots(2,1) =  dcmplx(bfield(ibf,1), bfield(ibf,2))  *facB
          bdots(2,2) = -bfield(ibf,3)                         *facB
          do  m1 = -l, l
            l1 = l1 + 1
            do  isp = 1, 2
C         ... hso_zz
              if (i1 > nf1s .and. i2 > nf2s) then
                tmp1 = sodb(7,k,isp,1)
                tmp2 = sodb(7,k,isp,2)
C         ... hso_zu
              elseif (i1 > nf1s) then
                tmp1 = sodb(8,k,isp,1)*v2(l,i2)+sodb(9,k,isp,1)*d2(l,i2)
                tmp2 = sodb(8,k,isp,2)*v2(l,i2)+sodb(9,k,isp,2)*d2(l,i2)
C         ... hso_uz
              elseif (i2 > nf2s) then
                tmp1 = v1(l,i1)*sodb(5,k,isp,1)+d1(l,i1)*sodb(6,k,isp,1)
                tmp2 = v1(l,i1)*sodb(5,k,isp,2)+d1(l,i1)*sodb(6,k,isp,1)
              endif
              ppiz(i1,i2,l1,l1,isp,1) = ppiz(i1,i2,l1,l1,isp,1) +
     .          bdots(isp,isp) * tmp1
              if (lbf == 1) then
              ppiz(i1,i2,l1,l1,isp,2) = ppiz(i1,i2,l1,l1,isp,2) +
     .          bdots(isp,3-isp) * tmp2
              endif
            enddo
          enddo
          enddo
          endif
        enddo
      enddo

      end
