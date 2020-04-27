      subroutine rhodif(rnew,rold,nr,nlm,a,rmax,diff)
C- Get charge densities differences in spheres
C ----------------------------------------------------------------------
Ci Inputs: rnew, rold (new and old densities)
Ci         nr, nlm, rmax (sphere properties)
Ci         h work array for Simpson weights
Co Outputs:
Co         diff rms charge difference
Cr Remarks:
Cr         Adapted from mixrhl by removing the linear mix
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nr,nlm
      double precision rnew(nr,nlm,2),rold(nr,nlm,2),diff,
     .                 a,rmax
C Local Variables
      integer ipr,nsp,lsp,isp,lmax,k,ilm,ll,ir,l,intopt,nglob
      double precision charl(11),b,ea,r,wgt,rpb,sam,sum,char,fac,vsph,
     .                 try,diffl(11)
      double precision rofi(nr),wt(nr),h(nr)
      double precision dexp,dsqrt

C --- Setup ---
      call getpr(ipr)
      nsp = lsp() + 1
      lmax = ll(nlm)
      vsph = 4.1887902d0*rmax**3
      do  k = 1,lmax + 1
        diffl(k) = 0.d0
        charl(k) = 0.d0
      enddo
C ... Weights for quadrature on shifted log mesh
      intopt = 10*nglob('lrquad')
      call radmwt(intopt,rmax,a,nr,rofi,wt)
      do  ir = 2, nr
        wt(ir) = wt(ir)/rofi(ir)**2
      enddo

C --- Loop over ilm ---
      do  isp = 1,nsp
        do  ilm = 1,nlm
          l = ll(ilm)
          sum = 0d0
          sam = 0d0
          do  ir = 2, nr
            try = rnew(ir,ilm,isp) - rold(ir,ilm,isp)
            sum = sum + wt(ir)*try*try
            sam = sam + wt(ir)*rnew(ir,ilm,isp)**2
          enddo
          charl(l+1) = charl(l+1) + sam*nsp
          diffl(l+1) = diffl(l+1) + sum*nsp
        enddo
      enddo

C --- Make diff, write out ---
      fac = vsph
C     fac = vsph turns rms density into a "charge" in sphere
      diff = 0d0
      char = 0d0
      do   k = 1, lmax+1
        diff = diff + diffl(k)
        char = char + charl(k)
        diffl(k) = fac*dsqrt(diffl(k)/vsph)
        charl(k) = fac*dsqrt(charl(k)/vsph)
      enddo
      diff = fac*dsqrt(diff/vsph)
      char = fac*dsqrt(char/vsph)
      if(ipr >= 50) then
        write(6,201)
  201   format('   l           rms charge      rms difference')
        write(6,202) (k-1,charl(k),diffl(k),k = 1,lmax + 1)
  202   format(i4,4x,2f17.6)
        write(6,203) char,diff
  203   format(' total  ',2f17.6)
      endif
      if(ipr >= 40) write(6,205) diff/fac,diff
  205 format(' rhodif: rms density diff=',f11.6,
     .  '    charge diff=',f11.6)
      return
      end
