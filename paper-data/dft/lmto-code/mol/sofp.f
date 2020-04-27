      subroutine sofp(mode,phiz,phi,phid,nr,nsp,lmxs,lmx,v,dv,ez,enu,z,
     .  ri,wi,wk,sop,sopz)
C- Radial matrix elements between orbitals of different spin
C ---------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make spin-orbit parameters, i.e. matrix elements
Ci         :only <phi|so|phi> <phi|so|phidot> <phidot|so|phidot>
Ci         :2 make spin-orbit parameters, i.e. matrix elements
Ci         :  <phi|so|phi> <phi|so|phidot> <phidot|so|phidot>
Ci         :   and
Ci         :  <phiz|so|phiz> <phiz|so|phi> <phiz|so|phidot>
Ci         :4 orthonormalize phi,phidot in separate spin channels
Ci         :5 1+4 above
Ci         :6 2+4 above
Ci   phi   :radial wave function * r
Ci   phid  :energy derivative of phi
Ci   phiz  :local orbital wave function
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmxs  :leading dimension of sop is 0:lmxs
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   v     :electron part of potential: complete potential is v(r)-2*z/r.
Ci   dv    :(mode=1) radial derivative of potential, dV/dr with V=(V+ + V-)/2
Ci         :(mode=2) magnetic field
Ci         :(mode=4) not used
Ci   enu   :enu's for making charge density
Ci   ez    :enu's for local orbitals
Ci   z     :nuclear charge
Ci   ri    :radial mesh
Ci   wi    :weights for integration on mesh
Ci   wk    :work array of length nr*4
Co  Outputs
Co   sop   :sop(l,is1,is2,i=1..3) : matrix elements between orbitals
Co         :of spin is1 and is2 for quantum number l.
Co         :Three types of integrals are calculated for i=1..3:
Co         :<phi SO phi>  <phi SO phidot>  <phidot SO phidot>
Co   sopz  :sopz(l,is1,is2,i=1..3) : matrix elements between local orbitals
Co         :and orbitals of spin is1 and is2 for quantum number l.
Co         :Three types of integrals are calculated for i=1..3:
Co         :<phiz SO phiz>  <phiz SO phi>  <phiz SO phidot>.
Cr  Remarks
Cr   Adapted from soprm.f, instead of calculating the <|B|>, it
Cr   calculates the additional matrix elements of so,
Cr   when local orbitals are added to the basis.
Cr   so = 2/(c^2) dV/dr*(1/r), V(r)=-2*z/r+v(r)
Cr   Note: so=1/(2*m^2*c^2)*(dV/dr*1/r), m=.5, c=274 (at. Rydberg units)
Cr   H_so = so*s^ dot l^, s^=0.5d0*sigma (Pauli matrix).
Cb Bugs
Cb   should be merged with soprm
Cu Updates
Cu   25 Apr 05 A. Chantis adapted from soprm
C ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmxs,mode(0:lmxs),nr,lmx,nsp
      double precision phiz(nr,0:lmxs,nsp),phi(nr,0:lmxs,nsp),
     .  phid(nr,0:lmxs,nsp),ri(nr),
     .  wi(nr),z,sop(0:lmxs,nsp,nsp,3),v(nr,nsp),wk(nr,4),dv(nr),
     .  enu(0:8,nsp),ez(0:8,nsp),sopz(0:lmxs,nsp,nsp,3)
C ... Local parameters
      integer l,ir,is,is1,is2,ipr,mode0(0:lmx),stdo,lgunit,lmin
      double precision c,pa,r,r1,r2,dot3,vavg,eavg,eavgz,dva,xx,
     .  xxz,xxavg,wkz(nr,4)
C ... External calls
      external daxpy,dscal,getpr
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c
      data pa /1d0/

C --- Setup ---
      call getpr(ipr)
      stdo = lgunit(1)
      do l = 0, lmx
      mode0(l) = mod(mode(l),4)
      enddo
C     c = 274.071979d0
      lmin = 1

C --- Orthonormalize phi, phidot, neglecting small component ---

      if (ipr > 50)
     .  print '(/'' sofp: overlaps  phi*phi     phi*phidot'')'
      do  1  is = 1, nsp
      do  1  l = 0, lmx
        r1 = dot3(nr,phi(1,l,is),phi(1,l,is),wi)
        call dscal(nr,1/dsqrt(r1),phi(1,l,is),1)
        r2 = dot3(nr,phi(1,l,is),phid(1,l,is),wi)
        call daxpy(nr,-r2,phi(1,l,is),1,phid(1,l,is),1)
        if (ipr > 50) write(stdo,334) is,l,r1,r2
  334   format('  spin',i2,'  l=',i1,2f13.6)
    1 continue

C --- Matrix elements for each l ---
      do is = 1, 4
      wk(1,is) = 0d0
      wkz(1,is) = 0d0
      enddo

      do  10  l = 0, lmx
         eavg = (enu(l,1)+enu(l,nsp))/2
         if (mode0(l) == 2) eavgz = (ez(l,1)+ez(l,nsp))/2
         do  12  is2 = 1, nsp
         do  12  is1 = 1, nsp
          if (mode0(l) == 1 .or. mode0(l) == 2) then
             do  14  ir = 2, nr
               r = ri(ir)
               vavg = (v(ir,1)+v(ir,nsp))/2 - 2*z/r
               dva  = dv(ir) + 2*z/r**2
               xx = 1/r/(1d0+pa*(eavg-vavg)/c**2)**2
               if (mode0(l) == 2) then
               xxz = 1/r/(1d0+pa*(eavgz-vavg)/c**2)**2
               xxavg = 0.5d0*(xx+xxz)
               wkz(ir,1)=phiz(ir,l,is1)*dva
               wkz(ir,2)=phiz(ir,l,is1)*xxz
               wkz(ir,3)=phi(ir,l,is2)*xxavg
               wkz(ir,4)=phid(ir,l,is2)*xxavg
               endif
               wk(ir,1)=phi(ir,l,is1)*dva
               wk(ir,3)=phid(ir,l,is1)*dva
               wk(ir,2)=phi(ir,l,is2)*xx
               wk(ir,4)=phid(ir,l,is2)*xx
   14        continue
             sop(l,is1,is2,1)=dot3(nr,wk,wk(1,2),wi)*2d0/c**2
             sop(l,is1,is2,2)=dot3(nr,wk,wk(1,4),wi)*2d0/c**2
             sop(l,is1,is2,3)=dot3(nr,wk(1,3),wk(1,4),wi)*2d0/c**2
           if (mode0(l) == 2) then
             sopz(l,is1,is2,1)=dot3(nr,wkz,wkz(1,2),wi)*2d0/c**2
             sopz(l,is1,is2,2)=dot3(nr,wkz,wkz(1,3),wi)*2d0/c**2
             sopz(l,is1,is2,3)=dot3(nr,wkz,wkz(1,4),wi)*2d0/c**2
           endif
           else
             call rxi('sofp: bad mode:',mode)
           endif

   12     continue
   10 continue

C --- Printout ---
      if (ipr <= 50) return
      write(stdo,332) 'spin-orbit coupling'
  332 format(' sofp: matrix elements for perturbation from ',a/
     .  13x,'l',5x,'<phi || phi>',8x,'<dot || phi>',8x,'<dot || dot>')
      if (nsp == 1) then
        do  22  l = lmin, lmx
        write(stdo,333) '          ',
     .      l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
        if (mode0(l) == 2) then
       write(stdo,333) '          ',
     .      l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
         endif
 22      continue
      else
        do  20  l = lmin, lmx
        write(stdo,333) 'up   up   ',
     .      l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
        write(stdo,333) 'down down ',
     .    l,sop(l,2,2,1),sop(l,2,2,2),sop(l,2,2,3)
        write(stdo,333) 'up   down ',
     .    l,sop(l,1,2,1),sop(l,1,2,2),sop(l,1,2,3)
        write(stdo,333) 'down up   ',
     .    l,sop(l,2,1,1),sop(l,2,1,2),sop(l,2,1,3)
        write(stdo,333)
        if (mode0(l) == 2) then
        write(stdo,333) 'up   up   ',
     .      l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
        write(stdo,333) 'down down ',
     .    l,sopz(l,2,2,1),sopz(l,2,2,2),sopz(l,2,2,3)
        write(stdo,333) 'up   down ',
     .    l,sopz(l,1,2,1),sopz(l,1,2,2),sopz(l,1,2,3)
        write(stdo,333) 'down up   ',
     .    l,sopz(l,2,1,1),sopz(l,2,1,2),sopz(l,2,1,3)
        write(stdo,333)
        endif
   20   continue
      endif
  333 format(1x,a,i3,3f20.15)

      end

