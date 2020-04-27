      subroutine tetwtt(job,nfilox,nfiupx,nemlox,nemupx,npm,n1,n2,n3,
     .  nqbz,ipq,nwhis,nwgtx,ibjb,nhwtot,ihw,nhw,jhw,whw,nfilo,nfiup,
     .  nemlo,nemup,nfilmw,nempw,nsp,optmt,njdosw,jdosw,npol,nkp,
     .  wtthis,wtphis)
C- Combine partial DOS weights for total JDOS or optical matrix element
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0  Make total joint DOS, wtthis
Ci         :1  Make joint DOS * optmt
Ci         :2  Make joint DOS with Mulliken weighting factors
Ci         :10s digit
Ci         :1  Make wtphis, resolves wtthis into (occ,unocc) parts
Ci         :2  Make wtphis, resolves wtthis into k contributions
Ci         :3  Make wtphis, resolves wtthis both by (occ,unocc) and k
Ci   nfilox:first occupied band for which transitions are indexed by ibib
Ci         :ibjb(1,:,:,:) points to band nfilox.
Ci         :Note: it is an error for nfilox>nfilo
Ci   nfiupx:last  occupied band for which transitions are indexed by ibib
Ci   nemlox:first unoccupied band for which transitions are indexed by ibib
Ci         :ibjb(:,1,:,:) points to band nfilox.
Ci         :Note: it is an error for nemlox>nemlo
Ci   nemupx:last  unoccupied band for which transitions are indexed by ibib
Ci         :It may be less than nemup because all transitions for jb>nemupx
Ci         :lie outside the energy window.
Ci   npm   :1 : use time reversal symmetry
Ci         :2 : no time reversal symmetry
Ci   n1..n3:no. of divisions made along each reciprocal lattice vector
Ci   nqbz  :number of k-points in the 1st BZ
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   nwhis :number of energies in histogram
Ci   nwgtx :leading dimension of ihw,jhw,nhw
Ci   ibjb  :ijb = ibjb(ib,jb,k) = index to given ib,jb pair for this k point.
Ci         :in ihw,nhuw,jw,whw below
Ci   nhwtot:dimension of whw
Ci   ihw   :ihw(ijb,k) = index to first histogram bin encompassing
Ci                       (demin,demax) for a given ib,jb pair and k; see whw
Ci   nhw   :nhw(ijb,k) = number of histogram bins encompassing
Ci                       (demin,demax) for a given ib,jb pair and k
Ci   jhw   :jhw(ijb,kx) = index to first bin in whw for this pair; see whw
Ci   whw   :whw(i:j) histogram weights in bins i:j for given ib,jb,kx
Ci         :i = jhw(ijb,kx)
Ci         :j = jhw(ijb,kx)+nhw(ijb),kx)-1
Ci         :ibjb = ibjb(ib,jb,kx) = index to iwh,jhw,nhw for given
Ci         :(ib,jb) pair at kx
Ci   nfilo :dimensions, and sets lower occ range for optmt
Ci   nfiup :dimensions, and sets higher occ range for optmt
Ci   nemlo :dimensions, and sets lower unocc range for optmt
Ci   nemup :dimensions, and sets higher unocc range for optmt
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   optmt :matrix element <i | grad | j>
Ci   njdosw:number of color weights
Ci   jdosw :color weights for states i and j
Ci   npol  :lopt=T => npol=3
Ci         :Joint DOS w/out color weights => npol=1
Ci         :Joint DOS w/   color weights => npol=1+njdosw
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Co Outputs
Co   wtthis:is accumulated from whw; see Remarks
Co   wtphis:
Cl Local variables
Ci   lopt  :T => npol=3
Cr Remarks
Cr   This routine accumulates by the tetrahedron method this integral:
Cr     wtthis(ihis) = \int_hislow^hisup d\omega \times
Cr                    \int d^3k f(e_i(k))(1-f(e_j(q+k))) \times
Cr                              \delta(omg - e_j(q+k) + e_i(k)) \ times P
Cr   where
Cr     f(E) = Fermi distribution function
Cr     P    = 1 for joint DOS
Cr          = optmt(:,i,j,ik) if lopt is .true.
Cr
Cr   wtthis is accumulated from whw, which resolves
Cr   wtthis into (ib,jb) pairs and individual k-contributions.
Cr
Cu Updates
Cu   07 Sep 09  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nfilox,nfiupx,nemlox,nemupx,nfilmw,nempw,njdosw
      integer nfilo,nfiup,nemlo,nemup,nsp,npol,nkp,n1,n2,n3
      integer job,nqbz,npm,ipq(nqbz)
      integer nwgtx,nhwtot,nwhis
      double precision wtthis(nwhis,npol),whw(nhwtot)
      double precision jdosw(njdosw,nfilo:nfiup+nemup-nemlo+1,nkp)
      double precision optmt(3,nfilo:nfiup,nemlo:nemup,nsp,nkp),
     .                 wtphis(nwhis,npol,nfilmw,nempw,nkp)
      integer ihw(nwgtx,nqbz,npm), ! omega pointer
     .        nhw(nwgtx,nqbz,npm), ! number of bins
     .        jhw(nwgtx,nqbz,npm), ! index to whw(*)
     .        ibjb(nfilox:nfiupx,nemlox:nemupx,nqbz,npm)
C ... Local parameters
      logical lopt,lmull
      integer iq1,ipol,ib,jb,jpm,lpart,i1,i2,i3,iq,ibw,jbw,iqw
      integer ini,ied,ioff,ijpr
      double precision fac1,fac2

      if (n1*n2*n3 /= nqbz) call rx('mismatch n1*n2*n3 with nqbz')
      lpart = mod(job/10,10)
      lopt  = mod(job,10) == 1
      lmull = mod(job,10) == 2
      if (lopt .and. npol /= 3) then
        call rxi('tetwtt: illegal value, npol =',npol)
      elseif (lmull .and. npol /= njdosw+1) then
        call rxi('tetwtt: illegal value, npol =',npol)
      endif
      wtthis = 0

C ... Loop over iq1=1...n1*n2*n3 qp, extracting irr iq for each point.
C     Though the loops execute nkp*n1*n2*n3 times, only n1*n2*n3
C     survive the if(ipq(iq1) == iq) conditional
      do  iq = 1, nkp
        iq1 = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
        iq1 = iq1+1
        if (ipq(iq1) == iq) then

        do  ib = nfilo, nfiup
        do  jb = nemlo, nemup
        do  jpm = 1, npm

C       No index to this ibjb pair ... skip
          if (ib < nfilox .or. ib > nfiupx .or.
     .      jb < nemlox .or. jb > nemupx) then
            cycle
          endif
          ibw = ib-nfilo+1     ! index to wtphis
          jbw = jb-nemlo+1     ! index to wtphis
          if (lpart == 2) then
            ibw = 1 ; jbw = 1
          endif
          iqw = iq; if (lpart < 2) iqw = 1

          ijpr = ibjb(ib,jb,iq1,jpm)
          if (ijpr == 0) cycle
          ini = ihw(ijpr,iq1,jpm)
          ied = ihw(ijpr,iq1,jpm) + nhw(ijpr,iq1,jpm)-1
          ioff= jhw(ijpr,iq1,jpm)
          if (ied > nwhis)
     .    call rx('something wrong with dimensioning of wtthis')
          if (lopt) then
            do  ipol = 1, npol
              wtthis(ini:ied,ipol) = wtthis(ini:ied,ipol) +
     .          whw(ioff:ioff+nhw(ijpr,iq1,jpm)-1) *
     .          optmt(ipol,ib,jb,1,iq)
              if (lpart > 0) then
                wtphis(ini:ied,ipol,ibw,jbw,iqw) =
     .          wtphis(ini:ied,ipol,ibw,jbw,iqw) +
     .            whw(ioff:ioff+nhw(ijpr,iq1,jpm)-1) *
     .            optmt(ipol,ib,jb,1,iq)
              endif
            enddo
          else
          wtthis(ini:ied,1) = wtthis(ini:ied,1) +
     .      whw(ioff:ioff+nhw(ijpr,iq1,jpm)-1)
          if (lmull) then
            do  ipol = 1, njdosw
              fac1 = jdosw(ipol,ib,iq)
              fac2 = jdosw(ipol,nfiup+jb-nemlo+1,iq)
              wtthis(ini:ied,1+ipol) = wtthis(ini:ied,1+ipol) +
     .          whw(ioff:ioff+nhw(ijpr,iq1,jpm)-1) * fac1*fac2
C              if (iq == 1) print 544, ib,jb,fac1,fac2
C  544         format(2i4,2f12.5)
            enddo
          endif
          if (lpart > 0) then
            wtphis(ini:ied,1,ibw,jbw,iqw) =
     .      wtphis(ini:ied,1,ibw,jbw,iqw) +
     .        whw(ioff:ioff+nhw(ijpr,iq1,jpm)-1)
          endif
          if (lmull) then
            do  ipol = 1, njdosw
              fac1 = jdosw(ipol,ib,iq)
              fac2 = jdosw(ipol,nfiup+jb-nemlo+1,iq)
              if (lpart > 0) then
                wtphis(ini:ied,1+ipol,ibw,jbw,iqw) =
     .          wtphis(ini:ied,1+ipol,ibw,jbw,iqw) +
     .            whw(ioff:ioff+nhw(ijpr,iq1,jpm)-1) * fac1*fac2
              endif
            enddo
          endif
          endif
        enddo                   ! end of jpm loop
        enddo                   ! end of jb loop
        enddo                   ! end of ib loop
        endif

        enddo ! end of i1 loop
        enddo ! end of i2 loop
        enddo ! end of i3 loop
      enddo                     ! end of iq loop

      end
