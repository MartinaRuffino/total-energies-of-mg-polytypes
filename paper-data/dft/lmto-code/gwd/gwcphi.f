      subroutine gwcphi(s_site,s_spec,isp,nsp,nlmax,ndham,nphimx,nev,nbas,
     .  ipb,lmxax,nlindx,ndima,ppnl,aus,cphi,cphin)
C- Reorder (phi,phidot) project onto MT sphere for GW input; make overlap
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  p pz lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus : must be at least hamiltonian rank
Ci   nphimx:dmensions aus: global max # of partial waves of a given l
Ci   nev   :number of eigenvectors to accumulate cphi
Ci   nbas  :number of sites in list
Ci   ipb   :index to true basis (excluding floating orbitals)
Ci         :given site index including those orbitals
Ci   lmxax :dimensions nlindx
Ci   nlindx:nlindx(3,0:lmxax,1:nat) hold offsets to channels in cphi.
Ci         :wf is projected onto a sum of (phi,phidot,LO) for each l and m.
Ci         :1st index is for phi,phidot,LO, 2nd index ranges over augmentation l
Ci         :Sites 1:nat, i.e. 1..nbas but excluding "floating orbital" sites
Ci         :Let ichan = nlindx(i,l,ia)  Channels for phi (i=1), dot(i=2) or phiz (i=3)
Ci         :are contained in (ichan+1:ichan+2*l+1) for given l, ia, with m=1..2l+1
Ci   ndima :leading dimension of cphi
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci         :This routine uses, in ppnl(i,1:lmxa+1,1:nsp,1:nbas),
Ci         :ppnl(2) = s00 = <phi | phi>
Ci         :ppnl(7) = s11 = <phidot | phidot>
Ci         :ppnl(8)  = szz = <gz|gz> for local orbital
Ci         :ppnl(9)  = overlap <phi|gz> = s0z
Ci         :ppnl(10) = overlap <phidot|gz> = s1z
Ci   aus   :values of (phi,phidot) at MT sphere boundary; see makusq
Co Outputs
Co   cphi :coefficients to phi,phidot,phiz following Kotani conventions
Co        :cphi(ichan,iv) :
Co           ichan = orbital channel
Co           for a given site, l, and m and one of (phi,dot,phiz); see nlindx.
Co           Channels ordered by first by L, then by site, then by orbital-type
Co           iv = eigenvector
Co   cphin:diagonal matrix elements, one for each eigenvector.
Co        :cphin(1,iv) = <cphi(iv) | overlap | cphi(iv)>
Co        :cphin(2,iv) = <cphi(iv) | overlap-from-phi-only | cphi(iv)>
Cr Remarks
Cr   gwcphi converts the projection of an eigenvector into coefficients
Cr   of (phi,phidot) pairs, or with local orbitals, (phi,phidot,phiz)
Cr   into format for GW codes.  The overlap matrix between
Cr   the original functions is (see ppnl above)
Cr
Cr             (s00   0      s0z  )
Cr      S =  = ( 0   s11     s1z  )
Cr             (s0z  s1z     szz  )
Cr
Cr *Linear transformation of augmentation orbitals.  This code was
Cr  originally designed to make a linear transformation of atomic
Cr  orbitals (phi,phidot,phiz) to an orthonormal set, and generate
Cr  cphi for the orthonormalized set.  Now it generates cphi for the
Cr  original (non-orthonormal) orbitals.  To compute the sphere
Cr  contribution to the normalization, the sphere contribution to matrix
Cr  elements between eigenvectors is returned in cphin.
Cr
Cr  Vector aus corresponds to coefficients of augmented functions
Cr  {phi,phidot,gz}, which are not orthnormal.  Let us choose a linear
Cr  transformation L that transforms new functions {f1~..fn~} back to
Cr  the original functions {f1...fn}.  L^-1 transforms {f1...fn} to
Cr  {f1~..fn~}.
Cr
Cr  A sum of functions sum_i C_i f_i can be expressed as
Cr  sum_i C~_i f~_i by provided C~ satisfy
Cr    sum_j C_j f_j = sum_i C~_i f~_i = sum_ij C~_i L^-1_ij f_j
Cr                  = sum_j (sum_i L+^-1_ji C~_i) f_j
Cr  Writing as vectors
Cr    C+ f> = C+ L L^-1 f> = (L+ C)+ (L^-1 f>)
Cr  Therefore
Cr    C~ = L+ C  and  f~> = (L^-1 f>) and then C+ f> = (C~)+ f~>
Cr
Cr  If S_ij is the overlap in the old basis, in the new basis it is
Cr    < f~_i f~_j> = sum_km L^-1_ik S_km L^-1_mj
Cr             S~  = L^-1 S L^-1+
Cr             S   = L S~ L+
Cr
Cr *Choice of L to orthormal basis.  Then S~ = 1 and
Cr    S = L L+
Cr
Cr  1. No local orbital; basis consists of orthogonal (phi,phidot)
Cr     See above for overlap matrix.
Cr
Cr          (sqrt(s00)  0          )          (1/sqrt(s00)  0          )
Cr     L  = (                      )   L^-1 = (                        )
Cr          ( 0          sqrt(s11) )          ( 0          1/sqrt(s11) )
Cr
Cr  2. Local orbital; basis consists of (phi,phidot,gz)
Cr     See above for overlap matrix.
Cr
Cr     Let D = sqrt(szz - s0z^2/s00 - s1z^2/s11)
Cr
Cr             (sqrt(s00)          0            0)
Cr             (                                 )
Cr      L    = ( 0              sqrt(s11)       0)
Cr             (                                 )
Cr             (s0z/sqrt(s00)  s1z/sqrt(s11)    D)
Cr
Cr
Cr             (1/sqrt(s00)         0          0 )
Cr             (                                 )
Cr      L^-1 = ( 0              1/sqrt(s11)    0 )
Cr             (                                 )
Cr             (-s0z/s00/D     -s1z/s11/D     1/D)
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu    5 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   25 Apr 02 Returns cphi as coefficients to original nonlocal
Cu             orbitals, and cphin as matrix elements of evecs.
Cu             Altered argument list.
Cu   19 Feb 02 Extended to local orbitals.
Cu   28 Mar 01 written by MvS.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nsp,nlmax,ndham,nphimx,nbas,nev,lmxax,ndima,ipb(nbas),
     .  nlindx(3,0:lmxax,nbas)
      integer n0,nppn
      parameter (n0=10,nppn=12)
      double precision ppnl(nppn,n0,nsp,*),cphin(2,nev)
      double complex aus(nlmax,ndham,nphimx,nsp,*),cphi(ndima,nev)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      double precision lmat(3,3),pnu(n0,2),pnz(n0,2)
      integer lmxa,ichan,ib,is,iv,ilm,l,im,k,ia
      double precision s00,s11,szz,s0z,s1z,D
      double complex au,as,az,sqrsz(3)
C#ifdefC DEBUG
C      integer ndiml,offp,offd,offl,ilml,indx,i,iperm(3*nlmax),
C     .  nindx(0:n0),fopng
C      double precision, allocatable:: ovl(:,:)
C      double complex, allocatable:: cphil(:,:)
C#endif

      call dpzero(lmat,9)
      call dpzero(cphin,2*nev)
C     ichan0 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        ia = ipb(ib)
        pnu = s_spec(is)%p
        pnz = s_spec(is)%pz
        lmxa = s_spec(is)%lmxa
        if (lmxa == -1) cycle

C   ... Overlap matrix for this site
C#ifdefC DEBUG
C        ndiml = 2*(lmxa+1)**2  ! assumes phi,phidot in every channel
C        offp = 0
C        offd = (lmxa+1)**2
C        offl = 2*(lmxa+1)**2
C        do  l = 0, lmxa
C          nindx(l) = 2
C          if (pnz(l+1,1) /= 0) then
C            ndiml = ndiml + (2*l+1)
C            nindx(l) = 3
C          endif
C        enddo
C        allocate(ovl(ndiml,ndiml),cphil(ndiml,nev))
C        call dpzero(ovl,ndiml**2)
C        ilm = 0
C        do  l = 0, lmxa
C          s00 = ppnl(2,l+1,isp,ib)
C          s11 = ppnl(7,l+1,isp,ib)
C          do  im = 1, 2*l+1
C            ilm = ilm+1
C            ovl(ilm+offp,ilm+offp) = s00
C            ovl(ilm+offd,ilm+offd) = s11
C            cphil(ilm+offp,1:nev) = aus(ilm,1:nev,1,isp,ib)
C            cphil(ilm+offd,1:nev) = aus(ilm,1:nev,2,isp,ib)
C          enddo
C        enddo
C        ilm = 0
C        ilml = 0
C        indx = 0
C        do  l = 0, lmxa
C          do  im = 1, 2*l+1
C            ilm = ilm+1
C
CC           Permutation of indices to spex style
CC           This generates permutation to cphi actually passed
CC           indx should be initialized before 1st atom only
CC            do i = 1, nindx(l)
CC              indx = indx+1
CC              ichan = nlindx(i,l,ia) + im
CC              iperm(indx) = ichan
CC            enddo
C
CC           This is permutation to cphil, ovl for phi,dot
C            indx = indx+1
C            iperm(indx) = ilm+offp
C            indx = indx+1
C            iperm(indx) = ilm+offd
C
C            if (pnz(l+1,1) /= 0) then
C              ilml = ilml+1
C
C              indx = indx+1
C              iperm(indx) = ilml+offl
C
C              szz = ppnl(8,l+1,isp,ib)
C              s0z = ppnl(9,l+1,isp,ib)
C              s1z = ppnl(10,l+1,isp,ib)
C              ovl(ilml+offl,ilml+offl) = szz
C              cphil(ilml+offl,1:nev) = aus(ilm,1:nev,3,isp,ib)
C              ovl(ilm+offp,ilml+offl) = s0z
C              ovl(ilml+offl,ilm+offp) = s0z
C              ovl(ilm+offd,ilml+offl) = s1z
C              ovl(ilml+offl,ilm+offd) = s1z
C            endif
C            enddo
C        enddo
C
C        print *, 'writing file perm ...'
C        i = fopng('perm',-1,0)
C        write(i,"(i3)") iperm(1:indx)
C        call fclose(i)
C        call prmx('ovl',ovl,ndiml,ndiml,ndiml)
C        call zprm('cphi',2,cphil,ndiml,ndiml,nev)
C        deallocate(ovl,cphil)
C#endif

        do  iv = 1, nev
C         ichan = ichan0
          ilm = 0
          do  l = 0, lmxa

          k = l+1

          s00 = ppnl(2,k,isp,ib)
          s11 = ppnl(7,k,isp,ib)
          lmat(1,1) = sqrt(s00)
          lmat(2,2) = sqrt(s11)
          if (pnz(k,1) /= 0) then
            szz = ppnl(8,k,isp,ib)
            s0z = ppnl(9,k,isp,ib)
            s1z = ppnl(10,k,isp,ib)
            D = sqrt(szz - s0z**2/s00 - s1z**2/s11)
            lmat(3,1) = s0z/sqrt(s00)
            lmat(3,2) = s1z/sqrt(s11)
            lmat(3,3) = D
          else
            lmat(3,1) = 0
            lmat(3,2) = 0
            lmat(3,3) = 0
          endif

          if (.not. (pnz(k,1) == 0 .eqv. nlindx(3,l,ia) < 0)) then
             call rx('gwcphi: nlindx mismatch')
           endif

          az = 0
          do  im = 1, 2*l+1
            ilm = ilm+1
            au = aus(ilm,iv,1,isp,ib)
            as = aus(ilm,iv,2,isp,ib)
            if (nphimx > 2) az = aus(ilm,iv,3,isp,ib)

            sqrsz(1) = lmat(1,1)*au + lmat(3,1)*az
            sqrsz(2) = lmat(2,2)*as + lmat(3,2)*az
            sqrsz(3) =                lmat(3,3)*az
            cphin(1,iv) = cphin(1,iv) +
     .                    dconjg(sqrsz(1))*sqrsz(1) +
     .                    dconjg(sqrsz(2))*sqrsz(2) +
     .                    dconjg(sqrsz(3))*sqrsz(3)
            cphin(2,iv) = cphin(2,iv) +
     .                    dconjg(sqrsz(1))*sqrsz(1)



C            cphin(2,iv) = cphin(2,iv) +
C     .                    dconjg(lmat(1,1)*au)*lmat(1,1)*au

C           ichan = ichan+1
            ichan = nlindx(1,l,ia) + im
            cphi(ichan,iv) = au
            ichan = nlindx(2,l,ia) + im
            cphi(ichan,iv) = as
            if (nlindx(3,l,ia) >= 0) then
              ichan = nlindx(3,l,ia) + im
              cphi(ichan,iv) = az
            endif

          enddo
          enddo
        enddo
C       ichan0 = ichan

C#ifdefC DEBUG
C        call prmx('cphin',cphin(1,:),nev,nev,1)
CC       cphin = 0
C#endif

      enddo



      end
