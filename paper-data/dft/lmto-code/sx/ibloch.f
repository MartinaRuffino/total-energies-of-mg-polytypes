      subroutine ibloch(lcmplx,isp,nsp,qp,wtqp,ikp,plat,offH,indxsh,is1,
     .  is2,iax,ldima,ldimb,sll,s,nds)
C- Accumulate unsymmetrized inverse Bloch transform from one q
C ----------------------------------------------------------------------
Ci Inputs
Ci   lcmplx : 0 if r.s transform is real, 1 if complex
Ci   nsp,isp: 2 if s is spin-split, in which case isp=1 is spin index.
Ci   qp,wtqp,ikp: table of k-points (units of 2*pi/a) and weights,
Ci           (sum wtqp=2), and which kp sll corresponds to this call.
Ci   plat,iax: primitive lattice vectors and neighbor table; see strscr
Ci   offH :   Table of orbitals and ham. offsets info; see makidx
Ci   ldima,ldimb:  dimensions of sll
Ci   is1,is2: sum sll into pairs between is1,is2
Ci   sll :    Accumulate inverse transform of sll for this qp.
Co Outputs
Co   s :      The inverse Bloch transform is accumulated here.
Cb Bugs
Cb   Should be replaced with newer version
Cr Remarks:
Cr   ibloch produces an unsymmetrized transform, the inverse of the one
Cr   defined in subroutine bloch:
Cr     s(r1,l1,T+r2,l2) = sum_k w_k s(k;r1,l1,r2,l2) exp(-i k . T)
Cr   It is correct only if the k-sum is taken over the full Brillouin
Cr   zone.  If the irreducible BZ is used, the unsymmetrized transform
Cr   must be symmetrized; see rsmsym.f.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nkap0,n0H,niax
      parameter (nkap0=4,n0H=5,niax=10)
      integer lcmplx,nds,isp,nsp,ikp,is1,is2,ldima,ldimb,
     .  offh(n0H,nkap0,1),indxsh(1),iax(niax,is2)
      double precision sll(ldima,ldimb,2),s(nds,nds,is2),
     .  qp(3,ikp),wtqp(ikp),plat(3,3)
C Local parameters
      integer ia,ib,isite,j,k,i,off1,ksite,n0,
     .  lmb,ila,la,ma,offa,norba,
     .  lma,ilb,lb,mb,offb,norbb
      parameter (n0=10)
      integer ltaba(n0*nkap0),ltabb(n0*nkap0),
     .        ktab(n0*nkap0),offl(n0*nkap0)

      double precision twopi,TdotK,cosT,sinT,wc,ws
C ... debugging
*     integer is,ia0,ib0
C ... Convert offset for array starting at 1 to complex equivalent
      off1(lcmplx,i) = (lcmplx+1)*(i-1) + 1

*     is = 20

      twopi = 8*datan(1d0)
      if (nsp == 2) call rx('ibloch: nsp=2 not implemented')

C --- Add contribution to inverse Bloch sum to each R,R' pair ---
      do  100  isite = is1, is2

C   ... (wc,ws) = wtqp*exp(-i k . T)
        TdotK = 0
        do  30  j = 1, 3
        do  30  k = 1, 3
   30   TdotK = TdotK + twopi*qp(j,ikp)*plat(j,k)*iax(2+k,isite)
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
        wc  =  wtqp(ikp)/2*cosT
        ws  = -wtqp(ikp)/2*sinT

C   ... Get offsets and table of orbitals
        ia = iax(2,isite)
        ib = iax(1,isite)
C       uses norba,ltaba,norbb,ltabb
        call orbl(ia,0,ldima,indxsh,norba,ltaba,ktab,offa,offl,i)
        call orbl(ib,0,ldimb,indxsh,norbb,ltabb,ktab,offb,offl,i)
        ksite = off1(lcmplx,isite)
*       if (isite == is) ia0 = offH(1,ia)
*       if (isite == is) ib0 = offH(1,ib)

C   --- Accumulate inverse transform ---
        offb = offH(1,1,ib)
        do  40  ilb = 1, norbb
        lb = ltabb(ilb)
        lmb = lb*lb
        do  40  mb = -lb, lb
        offb = offb+1
        lmb = lmb+1

C   ... Real transform ...
        if (lcmplx == 0) then
          offa = offH(1,1,ia)
          do  42  ila = 1, norba
          la = ltaba(ila)
          lma = la*la
          do  42  ma = -la, la
          offa = offa+1
          lma = lma+1
            s(lma,lmb,ksite) = s(lma,lmb,ksite) +
     .      sll(offa,offb,1)*wc - sll(offa,offb,2)*ws
   42     continue
C   ... Complex transform (loops repeated to optimize speed) ...
        else
          offa = offH(1,1,ia)
          do  44  ila = 1, norba
          la = ltaba(ila)
          lma = la*la
          do  44  ma = -la, la
          offa = offa+1
          lma = lma+1
            s(lma,lmb,ksite) = s(lma,lmb,ksite) +
     .      sll(offa,offb,1)*wc - sll(offa,offb,2)*ws
            s(lma,lmb,ksite+1) = s(lma,lmb,ksite+1) +
     .      sll(offa,offb,2)*wc + sll(offa,offb,1)*ws
   44     continue
        endif
   40   continue

*        if (is == isite) then
*         print 335, isite, sll(ia0+1,ib0+1,1),sll(ia0+1,ib0+1,2),
*     .    s(1,1,ksite),s(1,1,ksite),cosT,sinT
*        endif
*  335 format('i',i4,6f10.5)

  100 continue

*     print 335, off1(lcmplx,20),s(1,1,off1(lcmplx,20))
      end

C      subroutine snot(plat,nl,nbas,lmx,nbasp,ipl,npl,iclspp,pgplp,indxsh
C     .  ,nstar,qstar,rmat,ldim,ldim0,sk,s_ibl)
CC- Test ibloch.f
C      implicit none
C      integer nl,nbas,nbasp,ldim,ldim0,npl,ipl,nstar,
C     .  lmx(*),iclspp(nbasp),indxsh(1),pgplp(6,0:4),s_ibl(4)
C      double precision plat(3,3),qstar(3,nstar),sk(ldim0,ldim,2),
C     .  pti(2,ldim0),rmat(nl,nl,nl,nl,1)
CC local variables
C      integer nitab,is1,is2,nbaspp,nbaspl,
C     .  oiclpl,oiax,ontab,os,oioffH,oalph,jpl,ib1,isum,ib2,osr
C      logical ndnfld,ccor,spncpl,lerr
C      double precision ekap(2)
CC heap:
C      integer w(1)
C      common /w/ w
C
CC --- Setup ---
CC     call wkprnt(1)
C      nbaspp = 2*nbasp - nbas
C
CC --- Read all real-space strx from disk, 2D Bloch transform ---
C      call iostrx(8,'STR',nl,nbasp,1,ekap,itral,ckbas,lmaxw,nitab,oalph,
C     .  oiax,ontab,os)
C      call defi(oioffH, nbaspp)
C      call defi(oiclpl,-nbaspp)
C      call pghof(ipl,npl,pgplp,nbas,lmx,iclspp,w(ontab),w(oiax),
C     .    nl,ldim,ldim0,is1,is2,nbaspl,w(oiclpl),w(oioffH),indxsh)
C      call blchpl(qstar,nl**2,plat,lmx,iclspp,indxsh,w(oiax),w(os),
C     .  ldim,ldim0,is1,is2,w(oioffH),.false.,sk,sk)
CC --- Inverse Bloch transform ---
C      jpl = max(min(ipl,npl-1),0)
C      stop 'pgplp'
C      ib1 = isum(jpl,pgplp(1,0),6)+1
C      ib2 = ib1-1 + pgplp(1,jpl)
C      ontab = s_ibl(2)
C      oiax  = s_ibl(3)
C      osr   = s_ibl(4)
C      is1 = s_str%npr(ib1)+1
C      is2 = s_str%npr(ib2+1)
C      call ibloch(nstar,qstar,rmat,2,nl**2,plat,lmx,iclspp,w(oiax)
C     .  ,is1,is2,ldim,ldim0,w(oioffH),sk,w(osr),lerr)
C      if (lerr) call rx('problem with ibloch')
C      call rlse(oalph)
C      end
