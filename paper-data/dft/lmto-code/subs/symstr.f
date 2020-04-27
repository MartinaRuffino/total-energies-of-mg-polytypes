      subroutine symstr(mode,nds,nsites,iax,nsp,nkap,sflg,s,sc,asym)
C- Symmetrize structure constants
C ----------------------------------------------------------------
Ci Inputs
Ci   mode  :0 s is real
Ci         :1 s is complex
Ci         :2 s is complex; symmetrize s+_ij = s_ji
Ci   nds   :leading dimension of s
Ci   nsites:number of pairs
Ci   iax   :neighbor table containing pair information
Ci         :iax(*,i) contains information for pair i
Ci         :This routine uses:
Ci         :iax(6,i): index to conjugate (jb,ib) pair matching (ib,jb)
Ci         :iax(8,i): if nonzero, points to an pair which is to
Ci         :        : substitute for pair i.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cio Inputs/Outputs
Cio  sflg  :marks which pairs in s have been symmetrized.
Cio        :On output, any pairs symmetrized by symstr are marked.
Cio  s,sc  :real-space hamiltonian or structure-constant matrix.
Cio        :Only s is used if mode=0 or sc is used if mode=1.
Cio        :Any pair i which has a counterpart j is symmetrized.
Cio        :*If neither pair i not pair j=iax(6,i) is symmetrized,
Cio        : both j and i are symmetrized by averaging the two.
Cio        :*If j=iax(6,i) is already symmetrized, i is copied from j.
Cio        :*If i is already symmetrized, j is copied from i.
Co Outputs
Co   asym  :maximum asymmetry found when symmetrizing
Cr Remarks
Cu Updates
Cu   05 Aug 06 Passes parms for 2-kappa (just return for now)
Cu   30 Mar 03 complex case (mode>0) can symm. s_ij=s_ji or s+_ij=sji
Cu   17 Jul 02 Redesigned for more general cases (complex, nsp)
Cu             Altered argument list.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,nds,nsites,nsp,sflg(nsites),niax,nkap
      parameter (niax=10)
      integer iax(niax,nsites)
      double precision s(nds,nds,nsp,nkap,nkap,nsites),asym
      double complex sc(nds,nds,nsp,nkap,nkap,nsites)
C Local parameters
      integer i,j,lm1,lm2,isp,ii,jj,ik,jk
      double precision tmp
      double complex tmpc

C Cludge for now ...
      asym = 0
      if (nkap /= 1) return
      ik = 1
      jk = 1

      call sanrg(.true.,mode,0,2,'symstr:','mode')

      do  i = 1, nsites
        j = iax(6,i)
        if (j == 0) cycle

        ii = iax(8,i)
        if (ii == 0) ii = i
        jj = iax(8,j)
        if (jj == 0) jj = j

c       print *, i,ii,' ',j,jj,' ',sflg(ii),sflg(jj)
C   ... Symmetrize ii,jj as (ii+jj)/2 if neither was symmetrized
        if (sflg(ii) == 0 .and. sflg(jj) == 0) then
          if (mode == 0) then
            do  isp = 1, nsp
              do  lm1 = 1, nds
                do  lm2 = 1, nds
                  tmp = (s(lm1,lm2,isp,ik,jk,ii)
     .                  +s(lm2,lm1,isp,ik,jk,jj))/2
                  asym = max(asym,abs(s(lm1,lm2,isp,ik,jk,ii)-tmp))
                  asym = max(asym,abs(s(lm2,lm1,isp,ik,jk,jj)-tmp))
C                  if (asym > 1d-3) then
C                    print *, 'asym=',asym
C                  endif
                  s(lm1,lm2,isp,ik,jk,ii) = tmp
                  s(lm2,lm1,isp,ik,jk,jj) = tmp
                enddo
              enddo
            enddo
          else if (mode == 1) then
            do  isp = 1, nsp
              do  lm1 = 1, nds
                do  lm2 = 1, nds
                  tmpc = (sc(lm1,lm2,isp,ik,jk,ii)
     .                   +sc(lm2,lm1,isp,ik,jk,jj))/2
                  asym = max(asym,abs(sc(lm1,lm2,isp,ik,jk,ii)-tmpc))
                  asym = max(asym,abs(sc(lm2,lm1,isp,ik,jk,jj)-tmpc))
C                  if (asym > .004) then
C                    print *, i,lm1,lm2
C                  endif
                  sc(lm1,lm2,isp,ik,jk,ii) = tmpc
                  sc(lm2,lm1,isp,ik,jk,jj) = tmpc
                enddo
              enddo
            enddo
          else
            do  isp = 1, nsp
            do  lm1 = 1, nds
            do  lm2 = 1, nds
            tmpc = (sc(lm1,lm2,isp,ik,jk,ii)+
     .       dconjg(sc(lm2,lm1,isp,ik,jk,jj)))/2
            asym = max(asym,abs(sc(lm1,lm2,isp,ik,jk,ii)-tmpc))
            asym =max(asym,abs(dconjg(sc(lm2,lm1,isp,ik,jk,jj))-tmpc))
            sc(lm1,lm2,isp,ik,jk,ii) = tmpc
            sc(lm2,lm1,isp,ik,jk,jj) = dconjg(tmpc)
            enddo
            enddo
            enddo
          endif
C     ... Flag this s(ii), s(jj) as symmetrized
          sflg(ii) = 1
          sflg(jj) = 1
C   ... Symmetrize ii from jj if ii not symmetrized
        elseif (sflg(ii) == 0) then
          do  isp = 1, nsp
            do  lm1 = 1, nds
              do  lm2 = 1, nds
                if (mode == 0) then
                  s(lm1,lm2,isp,ik,jk,ii) = s(lm2,lm1,isp,ik,jk,jj)
                else if (mode == 1) then
                  sc(lm1,lm2,isp,ik,jk,ii) = sc(lm2,lm1,isp,ik,jk,jj)
                else
                  sc(lm1,lm2,isp,ik,jk,ii) =
     .              dconjg(sc(lm2,lm1,isp,ik,jk,jj))
                endif
              enddo
            enddo
          enddo
C     ... Flag this s(ii) as symmetrized
          sflg(ii) = 1
C   ... Symmetrize jj from ii if jj not symmetrized
        elseif (sflg(jj) == 0) then
          do  isp = 1, nsp
            do  lm1 = 1, nds
              do  lm2 = 1, nds
                if (mode == 0) then
                  s(lm2,lm1,isp,ik,jk,jj) = s(lm1,lm2,isp,ik,jk,ii)
                else if (mode == 1) then
                  sc(lm2,lm1,isp,ik,jk,jj) = sc(lm1,lm2,isp,ik,jk,ii)
                else
                  sc(lm2,lm1,isp,ik,jk,jj) =
     .              dconjg(sc(lm1,lm2,isp,ik,jk,ii))
                endif
              enddo
            enddo
          enddo
C     ... Flag this s(jj) as symmetrized
          sflg(jj) = 1
        endif
      enddo

C     call info(60,0,0,' symstr: max asymmetry=%d',asym,0)
      end
