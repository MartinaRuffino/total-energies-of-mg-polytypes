      subroutine asxsig(aintra,iqin,plat,nl,nsp,nbas,offH,indxsh,
     .  nkp,wtkp,lshft,ipq,igstar,ldim,lpdim,eval,nbmx,evec,bas,istab,g,
     .  ag,efermi,wscrk,w0,w1,w2,z,h3,wtr,sigma)
C- ASA screened exchange potential from nonlocal, static linear response
C ----------------------------------------------------------------------
Ci Inputs
Ci   eval,evec
Ci   w0,w1,w2  W for q-> 0
Ci   wscrk     W(q)
Ci   wtr       work array to hold W(rotated q)
Ci   z         work array to hold evec(rotated q)
Co Outputs
Co   sigma
Cr Remarks
Cr   Given q-dependent dW_ab = diff between screened coulomb potential
Cr   generated from nonlocal P and local P,
Cr
Cr      X           i=occ     i     i
Cr     V_ab(q) = sum(k)      z (k) z (k) dW (k-q)
Cr                            a     b      ab
Cr
Cr   clean up complex arithmetic; replace indxsh with offH in wtrans
Cr   For version 1.  See asxsgm for v2.
C ----------------------------------------------------------------------
      implicit none
      logical lshft(3)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer iqin,nbmx,nbas,nl,nsp,ldim,nkp,lpdim,istab(nbas,*),
     .  offH(n0H,nkap0,*),indxsh(*),ipq(*),igstar(0:*)
      double precision plat(3,3),bas(3,*),eval(nbmx,nkp),efermi,
     .  g(3,3,*),ag(3,*),evec(ldim,ldim,2,nkp),wtkp(1),z(ldim,ldim,2),
     .  h3(ldim,ldim,2),w0(lpdim,lpdim),w1(lpdim,lpdim),w2(lpdim,lpdim),
     .  wtr(lpdim,lpdim,2),wscrk(lpdim,lpdim,2,nkp),sigma(ldim,ldim,2)
      logical aintra
C Local variables
      integer iq0,j1,j2,j3,iq1,iq2,i1,i2,i3,ii1,ii2,ii3,jpq0,ipq1,
     .  jpq1,ipq2,jpq2,nu1,iprint,
     .  k,getdig,nkxyz(3),is(3),ifac(3),jj1,jj2,jj3,nlx,nuocc
      parameter (nlx=4)
      double precision wght1,rb(3,3),qb(3,3),qlat(3,3),q1(3),q2(3),qk,
     .  rmat(nlx**4),pi,q0m,volbz,wght2,wght3
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C --- Setup ---
      call rxx(nl > nlx,'asxsig: increase nlx')
      call rxx(nsp == 2,'asxsig not set up for spin pol')
      k = iabs(igstar(0))
      call cmpi123(0,nkxyz(1),nkxyz(2),nkxyz(3),k)
C ... Make volbz
      call dinv33(plat,1,qlat,wght1)
      volbz = 1/wght1
C ... Make is,ifac,qb,wght1
      call pshpr(1)
      call bzmsh0(plat,lshft,0,nkxyz(1),nkxyz(2),nkxyz(3),is,ifac,rb,qb)
      call poppr
      wght1 = 2d0/nkxyz(3)/nkxyz(2)/nkxyz(1)
      pi = 4d0*datan(1d0)
      q0m = (volbz/4d0*3d0/pi/nkxyz(3)/nkxyz(2)/nkxyz(1))**(1d0/3d0)
      wght2 = wght1*3d0/2d0/q0m
      wght3 = wght1*3d0/q0m**2
      call dpzero(sigma,ldim**2*2)

c ---- identify iqin ----------
      iq0 = 0
      do  10  j3 = 1, nkxyz(3)
      do  10  j2 = 1, nkxyz(2)
      do  10  j1 = 1, nkxyz(1)
        iq0 = iq0+1
        if (iqin /= ipq(iq0)) goto 10
        jpq0 = igstar(iq0)
        if (jpq0 /= 1) goto 10

C --- k-integration for vsxnl(q) ---
        iq1 = 0
        do  20  i3 = 1, nkxyz(3)
        ii3 = i3+j3-1
        if (ii3 > nkxyz(3)) ii3 = ii3-nkxyz(3)
        do  20  i2 = 1, nkxyz(2)
        ii2 = i2+j2-1
        if (ii2 > nkxyz(2)) ii2 = ii2-nkxyz(2)
        do  20  i1 = 1, nkxyz(1)
        ii1 = i1+j1-1
        if (ii1 > nkxyz(1)) ii1 = ii1-nkxyz(1)

C   --- Potential for q1=k ---
        iq1 = iq1+1
        ipq1 = ipq(iq1)
        jpq1 = igstar(iq1)
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)
C       call yprm('w',2,wscrk(1,1,1,ipq1),lpdim*lpdim,lpdim,lpdim,lpdim)
        call wtrans(ldim,lpdim,nl,nbas,bas,indxsh,istab(1,jpq1),
     .    g(1,1,jpq1),ag(1,jpq1),q1,aintra,wscrk(1,1,1,ipq1),wtr)
C       call yprm('wtrans',2,wtr,lpdim*lpdim,lpdim,lpdim,lpdim)

C   --- Wave function for q2=k+q ---
        iq2 = ii1+nkxyz(1)*((ii2-1)+nkxyz(2)*(ii3-1))
        ipq2 = ipq(iq2)
        jpq2 = igstar(iq2)
        q2(1) = qk(1,ii1,ii2,ii3)
        q2(2) = qk(2,ii1,ii2,ii3)
        q2(3) = qk(3,ii1,ii2,ii3)
        call rotwf(0,nl,nbas,1,bas,offH,indxsh,istab(1,jpq2),g(1,1,jpq2),
     .    ag(1,jpq2),q2,rmat,0,ldim,ldim,ldim,evec(1,1,1,ipq2),z)

C   --- Count number of occupied bands ---
        do  26  nu1 = 1, ldim
        nuocc = nu1-1
   26   if(eval(nu1,ipq1) > efermi) goto 27
   27   continue
C   --- h3 = sum_nu=occ z_RL1,nu (z+)_RL2,nu ---
        call yygemm('N4','C',ldim,ldim,nuocc,1d0,
     .    z,z(1,1,2),ldim,z,z(1,1,2),ldim,0d0,h3,h3(1,1,2),ldim)

C   --- Accumulate sum_nu=occ z_RL1,nu (z+)_RL2,nu W_Rl1,Rl2 into Sigma
        call psxsig(aintra,nbas,indxsh,ldim,lpdim,w0,w1,w2,
     .    iq1,wght1,wght2,wght3,h3,wtr,sigma)

   20 continue
   10 continue

      if (iprint() >= 90) then
C        write(*,*) ' --- Re(sigma) -----'
C        write(*,123) ((sigma(i,j,1),j=1,10),i=1,10)
C        write(*,*) ' --- Im(sigma) -----'
C        write(*,123) ((sigma(i,j,2),j=1,10),i=1,10)
        call yprm('sigma',2,sigma,ldim*ldim,ldim,ldim,ldim)
      endif
C 123 format(10f8.3)
      end

      subroutine psxsig(aintra,nbas,indxsh,ldim,lpdim,w0,w1,w2,
     .  iq1,wght1,wght2,wght3,h3,wtr,sigma)
C- Kernel called by asxsig
      implicit none
      integer nbas,ldim,lpdim,indxsh(*),iq1
      double precision h3(ldim,ldim,2),wtr(lpdim,lpdim,2),
     .  w0(lpdim,lpdim),w1(lpdim,lpdim),w2(lpdim,lpdim),
     .  sigma(ldim,ldim,2),wght1,wght2,wght3
      logical aintra
C Local variables
      integer mxorb,nglob,lmxa,
     .  nl1,ib1,off1,lmr1,nlm1,l1,m1,
     .  nl2,ib2,off2,lmr2,nlm2,l2,m2
      double precision xxr,xxi

      mxorb = nglob('mxorb')
      lmxa = nglob('nl')-1

C --- Loop over orbitals ib2, extracting Rl(ib2) for wtr ---
      nl2 = 0
      do  40  ib2 = 1, nbas
        lmr2 = mxorb*(ib2-1)
        do  42  l2 = 0, lmxa
          off2 = indxsh(lmr2+1) - 1
          lmr2 = lmr2 + 2*l2+1
          if (off2 >= ldim) goto 42
          nlm2 = 2*l2+1
          nl2 = nl2+1
C       --- Loop over orbitals ib1, extracting Rl(ib1) for wtr ---
            nl1 = 0
            do  50  ib1 = 1, nbas
              lmr1 = mxorb*(ib1-1)
              do  52  l1 = 0, lmxa
                off1 = indxsh(lmr1+1) - 1
                lmr1 = lmr1 + 2*l1+1
                if (off1 >= ldim) goto 52
                nlm1 = 2*l1+1
                nl1 = nl1+1
                if (iq1 == 1) then
                  if (aintra) then
                    xxr = w0(nl1,nl2)*wght1 + w2(nl1,nl2)*wght3
                    xxi = w1(nl1,nl2)*wght2
                  else
                    xxr = w0(ib1,ib2)*wght1 + w2(ib1,ib2)*wght3
                    xxi = w1(ib1,ib2)*wght2
                  endif
                else
                  if (aintra) then
                    xxr = wtr(nl1,nl2,1)*wght1
                    xxi = wtr(nl1,nl2,2)*wght1
                  else
                    xxr = wtr(ib1,ib2,1)*wght1
                    xxi = wtr(ib1,ib2,2)*wght1
                  endif
                endif
C           ... Accumulate z1_Rlm1 W_Rl1,Rl2 (z+)_Rlm2  into Sigma
                do  58  m2 = 1, nlm2
                do  58  m1 = 1, nlm1
                  sigma(m1+off1,m2+off2,1) = sigma(m1+off1,m2+off2,1) +
     .             h3(m1+off1,m2+off2,1)*xxr - h3(m1+off1,m2+off2,2)*xxi
                  sigma(m1+off1,m2+off2,2) = sigma(m1+off1,m2+off2,2) +
     .             h3(m1+off1,m2+off2,2)*xxr + h3(m1+off1,m2+off2,1)*xxi
   58           continue
   52         continue
   50       continue
   42   continue
   40 continue
        end
