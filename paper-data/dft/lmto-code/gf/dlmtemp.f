      subroutine dlmtemp(mode,nclass,nthet,theta,wtg,tdlm,etrms,temp)
      implicit none
Ci Inputs
Ci   wtg :Gibbs weights including Euler quardature factor, per class
Ci       :wtg = 0 for DLM dummy classes
Cr Remarks
Cr    mode = 0 : use total energy from etrms(2,:)
Cr    mode = 1 : use sumev from etrms(16,:)
C ... Passed parameters
      integer mode,nclass,nthet(nclass)
      double precision theta(*),wtg(*),tdlm,etrms(22,*),temp(*)
C ... Local parameters
      integer ic,idc,idlm,ie
      double precision x,wt
      integer nmax
      parameter (nmax=4)
      double precision amat(nmax,nmax),bvec(nmax),iwk(nmax)
      integer lgunit,n1,n2,ith,info
      double precision pl1,pl2,plegn
      double precision boltz
      parameter (boltz=8.6173324d-5/13.605698d0)
      integer mpipid,procid,master
c     pl0(x) = 1d0
c     pl1(x) = x
c     pl2(x) = 0.5d0*(3d0*x*x-1d0)
c     pl3(x) = 0.5d0*(5d0*x*x*x-3d0*x)

      procid = mpipid(1)
      master = 0

      if (procid == master) write(lgunit(1),500)
      if (mode == 0) then
        ie = 2
        if (procid == master) write(lgunit(1),501)
      elseif (mode == 1) then
        ie = 16
        if (procid == master) write(lgunit(1),502)
      else
        call rx1('Invalid mode %i',mode)
      endif

      idc = 0
      idlm = 0
      do  ic = 1, nclass
        if (nthet(ic) < 2) cycle
        idlm = idlm + 1
        call dpzero(amat,nmax*nmax)
        call dpzero(bvec,nmax)
C ...   Coefficients for P_n[cos(th)]
        if (nthet(ic) == 2) then
          idc = idc + 1
          if (idc == 1 .and. procid == master) write(lgunit(1),503)
          bvec(1) = etrms(ie,nclass+idc)
          bvec(2) = bvec(1)
          idc = idc + 1
          bvec(1) = bvec(1) + etrms(ie,nclass+idc)
          bvec(2) = bvec(2) - etrms(ie,nclass+idc)
          bvec(1) = 0.5d0*bvec(1)
          bvec(2) = 0.5d0*bvec(2)
        else
          do  ith = 1,nthet(ic)
            idc = idc + 1
            if (idc == 1 .and. procid == master) write(lgunit(1),503)
            x = cos(theta(idc))
            wt = wtg(nclass+idc)
            do  n1 = 1, nmax
              pl1 = plegn(n1-1,x)
              bvec(n1) = bvec(n1) + wt * pl1 * etrms(ie,nclass+idc)
              do  n2 = 1, nmax
                pl2 = plegn(n2-1,x)
                amat(n1,n2) = amat(n1,n2) + wt * pl1 * pl2
              enddo
            enddo
          enddo
C         call yprm('amat',1,amat,0,nmax,nmax,nmax)
          call dgefa(amat,nmax,nmax,iwk,info)
          if (info /= 0) call rx('GETTEMP: Matrix singular')
          call dgesl(amat,nmax,nmax,iwk,bvec,0)
        endif
        temp(idlm) = bvec(2)/tdlm/boltz
        if (procid == master)
     .    write(lgunit(1),504) ic, (bvec(n1),n1=1,nmax),temp(idlm)
      enddo

 500  format(/' Legendre coefficients and temperatures for DLM classes')
 501  format(' (Using total energies)')
 502  format(' (Using single-electron energies)')
 503  format('   Class',8x,'E0',10x,'E1',8x,'E2',8x,'E3',9x,'T, K')
 504  format(3x,i3,1X,F15.6,3F10.6,F12.2)
      end
