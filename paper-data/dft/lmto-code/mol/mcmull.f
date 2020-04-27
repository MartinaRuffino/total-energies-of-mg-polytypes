      subroutine mcmull(nsp,nhs,nbas,atid,z,n0,pnu,pnz,lmxa,ips,nel,
     .                  ioffb,ewt,t,s,c,pos,qmull,qc)
C- Simple atom resolved Mulliken charges
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nsp, nhs, nbas, atid, nel,ioffb,
Ci   z,n0,pnu,pnuz,lmxa
Ci   ewt: occupancies
Ci   t: eigenvectors
Ci   s: overlap
Ci   c: n x n work array (probably not really needed)
Co Outputs:
Co   qmull: Mulliken charges
Co   qc:    core and nuclear charges
Cr Remarks
Cu Updates
Cu   22 Jan 2007 (ATP)
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer n0,nsp,nhs,nbas,nel
      integer lmxa(1),ips(nbas),ioffb(nbas,nel)
      double precision z(1),c(nhs,nhs),t(nhs,nhs),s(nhs,nhs),
     .                 pnu(n0,2,nbas),pnz(n0,2,nbas),
     .                 pos(3,nbas),qmull(nbas,nsp),qc(nbas,2),ewt(nhs,1)
      character*8 atid(nbas)
C Local Variables
      integer lmax,isp,n,i1,i2,ie,ib,ilmi,nlmi,ncore
      integer i,iprint
      double precision wt,mmom,dsum,qcor(2),z1,dip(3),q,qtot
      double precision Y0,D,Cm,pi,datan,dsqrt
      data D /2.541748/ Cm /8.478358/

C For MPI ...
      integer mpipid,procid,master,numprocs
      logical MPI,mlog,cmdopt
      character*80 outs
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      pi = 4d0*datan(1d0)
      Y0 = 1/dsqrt(4d0*pi)
      call dcopy(nbas*nsp,0d0,0,qmull,1)

      do  isp = 1, nsp

        if (nsp == 2) then
C --- get eigenvectors from disc ---
          if (procid == master) then
            if (isp == 1) rewind 95
            read (95) wt
            call dpdump(t,nhs*nhs,95)
          endif
          call mpibc1(wt,1,4,mlog,'mcmull','wt')
          call mpibc1(t,nhs*nhs,4,mlog,'mcmull','t')
        endif

C --- make c = S t ---
        n = nhs
        call dmpy(s,n,1,t,n,1,c,n,1,n,n,n)

C --- make q_R = \sum_n(occ.) \sum_Le  t_RLe c_RLe ---
        do ie = 1, nel
          do ib = 1, nbas
            nlmi=ioffb(ib+1,ie)-ioffb(ib,ie)
            i1 = ioffb(ib,ie) + 1
            i2 = i1 + nlmi - 1
            if (i2 > nhs) call rx('Bug in mcmull. i2>nhs')
            do  n = 1, nhs
              wt = ewt(n,isp) * 2d0/nsp
              do  ilmi = i1, i2
                qmull(ib,isp) = qmull(ib,isp)
     .                        + t(ilmi,n) * c(ilmi,n) * wt
              enddo
            enddo
          enddo
        enddo

      enddo

C --- get core charges ---
      do  ib = 1, nbas
        z1 = z(ips(ib))
        lmax = lmxa(ips(ib))
        call pshpr(1)
        call getcor(2,z1,0d0,pnu(1,1,ib),pnz(1,1,ib),0d0,lmax,0d0,0d0,0,
     .              0,qcor,0d0,0d0,0d0,ncore,0d0,0d0)
        call poppr
        qc(ib,1) = qcor(1)
        qc(ib,2) = -z1
      enddo

C --- printout ---
      if (iprint() < 30) return
      if (nsp == 1) then
        write (*,1)
      else
        write (*,2)
      endif
      do  ib = 1, nbas
        if (nsp == 1) then
          write (*,3) ib,atid(ib),qmull(ib,1),qc(ib,1),qc(ib,2),
     .                qmull(ib,1)+qc(ib,1)+qc(ib,2)
        else
          write (*,3) ib,atid(ib),qmull(ib,1),qmull(ib,2),
     .                qc(ib,1),qc(ib,2),
     .                qmull(ib,1)+qmull(ib,2)+qc(ib,1)+qc(ib,2)
        endif
      enddo
      call dcopy(3,0d0,0,dip,1)
      qtot = 0d0
      do  ib = 1, nbas
        q = dsum(nsp,qmull(ib,1),nbas) + dsum(2,qc(ib,1),nbas)
        qtot = qtot + q
        do  i = 1, 3
          dip(i) = dip(i) + q * pos(i,ib)
        enddo
      enddo
      if (nsp == 1) then
        write (*,4) dsum(nbas,qmull,1),qtot
      else
        mmom = dsum(nbas,qmull(1,1),1) - dsum(nbas,qmull(1,2),1)
        write (*,5) dsum(nbas*2,qmull,1),mmom
      endif
      write (*,6) qtot*Y0,(Y0*dip(i),i=1,3),-qtot,(-dip(i),i=1,3),
     .  (-D*dip(i),i=1,3),(-Cm*dip(i),i=1,3)

    1 format (/' MCMULL: Mulliken charges (e=1) ...'/
     .         '  atom  species',8x,'valence',4x,'(core)  (nucleus)',
     .         '   total')
    2          format (/' MCMULL: Mulliken charges ...'/
     .         '  atom  species',12x,'valence',6x,'(core)  (nucleus)',
     .         '   total'/
     .         21x,'spin up       down ')
    3          format (2x,i3,5x,a8,3x,5f10.6)
    4 format ('       Total charge: ', f10.6,20x,f10.6)
    5 format ('       Total charge: ', f10.6,' moment: ',f10.6)
    6 format(/'    Total dipole moment:'/
     .       22x,'q',10x,'px',10x,'py',10x,'pz'
     .       /14x,4f12.6/'    (e=-1)',4x,4f12.6/
     .         '    Debye (10^-18 esu)',4x,3f12.6/
     .         '    10^-30 Cm',13x,3f12.6/)

      end

