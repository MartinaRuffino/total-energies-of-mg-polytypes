      subroutine percos(nhs,isp,nsp,nstate,h,z,evl)
C- Undo hamiltonian fix for a given spin by pertubation
C      implicit none
      integer nhs,nstate,isp,nsp,nstx
      double precision h(nhs,nhs),z(nhs,nhs),evl(nhs,nsp),xx

C For MPI ...
      integer mpipid,procid,master,numprocs
      logical MPI,mlog,cmdopt
      character*80 outs
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      nstx = min0(nstate+5,nhs)
      if (procid == master) then
        rewind 94
        read (94) nhs
        call dpdump(h,nhs*nhs,94)
      endif
      call mpibc1(nhs,1,2,mlog,'percos','nhs')
      call mpibc1(h,nhs*nhs,4,mlog,'percos','h')
C ... nsp=2: get evecs appropriate to this isp
      if (nsp == 2) then
        if (procid == master) then
          if (isp == 1) rewind 95
          read(95) xx
          call dpdump(z,nhs*nhs,95)
        endif
        call mpibc1(z,nhs*nhs,4,mlog,'percos','z')
      endif
      call perco1(nhs,nstx,nstate,z,h,evl(1,isp))
c|    call dpdump(h,nhs*nhs,94)
c|    call perco2(nhs,nstx,nstate,z,h,vint,vintf,w(oevl))
      end
