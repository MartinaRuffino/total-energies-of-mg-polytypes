      subroutine rxcmro(dx,dy,dz,pos,nbas,rsm,rint,nxi,lxi,exi,n0,
     .   ips,ioff,nri,rhoi,cg,jcg,indxcg,cenx,ceny,cenz,z,k1,k2,l1,l2,
     .   rhot,grhot,g2rhot,ggrhot,ntot)
C- loops over all atoms, makes density and density derivatives
C  in one xy-plane
Cu Updates
Cu   21 Jul 07 (S. Lozovoi) New implementation: rho derivs for GGA

      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   cg(1),jcg(1),indxcg(1),
     .   rhoi(nri,*),pos(3,1),ipd(1),ioff(1),ips(1),
     .   rhot(k1:k2,l1:l2,*),
     .   grhot(k1:k2,l1:l2,3,*),
     .   g2rhot(k1:k2,l1:l2,5,*),
     .   ggrhot(k1:k2,l1:l2,*)
      real w(1)
      common /w/ w

C For MPI ...
      integer, dimension(:), allocatable :: bproc
      integer mpipid,lgunit,procid,pid,master,numprocs,length,ierr
      integer i1mach,iprint,ib1,ib2,len
      logical mlog,cmdopt,MPI
      character outs*80, shortname*10, datim*26
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      call mpiprn(shortname,length)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      call tcn('rxcmro')
      nsp = lsp()+1
      nsum = 0d0
      nat = 0
c --------- start loop over atoms ---------------
      if (MPI) then
        allocate (bproc(0:numprocs), stat=ierr)
        call dstrbp(nbas,numprocs,1,bproc(0))
        ib1 = bproc(procid)
        ib2 = bproc(procid+1)-1
      else
        ib1 = 1
        ib2 = nbas
      endif
      do  10  ib = ib1, ib2
        if (MPI .and. mlog) then
        if (ib == bproc(procid)) then
          call gettime(datim)
          call awrit4(' rxcmro '//datim//' Process %i of %i on '
     .        //shortname(1:length)//
     .        ' starting atoms %i to %i',' ',256,lgunit(3),
     .         procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        endif
        is = ips(ib)
        posx = pos(1,ib)
        posy = pos(2,ib)
        zz = z-pos(3,ib)
        ri = rint(is)
        npmx = 4d0*(ri*ri-zz*zz)/(dx*dy)
        npmx = max0(npmx,40)
c --------- gather x,y coordinates --------------
        call defrr(ox,    npmx)
        call defrr(oy,    npmx)
        call defrr(oz,    npmx)
        call rxcgth(posx,posy,cenx,ceny,dx,dy,ri,zz,np,
     .              w(ox),w(oy),w(oz))
        if (np == 0) goto 10
        nsum = nsum + np
        nat = nat + 1
        if (np > npmx) call rx('rxcmro: np gt npmx')
c --------- make and add density in xy-plane --------
        call defrr(orho,   npmx*nsp)
        if(lxcg() > 0) then
          call defrr(ogrho,  npmx*nsp*3)
          call defrr(og2rho, npmx*nsp*5)
          call defrr(oggrho, npmx*nsp)
        else
          ogrho   = 1
          og2rho  = 1
          oggrho  = 1
        endif

        ir = ioff(ib) + 1
        call rxcrh1(nxi(is),lxi(1,is),exi(1,is),rsm(is),rint(is),nri,
     .   rhoi(ir,1),cg,jcg,indxcg,np,w(ox),w(oy),w(oz),
     .   w(orho),w(ogrho),w(og2rho),w(oggrho),np)
        call rxcadd(posx,posy,cenx,ceny,dx,dy,ri,zz,np,w(orho),
     .    w(ogrho),w(og2rho),w(oggrho),
     .    k1,k2,l1,l2,rhot,grhot,g2rhot,ggrhot)
        call rlse(ox)
   10 continue
      if (MPI) then
        len = (k2 - k1 + 1)*(l2 - l1 + 1)*nsp
        if(lxcg() > 0) then
          call mpibc2(grhot,3*len,4,3,mlog,'rxcmro','grhot')
          call mpibc2(g2rhot,5*len,4,3,mlog,'rxcmro','grhot')
          call mpibc2(ggrhot,len,4,3,mlog,'rxcmro','grhot')
        endif
        call mpibc2(rhot,len,4,3,mlog,'rxcmro','rhot')
        deallocate (bproc, stat=ierr)
      endif

      ntot = ntot + nsum

      if(iprint() >= 50) write(6,817) z,k1,k2,l1,l2,nat,nsum
  817 format('   z=',f9.3,'   x:',2i5,'    y:',2i5,'   nat',i3,i8)
      call tcx('rxcmro')
      end
