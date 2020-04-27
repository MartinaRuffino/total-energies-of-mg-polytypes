      subroutine veltet(n1,n2,n3,eband,nkp,nbmax,nsp,ef,ntet,idtet,
     .  nfilo,nfiup,nfilm,velmt)
C ----------------------------------------------------------------------
      implicit none
      integer n1,n2,n3,nkp,nfilm,nbmax,nsp,idtet(0:4,*),ntet
      integer nfilo,nfiup
      double precision velmt(3,nfilm,2,*),ef
C Local variables
      integer ib,ibf,iq1,iq2,iq3,iq4,k,isp,npts,iii,itet,lgunit
      parameter (npts=10)
      double precision eband(nbmax,nsp,nkp),emin,emax,vel(3,2)
      double precision volwgt,wtm,wt,ds2,eigen(4),zos(npts+1)

      emin = ef
      emax = ef+1d0
      volwgt = dble(3-nsp)/(n1*n2*n3*6)
      ds2 = dsqrt(2d0)*1d-6

C --- Loop over spins and polarizations ---
      do  10  isp = 1, nsp
      do  10  k = 1, 3

        vel(k,isp) = 0d0

C   --- Loop over bands ---
        ibf = 0
        do  20  ib = nfilo, nfiup
        ibf = ibf+1
        do  22  iii = 1, 10
   22   zos(iii) = 0d0

C   --- Loop over tetrahedra ---
        do  11  itet = 1, ntet
          iq1 = idtet(1,itet)
          iq2 = idtet(2,itet)
          iq3 = idtet(3,itet)
          iq4 = idtet(4,itet)

C     ... Set up energies at 4 corners of tetrahedron for ib
          eigen(1) = eband(ib,isp,iq1)  + ds2*0d0
          eigen(2) = eband(ib,isp,iq2)  + ds2*1d0
          eigen(3) = eband(ib,isp,iq3)  + ds2*2d0
          eigen(4) = eband(ib,isp,iq4)  + ds2*3d0

          wtm=
     .      (velmt(k,ibf,isp,iq1)+velmt(k,ibf,isp,iq2) +
     .       velmt(k,ibf,isp,iq3)+velmt(k,ibf,isp,iq4))
          wt = volwgt*idtet(0,itet)*wtm

C         call slinz(wt,eigen,emin,emax,zos(1),npts)
          call sliny(wt,eigen,emin,emax,zos(1),npts)

   11   continue

        vel(k,isp) = vel(k,isp) + zos(1)

C ... End loops over bands, polarizations, spins
   20 continue
   10 continue

      print *, ' '
      do  42  isp = 1, nsp
   42 call awrit2(' VELTET:  spin %i  velocities = %3:1,4d',
     .  ' ',180,lgunit(1),isp,vel(1,isp))

      end
