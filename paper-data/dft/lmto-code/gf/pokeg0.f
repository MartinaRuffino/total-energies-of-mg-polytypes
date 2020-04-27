      subroutine pokeg0(mode,norb,iprm,ldg1,ldg2,lidim,nspc,nsp,gii,g)
C- Push/poke between a local matrix gii and global matrix g (H order)
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 insert data from gii into g
Ci         :1 extract from g into gii
Ci   norb  :leading dimension of gii, and number of orbitals in it
Ci   iprm  :permutation table for this site.
Ci         :If iprm(1) = 0 then table is not used.  Orbitals are ordered 1,2,3,...
Ci   ldg1  :leading dimension of g
Ci   ldg2  :second dimension of g
Ci   lidim :upper bound to lower+intermediate block
Ci   nspc  :number of coupled spins
Cio Inputs/Outputs
Cio  gii   :site-local matrix
Cio  g     :global matrix
Cu Updates
Cu   03 Dec 16 (MvS) Added nsp argument: reads/writes one or two spins
Cu   09 Jun 14 (Kirill) First created from earlier version
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,intent(in) :: mode,norb,iprm(norb),lidim,ldg1,ldg2,nspc,nsp
      complex(8),intent(inout) :: gii(norb,nspc,norb,nsp),
     .  g(ldg1,nspc,ldg2,nsp)
C ... Local parameters
      logical lordered
      integer n1,n2,ip1,ip2

      lordered = iprm(1) == 0

      do  n2 = 1, norb
        do  n1 = 1, norb
          ip1 = n1 ; ip2 = n2
          if (.not. lordered) then
            ip1 = iprm(n1) ; ip2 = iprm(n2)
          endif
          if (ip1 > lidim .or. ip2 > lidim) cycle
          if (mode == 0) then
            g(ip1,:,ip2,:) = gii(n1,:,n2,:)
          elseif (mode == 1) then
            gii(n1,:,n2,:) = g(ip1,:,ip2,:)
          endif
        enddo
      enddo

      end
