      subroutine afmevc(mode,nbas,ng,istab,iprmb,q,plat,pos,igapw,afmt,
     .  ndimh,napw,nev,z,zt)
C- Shift the spin-1 eigenvector to (AFM) spin-2 case, using given lattice translation
C ----------------------------------------------------------------------
Ci Inputs
ci   mocde :1: scale evecs by phase tpi * (pj-pi)
Ci   nbas  :size of basis
Ci   ng    :number of group operations
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ndimh :hamiltonian dimension
Ci   napw  :number of augmented PWs in basis
Ci   nev   :actual of eigenvectors
Ci   z     :eigenvectors for spin 1
Co Outputs
Ci   zt    :eigenvectors for AFM spin 2
Cl Local variables
Cl         :
Cr Remarks
Cr
Cr   Check sm densities of two spins identical.  If 30x30x30 divisions:
Cr   mch 1 -coll 1  -a:nr=30 1a 1a -roll:15,15 -a 1b \
Cr       2 -coll 16 -a:nr=30 2b 1b 2b --  -px
Cu Updates
Cu   04 Jan 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,ng,ndimh,napw,nev,igapw(3,napw)
      integer istab(nbas,ng+1),iprmb(*)
      double precision q(3),pos(3,nbas),afmt(3),plat(3,3)
      double complex z(ndimh,nev),zt(ndimh,nev)
C ... Local parameters
      integer n0,nkap0,ib,jb,iv,nlmto,ig
      parameter (n0=10,nkap0=4)
      integer offi,offj,ndimi,ndimj
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      double precision xx,qlat(3,3),qpgv(3)
      double complex phase,tpi
      parameter (tpi=(0d0,6.2831853071795862d0))

      nlmto = ndimh-napw

C     call shoist(0,istab,nbas,xx,xx,ng,1)
C     print *, '!!'; zt = 0

      do  ib = 1, nbas
C       Translation transforms jb into ib
        jb = istab(ib,ng+1)

C       List of orbitals, their l- and k- indices, and ham offsets
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,offi,offl,ndimi)
        call orbl(jb,0,nlmto,iprmb,norb,ltab,ktab,offj,offl,ndimj)
        if (ndimi /= ndimj) call fexit2(-1,111,' Exit -1 AFMEVC: '//
     .    'site %i mapped to site %i, which has inequivalent '//
     .    'hamiltonian',ib,jb)

        phase = 1
        if (mode == 1) then
          xx = (pos(1,ib)-pos(1,jb))*q(1) +
     .         (pos(2,ib)-pos(2,jb))*q(2) +
     .         (pos(3,ib)-pos(3,jb))*q(3)
          phase = exp(tpi * xx)
C         print *, 'ib,phase',ib,phase
        endif
        do  iv = 1, nev
          call zcopy(ndimi,z(1+offj,iv),1,zt(1+offi,iv),1)
          call zscal(ndimi,phase,zt(1+offi,iv),1)
        enddo

      enddo

      if (napw == 0) return
      call dinv33(plat,1,qlat,xx)
      do  ig = 1, napw
C       qpgv(:) = dimag(tpi) * (q + matmul(qlat, igapw(:,ig)))
C       print *, ig,sngl(qpgv(:))
C       qpgv(:) = q + matmul(qlat,igapw(:,ig))
        qpgv(:) = matmul(qlat,igapw(:,ig))
        xx = afmt(1)*qpgv(1) + afmt(2)*qpgv(2) + afmt(3)*qpgv(3)
        phase = exp(-tpi * xx)
C       print *, ig, xx, phase
        call zcopy(nev,z(ig+nlmto,1),ndimh,zt(ig+nlmto,1),ndimh)
        call zscal(nev,phase,zt(ig+nlmto,1),ndimh)
      enddo

C      call zprm('z1',2,z,ndimh,ndimh,nev)
C      call zprm('z2',2,zt,ndimh,ndimh,nev)

      end
