      subroutine phmbls(mode,ndz,nev,eval,iprm,wk,h,zd,z,zhz)
C- Make (zd+)*h*z; order in ascending order of eigenvalues
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :consists of a sequence of bits that prescribed independent
Ci         :functions
Ci         :0  =mode0   do nothing; just return
Ci         :1  =mode1   make zhz = zd+ * h * z.  zhd has dimension (ndzhz,ndzhz)
Ci         :2  =mode2   poke diagonal of zhz into eval
Ci         :4  =mode4   permute eval in ascending eval order
Ci         :8  =mode8   permute z in ascending eval order
Ci         :16 =mode16  permute zhz in ascending eval order
Ci         :32 =mode32  Overwrite input z with z^-1; use z^-1 in place of z
Ci         :64 =mode64  make zhz = zd+ * e * z
Ci         :some modes can be used in combination, e.g.
Ci         : mode32 + mode1  => (zd+)^-1 * h * (z)^-1
Ci         : mode32 + mode64 =>  reconstructs hamiltonian from evals,evecs.
Ci   ndz   :leading dimension of zd, z, zhz, wk, h
Ci   nev   :Modes (1,2,4,8,16) can operate in a reduced subspace of h,z.
Ci         :In that case, nev = rank of reduced subspace.
Ci   wk    :double complex work array of length ndz**2
Ci   h     :hamiltonian
Ci         :Not used unless mode1 set
Cio Inputs/Outputs
Cio  iprm  :integer work of length ndz.
Cio        :Space used if any permutations done (mode4,mode8,mode16)
Cio        :On output, holds permutations that order z and zhz in
Cio        :ascending eval order.
Cio  eval  :eigenvalues (dimensioned nev).
Cio        :eval is only accessed if mode contains one of these bits: 2,4,8,16,64
Cio        :  bit   purpose
Cio        :   1    -- not used
Cio        :   2    calculated from zHz
Cio        :   4    sorted
Cio        :   8    used to sort z
Cio        :  16    used to sort zHz
Cio        :  32    -- not used
Cio        :  64    used to make zd+ * e * z
Cio  zd    :(mode1,mode64 only) left evec to generate (zd+ h z).
Cio        :Usually zd and z point to the same address space;
Cio        :however, they may be different, e.g. in constructing
Cio        :off-diagonal spin blocks from spin diagonal z's
Cio  z     :approximate eigenvectors
Cio        :If mode16, z are permuted by ascending eval order
Cio  zhz   :zd+ * h * z, where z are evecs, h is hamiltonian
Cio        :If mode1, zhz is generated from z
Cio        :Otherwise zhz is input
Cio        :If mode8, output zhz is permuted by ascending eval order
Cio        :NB: zhz can occupy the same address space as h
Cio        :in which case h will be overwitten.
Cio        :Note, however, that zhz and h may be dimensioned differently
Cr Remarks
Cr
Cu Updates
Cu   16 Jan 07 New mode 32
Cu   19 May 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ndz,nev,iprm(ndz)
      double precision eval(ndz)
      double complex z(ndz,ndz),zd(ndz,ndz),zhz(ndz,ndz)
      double complex wk(ndz,ndz),h(ndz,ndz)
C ... Local parameters
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      integer i,j,k
      integer ipiv(ndz)

      if (mode == 0) return

C ... Invert z
      if (mod(mode/32,2) /= 0) then
        call zgetrf(ndz,ndz,z,ndz,ipiv,j)
        call rxx(j /= 0,'hambls: failed to invert Z')
        call zgetri(ndz,z,ndz,ipiv,wk,ndz**2,j)
C       call zprm('z^-1',2,z,ndz,ndz,ndz)
      endif

C ... zhz <- zd+ * h * z
C     Debug:  mc -f9f18.11 zd -cc -t h -x z -x
      if (mod(mode,2) == 1) then
        call zgemm('C','N',nev,nev,nev,one,zd,ndz,h,ndz,zer,wk,ndz)
        call zgemm('N','N',nev,nev,nev,one,wk,ndz,z,ndz,zer,zhz,ndz)
C       call zprm('z+ h z',2,zhz,ndz,nev,nev)
      endif

C ... Eigenvalue estimates
      if (mod(mode/2,2) == 1) then
        do  i = 1, nev
          eval(i) = zhz(i,i)
        enddo
      endif

C ... Find permutation iprm that orders eigenvalues in ascending order
      if (mod(mode/16,2) /= 0 .or. mod(mode/8,2) /= 0 .or. mod(mode/4,2) /= 0) then
        call dvheap(1,nev,eval,iprm,0d0,1)
      endif

C ... Permute evals eigenvalues in ascending order
      if (mod(mode/4,2) /= 0) then
        call dvprm(1,nev,eval,wk,iprm,1)
      endif

C ... Permute eigenvectors z in ascending eval order
      if (mod(mode/8,2) /= 0) then
      call zcopy(ndz**2,z,1,wk,1)
      do  i = 1, nev
        k = iprm(i)
        do  j = 1, ndz
          z(j,i) = wk(j,k)
        enddo
      enddo
      endif

C ... Permute zhz in ascending eval order
      if (mod(mode/16,2) /= 0) then
      call zcopy(ndz**2,zhz,1,wk,1)
      do  i = 1, nev
        k = iprm(i)
        do  j = 1, ndz
          zhz(j,i) = wk(iprm(j),k)
        enddo
      enddo
      endif

C ... zhz <- zd+ * e * z ... usu. used in conjunction w/ mode=32
      if (mod(mode/64,2) /= 0) then
        do  j = 1, ndz
        do  i = 1, nev
          wk(i,j) = eval(i)*z(i,j)
        enddo
        do  i = nev+1, ndz
          wk(i,j) = 0
        enddo
        enddo
C       call zprm('e * z',2,wk,ndz,ndz,ndz)
        call zgemm('C','N',ndz,ndz,ndz,one,zd,ndz,wk,ndz,
     .    zer,zhz,ndz)
C       call zprm('z+ h z',2,zhz,ndz,ndz,ndz)

      endif

      end
