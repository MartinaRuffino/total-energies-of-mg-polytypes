      subroutine lmorder(opt,lmax,l,m)
C- Ordering of l and m quantum numbers for spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0  Return in lmax, defaultlm = the default m ordering of spherical harmonics
Ci         :   Note: defaultlm may be modified  by calling setmorderdefault
Ci         :   defaultlm (specifies default ordering)
Ci         :       1         m = -l,-l+1, ..., l
Ci         :       0         m =  l,l-1, ..., -l
Ci         :       2         m =  l,l-1, ..., -l and also Ylm deviate from standard convention by (-1)^m.
Ci         :1  Return l(1:lmax+1)^2 and m(1:lmax+1)^2  ordered -l, -l+1, ..., l
Ci         :2  Return l(1:lmax+1)^2 and m(1:lmax+1)^2  ordered  l,  l-1, ...,-l
Ci         :<0 Return l and m, using defaultlm to select m ordering
Cio Inputs/Outputs
Cio   lmax :Input if opt=0.  Returns default Questaal convention.  l and m are not returned
Cio        :Output if opt nonzero.  l-cutoff for arrays l and m
Co Outputs
Co   l     :l(1:(lmax+1)**2) = l quantum number
Co   m     :l(1:(lmax+1)**2) = m quantum number or -m
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 May 18 New routine written to establish uniform convention for Questaal package
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,lmax,l((lmax+1)**2),m((lmax+1)**2)
C ... Local parameters
      integer il,im,ilm,lopt
      integer :: defaultlm = 0  ! Must always be 0, 1, 2, or 3.  See opt above

C ... Manage defaults
      lopt = opt
      if (opt <= 0) lopt = defaultlm
      if (opt == 0) then
        lmax = lopt
        return
      endif

C ... Make l and m for each ilm
      ilm = 0
      do  il = 0, lmax
        do  im = -il, il
          ilm = ilm+1
          l(ilm) = il
          m(ilm) = im
          if (lopt /= 1 .and. lopt /= 3) m(ilm) = -im
        enddo
      enddo

      return

C ... Set default ordering
      entry setmorderdefault(opt)
      call sanrg(.true.,opt,0,2,'setmorderdefault:','mode')
      defaultlm = opt

      end

C      subroutine pmorder(mode,ib1,ib2,nspc,nl,indxsh,ldimp,ldima,lds,lds2,ldi,ldj,s)
CC- Permute m ordering of structure matrix or vector
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :1s digit for row ordering
CCi         :0 do nothing
CCi         :1 permute m ordering
CCi         :10s digit for column ordering
CCi         :0 do nothing
CCi         :1 permute m ordering
CCi   nbas  :size of basis
CCi   nl    :(global maximum l) + 1
CCi   indxsh:permutation indices ordering orbitals in downfolding order
CCi   ldimp :hamiltonian block consists of orbitals betw. ldimp and ldima
CCi   ldima :offset to last orbital in block to be permuted
CCi         :Only elements i with ldimp<indxsh(i)<=ldima are rotated.
CCi   lds   :leading dimensions of s
CCi   lds2  :second dimensions of s
CCi   ldi   :number of rows to permute, when permuting columns
CCi   ldj   :number of columns to permute, when permuting rows
CCi   s     :
CCio Inputs/Outputs
CCo   s     :structure matrix or vector to be permuted. Permuted on output
CCl Local variables
CCr Remarks
CCu Updates
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,ldima,ldimp,lds,lds2,ldi,ldj,ib1,ib2,nspc,nl,indxsh(lds)
C      double complex s(lds,nspc,lds2,nspc)
CC ... Local parameters
C      logical lmordering
C      integer lmi,ib,i,ip,il,im,j,ilm,mode0,mode1,is1,is2,mxorb
C      double complex zx
C      procedure(integer) :: nglob
C
CC     call zprm('pmorder s before perm',2,s,lds,ldi,ldj)
C
C      mode0 = mod(mode,10)
C      mode1 = mod(mode/10,10)
C      lmordering = indxsh(1) == 0
C      mxorb = nglob('mxorb')
C
CC --- Permute Rows ---
C      if (mod(mode,10) == 1) then
CC       lmi = (ib1-1)*mxorb
C        lmi = 0
C        do  ib = ib1, ib2
C          ilm = 0
C          do  il = 0, nl-1
C          do  im = -il, il
C          lmi = lmi+1; ilm = ilm+1; i = lmi
C          if (im <= 0) cycle
C          if (.not. lmordering) i = indxsh(lmi)
C          ip = i - 2*im
C          if (i <= ldimp .or. i > ldima) cycle
C
CC         call info2(1,0,0,' ib l m i ip = %5:1,3i',[ib,il,im,i,ip],2)
C
C          do  j = 1, ldj
C            do  is1 = 1, nspc
C            do  is2 = 1, nspc
C              zx = s(i,is1,j,is2)
C              s(i,is1,j,is2) = s(i-2*im,is1,j,is2)
C              s(i-2*im,is1,j,is2) = zx
C            enddo
C            enddo
C          enddo
C          enddo
C          enddo
C        enddo
C      endif
C
CC     call zprm('pmorder s after row perm',2,s,lds,ldi,ldj)
C
CC --- Permute Columns ---
C      if (mod(mode/10,10) == 1) then
CC       lmi = (ib1-1)*mxorb
C        lmi = 0
C        do  ib = ib1, ib2
C          ilm = 0
C          do  il = 0, nl-1
C          do  im = -il, il
C          lmi = lmi+1; ilm = ilm+1; i = lmi
C          if (im <= 0) cycle
C          if (.not. lmordering) i = indxsh(lmi)
C          ip = i - 2*im
C          if (i <= ldimp .or. i > ldima) cycle
C
CC         call info2(1,0,0,' ib l m i = %5:1,3i',[ib,il,im,i,ip],2)
C
C          do  j = 1, ldi
C            do  is1 = 1, nspc
C            do  is2 = 1, nspc
C              zx = s(j,is1,i,is2)
C              s(j,is1,i,is2) = s(j,is1,i-2*im,is2)
C              s(j,is1,i-2*im,is2) = zx
C            enddo
C            enddo
C          enddo
C          enddo
C          enddo
C        enddo
C      endif
C
CC     call zprm('pmorder s after col perm',2,s,lds,ldi,ldj)
C
C      end

      subroutine pmorderx(mode,ib1,ib2,nspc,nl,iprmb,ldimp,ldima,lds,lds2,ldi,ldj,off,s)
C- Permute m ordering of structure matrix or vector
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit for row ordering
Ci         :0 do nothing
Ci         :1 permute m ordering
Ci         :10s digit for column ordering
Ci         :0 do nothing
Ci         :1 permute m ordering
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   iprmb:permutation indices ordering orbitals in downfolding order
Ci   ldimp :hamiltonian block consists of orbitals betw. ldimp and ldima
Ci   ldima :offset to last orbital in block to be permuted
Ci         :Only elements i with ldimp<iprmb(i)<=ldima are rotated.
Ci   lds   :leading dimensions of s
Ci   lds2  :second dimensions of s
Ci   ldi   :number of rows to permute, when permuting columns
Ci   ldj   :number of columns to permute, when permuting rows
Ci   off   :subtract row offset off when permuting rows
Ci         :subtract column offset off when permuting columns
Ci         :off should be the hamiltonian offset to site ib1 when
Ci         :indexing to s starts at 1
Cio   s    :structure matrix or vector to be permuted. Permuted on output
Cr Remarks
Cr    Either rows, or columns, or both may be permuted.
Cu Updates
Cu      May 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldima,ldimp,lds,lds2,ldi,ldj,ib1,ib2,nspc,nl,off,iprmb(lds)
      double complex s(lds,nspc,lds2,nspc)
C ... Local parameters
      logical lmordering
      integer lmi,ib,i,ip,il,im,j,mode0,mode1,is1,is2,mxorb
      double complex zx
      procedure(integer) :: nglob

C     call zprm('pmorder s before perm',2,s,lds,ldi,ldj)

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      lmordering = iprmb(1) == 0
      mxorb = 0; if (ib1>1) mxorb = nglob('mxorb')

C --- Permute Rows ---
      if (mode0 == 1) then
        lmi = (ib1-1)*mxorb
        do  ib = ib1, ib2
          do  il = 0, nl-1
          do  im = -il, il
          lmi = lmi+1
          if (im <= 0) cycle
          i = lmi-off; if (.not. lmordering) i = iprmb(lmi)-off
          ip = i - 2*im
          if (i <= ldimp .or. i > ldima) cycle

C         call info2(1,0,0,' ib l m i ip = %5:1,3i',[ib,il,im,i,ip],2)

          do  j = 1, ldj
            do  is1 = 1, nspc
            do  is2 = 1, nspc
              zx = s(i,is1,j,is2)
              s(i,is1,j,is2) = s(i-2*im,is1,j,is2)
              s(i-2*im,is1,j,is2) = zx
            enddo
            enddo
          enddo
          enddo
          enddo
        enddo
      endif

C     call zprm('pmorder s after row perm',2,s,lds,ldi,ldj)

C --- Permute Columns ---
      if (mode1 == 1) then
        lmi = (ib1-1)*mxorb
        do  ib = ib1, ib2
          do  il = 0, nl-1
          do  im = -il, il
          lmi = lmi+1
          if (im <= 0) cycle
          i = lmi-off; if (.not. lmordering) i = iprmb(lmi)-off
          ip = i - 2*im
          if (i <= ldimp .or. i > ldima) cycle

C         call info2(1,0,0,' ib l m i = %5:1,3i',[ib,il,im,i,ip],2)

          do  j = 1, ldi
            do  is1 = 1, nspc
            do  is2 = 1, nspc
              zx = s(j,is1,i,is2)
              s(j,is1,i,is2) = s(j,is1,i-2*im,is2)
              s(j,is1,i-2*im,is2) = zx
            enddo
            enddo
          enddo
          enddo
          enddo
        enddo
      endif

C     call zprm('pmorder s after col perm',2,s,lds,ldi,ldj)

      end
