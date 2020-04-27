      subroutine readsiginp(mode,ifi,nkp,nspin,nev,ilo,ndham,z,sigm)
C- Read static embedded DMFT sigma, convert to orbital basis
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 read nkp,nspin,nev and return
Ci         :1 read sig.inp
Ci         :10s digit
Ci         :0 sig.inp is real
Ci         :1 sig.inp is complex
Ci         :100s digit
Ci         :1 spin-average sig.inp so that sig.inp(+) = -sig.inp(-)
Ci         :1000s digit
Ci         :1 rotate sigm from eigenfunction to orbital basis
Ci   ifi   :file handle
Ci   ilo   :embedded sigm for bands ilo:ilo+nev-1
Ci   ndham:rank of hamiltonian and sigma in orbital basis
Cio Inputs/Ouputs
Cio The following are inputs unless mode=0, in which case they are output
Cio   nkp   :number of irreducible k-points (bzmesh.f)
Cio   nspin :number of spins (1 or 2)
Cio   nev   :actual number of eigenvectors generated
Co Outputs
Co   sigm   :embedded sigma, optionally rotated to orbital basis
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Jan 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ifi,nkp,nspin,nev,ilo,ndham
      complex(8) :: z(ndham,ndham,nspin),sigm(ndham,ndham,nspin)
C ... Local parameters
      integer i1,i2,i,j,isp
      real(8), allocatable :: sigemb(:,:,:,:) ! Local copy of embedded sigma
      complex(8), allocatable :: wk(:,:) ! Work array
      real(8) :: av(2)
      character id*2048  !, fnam*256, first*256
      procedure(logical) :: rdstrn,a2bin
      procedure(integer) :: rdm
      real(8), parameter :: ry2ev = 13.60569193d0

      if (mod(mode,10) == 0) then
        call info2(2,1,-1,' readsiginp mode=%i : count number of k-points',mode,0)
        nkp = 0; nspin = 1; nev = 0
        do  while (rdstrn(ifi,id,len(id),.false.))
          if (len(trim(id)) == 0) cycle
          call word(id,1,i1,i2)
          if (id(1:1) /= '#') then
            call words(id,i)
            if (nev == 0) nev = i
            if (nev /= i) call rx('mismatch reading sig.inp')
            cycle
          endif
          call word(id,2,i1,i2)
          if (i2 < i1) cycle
          call locase(id(i1:i2))
          if (id(i1:i2) == "kpt") nkp = nkp+1
          if (id(i1:i2) == "spin") then
            call word(id,4,i1,i2)
            if (id(i1:i2) == "2") nspin = 2
          endif
        enddo
        call info5(0,0,0,' ... found %i k-points, %i bands%?#n==2#, 2 spins##',nkp,nev,nspin,0,0)
        return

      elseif (mod(mode,10) == 1) then
        allocate(sigemb(nev,nev,2,2))   ! Read as imaginary following real
        call dpzero(sigemb,size(sigemb))
        isp = 0
        do  while (rdstrn(ifi,id,len(id),.false.))
          if (len(trim(id)) == 0) cycle
          call word(id,2,i1,i2)
          if (i2 < i1) cycle
          call locase(id(i1:i2))
          if (id(i1:i2) == "kpt") then
            isp = isp+1
            cycle
          endif
          if (id(i1:i2) /= "spin") call rx('mismatch reading sig.inp')
          call word(id,4,i1,i2)
          call rxx(i1 /= i2,'mismatch reading sig.inp')
          call rxx(iachar(id(i1:i1))-iachar('0') /= isp,'mismatch reading sig.inp')

C         Read real sigma only for now
          if (rdm(ifi,1000,nev*nev,' ',sigemb(1,1,1,isp),nev,nev) /= 1) call rxi(
     .      'readsiginp (abort): could not read sigma from file, mode',mode)
          if (isp == nspin) exit
          isp = isp+1

        enddo

C   ... Average up and down spins
        if (nspin == 2 .and. mod(mode/100,10) == 1) then
          do  i = 1, nev
            do  j = 1, nev
              av(:) = sigemb(i,j,:,1)/2 + sigemb(i,j,:,2)/2
              sigemb(i,j,:,1) = sigemb(i,j,:,1) - av(:)
              sigemb(i,j,:,2) = sigemb(i,j,:,2) - av(:)
            enddo
          enddo
        endif

C   ... eV to Ry
        call dscal(size(sigemb),1/ry2ev,sigemb,1)

C       print *, 'j'; j=0

C   ... Convert to complex form and copy to sigm
        call dpzero(sigm,2*size(sigm))
        do  i = 1, nspin
          call ztoyy(sigemb(1,1,1,i),nev,nev,nev,nev,0,1)
          call zmscop(0,nev,nev,nev,ndham,0,0,0,0,sigemb(1,1,1,i),sigm(ilo,ilo,i))
C         if (j == 1) call zprm('sigm(emb) eigenfunction basis one spin',2,sigm(1,1,i),ndham,ndham,ndham)
        enddo
        deallocate(sigemb)

C   ... Rotate to orbital basis
C       mch z1 -a z z -i -cc -t sige1 -x z -i -x sigo1 -- -px
        if (mod(mode/1000,10) == 1) then
          allocate(wk(ndham,ndham))
          do  i = 1, nspin
C           if (j == 1) call zprm('z',2,z(1,1,i),ndham,ndham,ndham)
            call phmbls(33,ndham,ndham,[0d0],[0],wk,sigm(1,1,i),z(1,1,i),z(1,1,i),sigm(1,1,i))
C           if (j == 1) call zprm('zi',2,z(1,1,i),ndham,ndham,ndham)
C           if (j == 1) call zprm('sigm(emb) orbital basis one spin',2,sigm(1,1,i),ndham,ndham,ndham)
          enddo
          deallocate(wk)
        endif

      else
        call rxi('readsiginp: invalid mode',mode)

      endif

      end
