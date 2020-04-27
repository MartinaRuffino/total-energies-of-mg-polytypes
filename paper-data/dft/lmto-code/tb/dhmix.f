      subroutine dhmix(neltst,nl,nsp,nbas,idxdn,dclabl,ipc,it,
     .                 nitmax,cnvg,mmix,nkill,beta,tm,tj,dh,delL,delta,
     .                 a,rms)
C- Mixing for TB-L and TB+U
C ----------------------------------------------------------------------
Ci Inputs:
Ci   neltst: total number of elements to mix (for dimensioning a)
Ci   a: the usual Anderson work array
Ci   dh: hamiltonian increments from TBESEL
Co Outputs:
Co   delta: mixed increments, added to H in TBADDH
Co   delL:  mixed increments for writing to disc in IODELL
Cr Remarks
Cr   delL holds the current increments to the non self consistent H
Cr   We think of tbesel as a black box: x_i goes in f(x_i) comes out;
Cr   amix then chooses x_i+1 to go back into tbesel and so on.
Cr   delta(..,i,1) holds the x_i and delta(..,i,2) holds the f(x_i).
Cr   On entry delta(..,i,1) holds the most recent x_i in delta(..,0,1)
Cr   the next most recent in delta(..,1,1) and so on.
Cr   On entry delta(..,i,2) holds the most recent f(x_i) in delta(..,0,1)
Cr   the next most recent in delta(..,1,1) and so on.
Cr   Firstly these are rolled forward so delta(..,0,.) are empty and the
Cr   least recent drop off. Then the latest deltas from tbesel are
Cr   copied into delta(..,0,2) --- these are the f(x_i) tbesel has
Cr   found from the last x_i. Then amix mixes and chooses a new x_i+1
Cr   which is copied into delta(..,i,1).
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer neltst,nl,nsp,nbas,it,nitmax,mmix,nkill,ipc(*),
     .        idxdn(0:nl-1,1)
      double precision cnvg,beta,tm,tj(*),delL(nl**2,nl**2,nbas,nsp),
     .  delta(nl**2,nl**2,nbas,nsp,0:mmix+1,2),dclabl(1),
     .  a(neltst,0:mmix+1,2),rms,dh(nl**2,nl**2,nbas,nsp)
C Local Variables
      integer nelts,npmix,isp,ib,ic,ilm,ilmp,i,j,ido,nmix,
     .  amix,ipr,nproc,mpipid
      integer iprint,i1mach
      real(8), allocatable :: norm(:)
      integer, allocatable :: kpvt(:)
      double precision b,d(9)
      character clabl*8,outs*20
C      logical kill,cmdopt
C      integer LOCIT
C      save LOCIT


      nproc = mpipid(0)
      nelts = nl**4*nbas*nsp

C --- kill "mix files" according to nkill or rms (see parmxp) ---
C      kill = .true.
C      if (cmdopt('--nomixkill',11,0,outs)) then
C        kill = .false.
C      endif

C      Problem with mod(it,nkill) if nkill is zero.
C      if (it == 1) then
C        LOCIT = 0
C      endif
C      if ( ( nkill < 0 .or. rms < 0d0 .or.
C     .     ( nkill > 0 .and. mod(it,nkill) == 0 ) )
C     .    .and. kill ) then
C        LOCIT = 1
C      else
C        LOCIT = LOCIT + 1
C      endif
      npmix = min(it-1,mmix)

C --- Roll back delta's ---
      do  i = mmix, 0, -1
        call dcopy(nelts,delta(1,1,1,1,i,1),1,delta(1,1,1,1,i+1,1),1)
        call dcopy(nelts,delta(1,1,1,1,i,2),1,delta(1,1,1,1,i+1,2),1)
      enddo

C --- Make new delta for this iteration ---
      call dcopy(nelts,0d0,0,delta(1,1,1,1,0,2),1)
      do  ib = 1, nbas
        do  isp = 1, nsp
          ic = ipc(ib)
          do  ilm = 1, nl**2
            do  ilmp = 1, nl**2
              if (it == 1) delta(ilm,ilmp,ib,isp,1,1)
     .                     = delL(ilm,ilmp,ib,isp)
              delta(ilm,ilmp,ib,isp,0,2) = dh(ilm,ilmp,ib,isp)
            enddo
          enddo
          if (iprint() > 60) then
            print *,' Delta previous ..'
            do  i = 1, nl**2
              write (*,300) (delta(i,j,ib,isp,1,1),j=1,nl**2)
            enddo
            print *,' Delta new ..'
            do  i = 1, nl**2
              write (*,300) (delta(i,j,ib,isp,0,2),j=1,nl**2)
            enddo
          endif
        enddo
      enddo
  300 format (10x,9f10.6)
C 310 format (10x,9g12.4)

C --- Build work array for amix ---
      do  i = 0, npmix
        call dcopy(nelts,delta(1,1,1,1,i+1,1),1,a(1,i,2),1)
        call dcopy(nelts,delta(1,1,1,1,i,2),1,a(1,i,1),1)
        if (i /= 0) then
          call daxpy(neltst,-1d0,a(1,i,2),1,a(1,i,1),1)
        endif
      enddo

C --- Mix; don't chatter about it ---
      b = beta
      call pshprt(0)
      ipr = iprint()
      ido = 0
      allocate(norm(mmix*mmix))
      allocate(kpvt(mmix))
      nmix = amix(neltst,npmix,mmix,ido,b,ipr,tm,norm,kpvt,a,tj,rms)
      call popprt
!       call rlse(onorm)
      deallocate(kpvt,norm)
C --- weird bug: processes do not always agree on rms
C      if (nproc > 1) then
C        call mpibc2(rms,1,4,3,.false.,' ',' ')
C        rms = rms / nproc
C      endif

C --- Get new deltas from work array ---
      call dcopy(nelts,a(1,0,2),1,delta(1,1,1,1,0,1),1)

C --- Copy delta's into delL for next iteration of MD or static ---
      if (it == nitmax .or. rms < cnvg)
     .  call dcopy(nelts,delta(1,1,1,1,0,1),1,delL,1)

C --- Printout ---
      if (iprint() < 20) return
      print 100
      call awrit6(
     .  ' Iteration %i. %i elements; mixed %i of %i, beta=%d, '
     .  //'rms diff: %g',' ',90,i1mach(2),it,neltst,nmix,npmix,b,rms)
      if (nmix > 0) write (*,110) (tj(i),i=1,nmix)
      if (iprint() < 40) return
      do  ib = 1, nbas
        do  isp = 1, nsp
          call dcopy(nl**2,delta(1,1,ib,isp,0,1),nl**2+1,d,1)
          ic = ipc(ib)
          call r8tos8(dclabl(ic),clabl)
          if (nsp == 1) then
            if (nl == 3) then
            call awrit1(' Atom '//clabl//'%cdiagonal increments:%9:1d',
     .        ' ',180,i1mach(2),d)
            elseif (nl == 2) then
            call awrit1(' Atom '//clabl//'%cdiagonal increments:%4:1d',
     .        ' ',180,i1mach(2),d)
            endif
          else
            if (nl == 3) then
            call awrit2(' Atom '//clabl//
     .          '%cspin %i diagonal increments:%9:1d',
     .          ' ',180,i1mach(2),isp,d)
            elseif (nl == 2) then
            call awrit2(' Atom '//clabl//
     .          '%cspin %i diagonal increments:%4:1d',
     .          ' ',180,i1mach(2),isp,d)
            endif
          endif
          if (iprint() > 50) then
            print *,' New DeltaH ..'
            if (nsp == 2) then
              if (isp == 1) print*, ' spin up:'
              if (isp == 2) print*, ' spin down:'
            endif
            do  i = 1, nl**2
              write (*,300) (delta(i,j,ib,isp,0,1),j=1,nl**2)
            enddo
          endif
        enddo
      enddo
  100 format(/' DHMIX mixing increments to on-site hamiltonian:')
C 115 format(/' DHMIX mixing increments to off-site hamiltonian:')
  110 format(' t_j :',10f8.4)
C 120 format(1028g12.4)

      end
