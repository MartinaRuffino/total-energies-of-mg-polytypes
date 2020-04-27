      subroutine gfp0io(mode,fn,nRLc,nsp,nspc,nkp,P0)
C- Writes out response matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 Write real part of P0
Ci         :1 write complex P0
Ci         :10s digit
Ci         :0 write in ascii format
Ci         :1 write in binary format
Ci         :2 read in ascii format
Ci         :3 read in binary format
Ci         :100s digit
Ci         :1 force constant potential shift P0 to have zero eigenvalue
Ci         :  (see prjrsp.f)
Ci   fn    :file name
Ci   nRLc  :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nkp   :number of k-points for which P0 is available
Ci   P0    :response matrix = 1/i G delta v G
Co Outputs
Co    P0 is written to disk
Cl Local variables
Cr Remarks
Cu Updates
Cu   24 May 02 Extended to write P0(q).  Altered argument list
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nRLc,nsp,nspc,nkp
      double complex P0(nRLc,nspc,nRLc,nkp,nsp)
      character*(*) fn
C ... Local parameters
      integer i,opt,ik,isp,iprint,nspx,nRLx,stdo
      character*80 outs
      procedure(integer) :: fopna,lgunit,rdm

C     nspx = number of spin channels
      stdo = lgunit(1)
      nspx = nsp
      if (nspc == 2) nspx = 1
      nRLx = nRLc * nspc
      call sanrg(.true.,mod(mode,10),0,1,'gfp0io:','1s digit mode')

      if (mod(mode/100,10) == 1) call prjrsp(nsp,nspc,2,nRLc,nkp,P0)

C     Binary read/write, open file as binary and set switch for ywrm
      if (mod(mod(mode/10,10),2) == 1) then
        opt = 1
        i = fopna(fn,-1,4)
C     ASCII read/write, open file as ascii and set switch for ywrm
      else
        i = fopna(fn,-1,0)
        opt = 0
      endif
      rewind i

C     call fshow

C --- Read response matrix ---
      if (mod(mode/10,10) >= 2) then
      do  ik = 1, nkp
      do  isp = 1, nspx
        if (mod(mode/10,10) == 2) then
          if (rdm(i,40,nRLc**2*2,' ',P0(1,1,1,ik,isp),nRLc,nRLc) /= 2)
     .      call rx('GFP0IO: failed to read response matrix')
        else
          call rx('gfp0io not ready for binary read')
        endif
      enddo
      enddo

      if (iprint() >= 30) then
        call info2(30,1,0,
     .    ' GFP0IO: read file '//trim(fn)//
     .    '%?#n>1#%-1j for %i k-points##',nkp,0)
      endif

      call fclose(i)
      return

C   99 continue
C      call rx('GFPIO: failed to read file ' // trim(fn))

C --- Write out response matrix ---
      else
      do  ik = 1, nkp
      do  isp = 1, nspx
        call awrit2('%xnsp=%i nspc=%i',outs,80,0,nsp,nspc)
        if (mod(mode,10) == 0) then
          call ztoyy(P0(1,1,1,ik,isp),nRLx,nRLx,nRLx,nRLx,1,0)
          call ywrm(opt,outs,1,i,'(5f15.9)',P0(1,1,1,ik,isp),0,nRLx,
     .      nRLx,nRLx)
          call ztoyy(P0(1,1,1,ik,isp),nRLx,nRLx,nRLx,nRLx,0,1)
        else
          call ywrm(opt,outs,3,i,'(5f15.9)',P0(1,1,1,ik,isp),0,nRLx,
     .      nRLx,nRLx)
        endif
      enddo
      enddo

      if (iprint() >= 30) then
        if (nkp == 1) write(stdo,333) fn
        if (nkp > 1) write(stdo,333) fn, nkp
  333 format(/' GFP0IO: wrote file ',a:' for',i3,' k-points')
      endif

      endif

      call fclose(i)

      end

