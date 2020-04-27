      subroutine wrhomt(mode,filnam,descr,ib,rhol,rhoc,rofi,nr,nlml,nsp)
C- Writes augmented charge density or potential to file for 1 sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
C          :1 write rhol only
C          :2 write rhoc only
C          :3 write rhol+rhoc
Ci   filnam:file name
Ci   descr :string describing rhol (informational, for output)
Ci   ib    :site index.  If ib>0, append as extension to file name
Ci   rhol  :sphere valence density, tabulated on a radial mesh
Ci         :rl = full charge density * r**2, written as:
Ci         :rl = sum_ilm rhol(ilm) Y_L(ilm)
Ci   rhoc  :sphere core density, tabulated on a radial mesh
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Co   Radial mesh rofi and sphere density rhol are written
Co   to binary file rhoMT.ib in mc-compatible format
Cl Local variables
Cr Remarks
Cr   File is written in format mc program can read.
Cr   To read, use, e.g.
Cr     mc -r:br:open rhoMT.1 -r:br rhoMT.1 -ccat
Cr   To integrate the l=0 density for the sphere charge:
Cr     mc -r:br:open rhoMT.1 -av:nr,1 rmt -a ri -r:br rhoMT.1 -a rhoin
Cr        ri rhoin -ccat -coll 1,2 -int 0 rmt -e2 x1 'x2*sqrt(4*pi)'
Cr   To integrate the l=0 magnetic moment, if lmxl=4
Cr     mc -r:br:open rhoMT.1 -av:nr,1 rmt -r:br rhoMT.1 -coll 1,1+5^2
Cr        -e1 '(x1-x2)*sqrt(4*pi)' -ccat -int 0 rmt
Cr   Integrate the core charge
Cr     mc -r:br:open rhoMT.1 -av:nr,1 rmt -r:br:open rhoMT.1 -pop
Cr        -r:br rhoMT.1 -ccat -int 0 rmt
Cr   A more sophisticated example: integrate rho*(v-2*z/r)
Cr     set rfile =
Cr     '-r:br:open rhoMT.1 -av:nr,1 rmt -a ri -r:br rhoMT.1 -a rin
Cr      -r:br:s=2 vtrue.1 -a v'
Cr     set v1 = 'ri v -ccat -rowl 2:nr -e1 x2-44/x1/0.282094791
Cr               -sub 0,nr,1,1 v -coll 2:nc -ccat -a v1'
Cr     set int = "ri -tog -ccat -int 0 rmt -coll 2:nc -csum"
Cr     mc $rfile rin $int
Cu Updates
Cu   04 Jul 10 Append core density to end of file
Cu   11 Jun 05 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ib,nr,nlml,nsp
      double precision rofi(nr), rhol(nr,nlml,nsp),rhoc(nr,nsp)
      character filnam*(*),descr*(*)
C ... Local parameters
      integer fopng,ip,jfi,flen
      parameter (flen=60)
      character fnam*(flen),strn*(flen)

C ... Get file name
      call words(filnam,ip)
      if (ip /= 1) call rxs('wrhomt : illegal file name: ',filnam)
      call word(filnam,1,jfi,ip)
      fnam = filnam
      if (ib > 0) then
        call bin2a(' ',0,0,ib,2,0,flen,fnam,ip)
      endif
      strn = descr
      call info0(30,0,0,
     .  ' writing sphere '//strn//'%a to file '//fnam//'%a ...')
      if (ip >= flen) call rxs('wrhomt : illegal file name: ',filnam)

C ... dump results
C     Open binary, preserving case
C     call pshpr(120)
      jfi = fopng(fnam,-1,8+4)
C     call poppr
      write (jfi) nr,1,1,0
      write (jfi) rofi
      write (jfi) nr,nlml*nsp,1,0,nsp
      if (mod(mode,2) == 1) write (jfi) rhol
      if (mode >= 2) write (jfi) nr,nsp,1
      if (mode >= 2) write (jfi) rhoc
      call fclose(jfi)

      end
