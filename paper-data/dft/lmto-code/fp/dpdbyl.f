      subroutine dpdbyl(a,nr,nlm,nlm0,nsp0,nsp,lbin,ifi)
C- Dumps or reads an array given as a(nr,nlm,nsp).
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlm   :number of components to be i/o from file.
Ci         :If nlm is greater than nlm0,
Ci         :higher components are read and discarded.
Ci   nlm0  :the array second dimension
Ci   nsp0  :number of spins contained in file
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio   a    :array a(1..nr,1..nlm,1..nsp) is read or written to file
Cr Remarks
Cu Updates
Cu   04 Feb 05 Spin-splits and unspin-polarized density
Cu   27 Apr 01 Added lbin switch
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lbin
      integer nr,nlm,nlm0,nsp0,nsp,ifi
      double precision a(nr,nlm0,nsp)
C ... Local parameters
      logical isanrg
      integer jfi,isp,ilm,i,nlmx
      double precision xx

C --- Output ---
      if (ifi < 0) then
        if (isanrg(nsp0,nsp,nsp,'dpdbyl writing rho','nsp',.true.)) stop
        jfi = -ifi
        do  isp = 1, nsp
          do  ilm = 1, nlm
            if (lbin) then
              write(jfi) (a(i,ilm,isp),i=1,nr)
            else
              call dfdump(a(1,ilm,isp),nr,-jfi)
            endif
          enddo
        enddo
        return
      endif

C --- Input ---
      jfi = ifi
      nlmx = min0(nlm,nlm0)

C ... loop over spins
      do  isp = 1, nsp

C ...   read the desired components
        do  ilm = 1, nlmx

          if (isp == 2 .and. nsp0 == 1) then
            call dscal(nr,0.5d0,a(1,ilm,1),1)
            call dpscop(a(1,ilm,1),a(1,ilm,2),nr,1,1,1d0)
          elseif (lbin) then
            read(jfi) (a(i,ilm,isp), i=1,nr)
          else
            call dfdump(a(1,ilm,isp),nr,jfi)
          endif
        enddo

C ...   read and discard higher components in file
        do  ilm = nlmx+1, nlm
          if (isp == 2 .and. nsp0 == 1) then
          elseif (lbin) then
            read(jfi) (xx, i=1,nr)
          else
            read(jfi,*) (xx, i=1,nr)
          endif
        enddo

C ...   zero out unset components in array
        do  ilm = nlmx+1, nlm0
          call dpzero(a(1,ilm,isp),nr)
        enddo

      enddo

      end


