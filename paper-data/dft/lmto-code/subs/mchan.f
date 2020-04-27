      subroutine mchan(lmdim,s_site,s_spec,nsp,nsites,
     .  lsites,ib,ilm,io,ichan,lchan)
C- set or get Mulliken channel for site ib and ilm index
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs:
Ci   lmdim,nsites,ib,ilm,io < 0 poke ichan into lchan
Ci                          > 0 get ichan from lchan
Ci                          = 0 print channel table
Co Outputs
Co   ichan :(io > 0) lchan(ilm,ib) for supplied (ilm,ib) pair
Co         :Otherwise, ichan is not set
Co   lchan :(io < 0) set lchan(ilm,ib) to ichan
Co         :Otherwise, lchan is not set
Cr Remarks
Cr    For the Mulliken decomposition it is convenient to keep a table
Cr    lchan(ilm,ib) which holds the DOS channel number associated with
Cr    site ib, and lm channel ilm. If the DOS is site projected then
Cr    ilm is 1; if l projected then ilm=1,lmax+1; if lm-projected then
Cr    ilm=1,(lmax+1)**2. The leading dimension lmdim hence depends on
Cr    the mode (see sumlst) in this way.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Mar 10 (MvS) handle case ichan>lmdim
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lmdim,nsp,nsites,lsites(nsites),lchan(lmdim,nsites),
     .        ib,ilm,io,ichan
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer i,j,lgunit,js,jb
      character clabl*8

      if (io < 0) then
        if (ilm <= lmdim) lchan(ilm,ib) = ichan
      elseif (io > 0) then
        if (ilm > lmdim) then
          ichan = 0
        else
          ichan = lchan(ilm,ib)
        endif
      else
        call awrit1('%N mchan: channel table, %i sites',
     .    ' ',256,lgunit(1),nsites)
        if (nsp == 2)
     .  call awrit0(' (each channel splits into two: up, down spin)',
     .    ' ',256,lgunit(1))
        call awrit0(' site  label          channels',' ',128,lgunit(1))
        do  j = 1, nsites
          jb = lsites(j)
          js = s_site(jb)%spec
          clabl = s_spec(js)%name
          write (lgunit(1),1) j, clabl, (lchan(i,j),i=1,lmdim)
        enddo
      endif
    1 format (1x,i3,5x,a8,1x,256i3)
      end
