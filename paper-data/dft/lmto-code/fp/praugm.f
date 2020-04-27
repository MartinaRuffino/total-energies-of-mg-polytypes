      subroutine praugm(s_spec,is)
C- Print species information
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name rmt rsma lmxa kmxt lmxl rg rsmv kmxv lfoca rfoca
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   is    :species index (use 0 to print info for all species)
Co Outputs
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer is
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer is1,is2,js,kmxt,lmxa,nglob,lgunit,stdo
      integer kmxv,lmxl,lfoca
      double precision rmt,rsma,rfoca,rg,rsmv
      character spid*8
C      parameter (n0=10)
C      integer idmod(n0),l
C      double precision pnu(n0*2)
C      character*1  s2*80, idcode(0:1)

      stdo = lgunit(1)
C     idcode(0) = ' '
C     idcode(1) = '*'
      is1 = is
      is2 = is
      if (is <= 0) then
        is1 = 1
        is2 = nglob('nspec')
      endif

      write (stdo,501)
  501 format(/' species data:  augmentation',27x,'density'/
     .     ' spec       rmt   rsma lmxa kmxa',5x,
     .     ' lmxl     rg   rsmv  kmxv foca   rfoca')
      do  js = is1, is2
        spid = s_spec(js)%name
        rmt = s_spec(js)%rmt
        rsma = s_spec(js)%rsma
        lmxa = s_spec(js)%lmxa
        kmxt = s_spec(js)%kmxt
        lmxl = s_spec(js)%lmxl
        rg = s_spec(js)%rg
        rsmv = s_spec(js)%rsmv
        kmxv = s_spec(js)%kmxv
        lfoca = s_spec(js)%lfoca
        rfoca = s_spec(js)%rfoca

C         pnu = s_spec(js)%p
C         idmod = s_spec(js)%idmod
C        write (s2,102) (pnu(l+1),idcode(idmod(l+1)),l=0,lmxa)
C  102   format(10(f5.2,a1))
C        write (stdo,300) spid,rmt,rsma,lmxa,kmxt,s2(1:36)
C  300   format(1x,a,f6.3,f7.3,i6,i6,2x,a)
        write (stdo,500) spid,rmt,rsma,lmxa,kmxt,
     .                   lmxl,rg,rsmv,kmxv,lfoca,rfoca
  500   format(1x,a,f6.3,f7.3,2i5,6x,i4,2f7.3,i6,i5,f8.3)

      enddo
      end

