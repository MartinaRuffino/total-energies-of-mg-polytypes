      subroutine lgupac(lgii,cp,lmem,ltype,lmod,lbulk,lsymm,lgen)
C- Packs or unpacks the compacted lgii labelling GF type
C ----------------------------------------------------------------
Ci Inputs
Ci   cp: character string indicating what to pack or unpack.
Ci       1st char 'p' or 'u' indicates whether to pack or unpack lgii.
Ci       If it is neither, lgupac writes description of status
Ci       to logical unit lmem; lmem=0 => string copied back to cp.
Ci       For packing and unpacking, caller may (un)pack everything or
Ci       do so selectively.  If no characters follow 'p' or 'u'
Ci       everything is (un)packed.  Otherwise:
Ci       Next character=a => lmem  (un)packed
Ci       Next character=t => ltype (un)packed
Ci       Next character=m => lmod  (un)packed
Ci       Next character=b => lbulk (un)packed
Ci       Next character=s => lsymm (un)packed
Ci       Next character=g => lgen  (un)packed
Ci       Eg 'pbs' packs lbulk and lsymm into lgii; 'p' packs everything
Cio  lgii:  information stored in packed form
Cio  lmem,ltype,lmod,lbulk,lgen
Ci   NB: For printout, lmem is used for the logical unit number.
Cr Remarks
Cr lgii contains in a packed format indicating the
Cr following attributes of the contents of gii, by bits:
Cr  0      (lmem)  1 memory allocated
Cr                 NB: for printout mode, lmem is file logical unit
Cr  1..2   (ltype) 0 no GF calculated
Cr                 1 left  semiinfinite GF
Cr                 2 right semiinfinite GF
Cr                 3 interface GF
Cr  3..5   (lmod)  0 contains GF
Cr                 1 contains GF^-1
Cr                 2 contains GF*S_(i,i+1)
Cr                 3 contains GF*S+_(i,i-1)
Cr                 4 contains S_(i-1,i)*GF
Cr                 5 contains S+_(i+1,i)*GF
Cr  6..9   (lbulk)
Cr         0 not generated from a bulk GF (ie from parts of s00,s0L,s0R)
Cr         1 bulk GF, generated from s00,s0L,s0L+
Cr         2 bulk GF, generated from s00,s0R+,s0R
Cr         3 bulk GF, s0L=s0R+
Cr         Note: "bulk" means potential is periodic in crystal or half-crystal
Cr         4 bulk-like GF but with perturbation to the left
Cr         8 bulk-like GF but with perturbation to the right
Cr         4,8 may be added in combination with 1..3
Cr  10     (lsymm) 1 s0L=s0R+
Cr  11.14  (lgen)
Cr         0 general; no special conditions apply
Cr         1 generating GF on lhs was created from its s00,s0L,s0L+
Cr         2 generating GF on lhs was created from its s00,s0R,s0R+
Cr         4 generating GF on rhs was created from its s00,s0R,s0R+
Cr         8 generating GF on rhs was created from its s00,s0L,s0L+
Cr         Combinations, eg 2+4 are also allowed
Cr  old
Cr  8      (lsymm) 1 s0L=s0R+
Cr  9..12  (lgen)
Cr         0 general; no special conditions apply
Cr         1 generating GF on lhs was created from its s00,s0L,s0L+
Cr         2 generating GF on lhs was created from its s00,s0R,s0R+
Cr         4 generating GF on rhs was created from its s00,s0R,s0R+
Cr         8 generating GF on rhs was created from its s00,s0L,s0L+
Cr         Combinations, eg 2+4 are also allowed
C ----------------------------------------------------------------
      implicit none
      integer lgii,lmem,ltype,lmod,lsymm,lbulk,lgen
      character*(*) cp
      integer ic,lc,l1,l2,l3,l4,l5,ifi,ipr
      logical all
      character outs*80

      call getpr(ipr)
      if (ipr < 60 .and. cp(1:1) == ' ') return

C ... Effective length is position of rightmost nonblank character
      lc = len(cp)
   10 if (cp(lc:lc) /= ' ') goto 12
      lc = lc-1
      if (lc > 1) goto 10
   12 continue
      all = lc <= 1

C --- Packing ---
      if (cp(1:1) == 'p') then
        if (all) then
          lgii = lmem + ltype*2 + lmod*8 + lbulk*64 +
     .            lsymm*1024+lgen*2048
        else
          ic = min(2,lc)
          if (cp(ic:ic) == 'a') then
            lgii = lgii + lmem - mod(lgii,2)
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 't') then
            lgii = lgii + (ltype - mod(lgii/2,4))*2
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 'm') then
            lgii = lgii + (lmod - mod(lgii/8,8))*8
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 'b') then
            lgii = lgii + (lbulk - mod(lgii/64,16))*64
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 's') then
            lgii = lgii + (lsymm - mod(lgii/1024,2))*1024
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 'g') then
            lgii = lgii + (lgen - mod(lgii/2048,16))*2048
            ic = min(ic+1,lc)
          endif
        endif
C --- Unpacking ---
      elseif (cp(1:1) == 'u') then
        if (all) then
          lgen =  lgii/2048
          lsymm = mod(lgii/1024,2)
          lbulk = mod(lgii/64,16)
          lmod =  mod(lgii/8,8)
          ltype = mod(lgii/2,4)
          lmem =  mod(lgii,2)
        else
          ic = min(2,lc)
          if (cp(ic:ic) == 'a') then
            lmem =  mod(lgii,2)
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 't') then
            ltype = mod(lgii/2,4)
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 'm') then
            lmod =  mod(lgii/8,8)
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 'b') then
            lbulk = mod(lgii/64,16)
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 's') then
            lsymm = mod(lgii/1024,2)
            ic = min(ic+1,lc)
          endif
          if (cp(ic:ic) == 'g') then
            lgen =  lgii/2048
            ic = min(ic+1,lc)
          endif
        endif
      else

C --- Printout ---
      lc = len(outs)
      outs = cp
      ifi = lmem
      l1 = mod(lgii/2,4)        ! ltype
      l2 = mod(lgii/8,8)        ! lmod
      l3 = mod(lgii/64,16)      ! lbulk
      l4 = mod(lgii/1024,2)     ! lsymm
      l5 = lgii/2048            ! lgen
      if (mod(lgii,2) == 0 .or. l1 == 0) then
        call awrit1('%a  no GF stored%?#n##, no memory allocated#',
     .    outs,lc,-ifi,mod(lgii,2))
        return
      endif
      call awrit6('%a '//
     .  '%?#n==1# inverse##'//
     .  '%?#n>3# perturbed##'//
     .  '%?#n# bulk##'//
     .  '%?#n==1# left semiinfinite##'//
     .  '%?#n==2# right semiinfinite##'//
     .  '%?#n==3# interface##',
     .  outs,lc,0,l2,l3,l3,l1,l1,l1)
      call awrit5('%a'//
     .  '%?#n<=1# GF##'//
     .  '%?#n==2# S0l*GF##'//
     .  '%?#n==3# GF*S0r##'//
     .  '%?#n#, gen from bulk gll##'//
     .  '%?#n#, gen from bulk grr##',
     .  outs,lc,0,l2,l2,l2,mod(l5,4),mod(l5/4,4))
      call awrit0('%a',outs,-lc,-ifi)

      endif

      end
      subroutine shopgf(mode,strn,nspc,ipl1,ipl2,pgplp,gii,strRL,gend,lgii)
C- Prints out contents of layer GF
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 print out attributes of GF
Ci         :2 print out attributes of GF and contents
Ci         :3 combination of 1+2
Ci   ipl1  :1st principal layer to include in printout
Ci   ipl2  :last principal layer to include in printout
Ci   lgii  :attributes of the gii generated; see lgupac for conventions
Ci   gend  :gend(1) = g(npl,-1)  (lpgf=5)
Ci         :gend(2) = g(-1,npl)  (lpgf=5)
Ci         :gend(3) = gll (Left surface GF at -1)        (lpgf=7)
Ci         :gend(4) = grr (right surface GF at npl)      (lpgf=7)
Co Outputs
Cr Remarks
Cr
Cu Updates
Cu   05 Dec 01 First created
C ----------------------------------------------------------------------
      use structures, only : s_lgf
      implicit none
C ... Passed parameters
      integer mode,nspc,ipl1,ipl2,pgplp(6,-1:ipl2),lgii(-2:ipl2)
      integer strRL(2)
      character *(*) strn
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*),gend(4)
C ... Local parameters
      integer ipl,lgunit,stdo,ldim0,ndg,ldim0x,ndgx,ltype,ld0x,k
      character*80 outs

      if (mode == 0) return
      if (strn /= ' ') call info(20,1,0,trim(strn),0,0)

      stdo = lgunit(1)
      do   ipl = ipl1, ipl2
        if (mod(mode,2) == 1) then
          call pshpr(80)
          outs = ' '
          write(outs,'(i5)') ipl
          call lgupac(lgii(ipl),outs,stdo,0,0,0,0,0)
          call poppr
        endif
        call lgupac(lgii(ipl),'ut',0,ltype,0,0,0,0)
        if (mode >= 2 .and. ltype > 0) then
          ldim0 = pgplp(4,ipl)
          ndg   = pgplp(3,ipl)
          ldim0x= nspc*ldim0
          ndgx  = nspc*ndg
C          call awrit2('%x gii(%i,%i)',outs,len(outs),0,ipl,ipl)
          call yprm('g',2,gii(ipl)%gll,ldim0x*ndgx,ldim0x,ldim0x,ndgx)
        endif

      enddo

      if (mode >= 2 .and. associated(gend(1)%gll)) then
        ld0x = strRL(1)
        ldim0x = nspc*pgplp(4,ipl2)
        ndgx =   nspc*pgplp(4,ipl1)
        call yprm('gRL',2,gend(1)%gll,ld0x*strRL(2),ld0x,ldim0x,ndgx)
      endif

      if (mode >= 2 .and. associated(gend(2)%gll)) then
        ld0x = strRL(2)
        ldim0x = nspc*pgplp(4,ipl1)
        ndgx  = nspc*pgplp(4,ipl2)
        call yprm('gLR',2,gend(2)%gll,ld0x*strRL(1),ld0x,ldim0x,ndgx)
      endif

      do  k  = 3, 4
      if (mode >= 2 .and. associated(gend(k)%gll)) then
        ipl = -1; if (k == 4) ipl = ipl2
        ldim0 = pgplp(4,ipl)
        ndg   = pgplp(3,ipl)
        ldim0x= nspc*ldim0
        ndgx  = nspc*ndg
        call awrit2('%x gs(%i)',outs,len(outs),0,ipl,ipl)
        call yprm(trim(outs),2,gend(k)%gll,ldim0x*ndgx,ldim0x,ldim0x,ndgx)
      endif
      enddo

      end
C      subroutine fmain
C
C      integer lgii,lmem,ltype,lmod,lbulk,lsymm,lgen
C      character*7 ch
C      ch = 'p'
C      lgii = 0
C      lmem   =  1
C      ltype  =  2
C      lmod   =  1
C      lbulk  =  2+4
C      lsymm  =  1
C      lgen   =  1
Cc     call lgupac(3211,'u',lmem,ltype,lmod,lbulk,lsymm,lgen)
C   10 read(*,'(a7)') ch
C      call lgupac(lgii,ch,lmem,ltype,lmod,lbulk,lsymm,lgen)
C      call lgupac(lgii,' testing    ',6,ltype,lmod,lbulk,lgen)
C*      goto 10
C      end
