      subroutine dvsxl(nbas,ldim,aintra,offH,indxsh,vsxl,vsxnl)
C- Subtract local potential dE_sx/d_n
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas,aintra
Ci   offH:      orbital downfolding and ham. offsets info
Ci   vsxnl:     nonlocal screened exchange
Ci   vsxl:      local part of vsxnl
Co Outputs
Ci   vsxnl:     vsxl subtracted from vsxnl
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
      logical aintra
      integer nkap0,n0H,n0
      parameter (nkap0=4,n0H=5,n0=10)
      integer ldim,nbas,offh(n0H,nkap0,nbas),indxsh(*)
      double precision vsxnl(ldim,ldim),vsxl(nbas)
      integer ib,ilb,norbb,ntorbb,lb,mb,offb
      integer ltabb(2*n0),ktab(2*n0),offl(2*n0)

C     call yprm('vsxnl(init)',1,vsxnl,0,ldim,ldim,ldim)

      do  10  ib = 1, nbas
C       uses norbb,ltabb
        call orbl(ib,0,ldim,indxsh,norbb,ltabb,ktab,offb,offl,ntorbb)
        offb = offH(1,1,ib)
        do  20  ilb = 1, norbb
        lb = ltabb(ilb)
        do  20  mb = -lb, lb
          offb = offb+1
          vsxnl(offb,offb) = vsxnl(offb,offb) - vsxl(ib)
   20   continue
   10 continue

C     call yprm('vsxnl(- vsxl)',1,vsxnl,0,ldim,ldim,ldim)
      end
