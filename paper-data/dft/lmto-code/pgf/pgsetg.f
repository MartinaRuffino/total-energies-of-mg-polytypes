      subroutine pgsetg(mode,ldmx,ndgn,ld0,ndg0,offi,nspc,gii,gin)
C- Copy gin to/from gii with different dimensions
C---------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy gii to gin
Ci         :1 copy gin to gii
Ci   ldmx  :leading dimension of gin, usually
Ci         :maximum leading dimension of gii(ipl) for ipl=-1,0,..,npl
Ci   ndgn  :column dimension of gin
Ci   ld0   :leading dimension and rank of (local) gii
Ci   ndg0  :column dimension of (local) gii (usuall ndg0=ld0)
Ci   offi  :offset in column dimension for local g
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Cio Inputs/Outputs
Cio  gii   :local g in format: gii(ld0,ndg0,2)
Cio        :source (mode=0) or destination (mode=1)
Cio  gin   :fixed dimension g in format: gin(lldmx,ndg0,2)
Cio        :source (mode=1) or destination (mode=0)
Cr Remarks
Cr    This routine convert g represented in gin, a fixed, layer-independent size
Cr    to/from gii sized for the local layer.
Cr
Cr    gin can be re-used as PL change but it inconvenient for matrix operations
Cr    Spin blocks in gii are contiguous : gii=gii(ld0*nspc,ndg0*nspc)
Cr
Cr                      ndg0*nspc
Cr            .----------------------------.
Cr            |                            |
Cr       ld0  |                            |
Cr      *nspc |                            |
Cr            |                            |
Cr            |                            |
Cr            .----------------------------.
Cr      While spin blocks of gin need not be:
Cr                 ndgn              ndgn
Cr            .------------.  .------------.
Cr       ldmx |            |  |            |
Cr            .------------.  .------------.
Cr
Cr            .------------.  .------------.
Cr       ldmx |            |  |            |
Cr            .------------.  .------------.
Cu Updates
Cu   28 Aug 16  Redesgined for greater generality
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldmx,ndgn,ld0,ndg0,nspc,offi
      double precision gii(ld0,nspc,ndg0,nspc,2),gin(ldmx,nspc,ndgn,nspc,2)
C ... Local parameters
      integer i,j,i1,j1

      if (mode == 0) call dpzero(gin,size(gin))
      if (mode /= 0) call dpzero(gii,size(gii))

C      call yprm('gin',2,gin,ldmx*nspc*ndgn*nspc,ldmx*nspc,ldmx*nspc,ndgn*nspc)
C      call yprm('gii',2,gii,ld0*nspc*ndg0*nspc,ld0*nspc,ld0*nspc,ndg0*nspc)

      do  j1 = 1, nspc
      do  i1 = 1, nspc

        if (mode == 0) then
          do  j = 1, ndg0
            do  i = 1, ld0
              gin(i,i1,j,j1,1) = gii(i,i1,j+offi,j1,1)
              gin(i,i1,j,j1,2) = gii(i,i1,j+offi,j1,2)
            enddo
          enddo
        else
          do  j = 1, ndg0
            do  i = 1, ld0
              gii(i,i1,j+offi,j1,1) = gin(i,i1,j,j1,1)
              gii(i,i1,j+offi,j1,2) = gin(i,i1,j,j1,2)
            enddo
          enddo
        endif

      enddo
      enddo

C      call yprm('gin',2,gin,ldmx*nspc*ndgn*nspc,ldmx*nspc,ldmx*nspc,ndgn*nspc)
C      call yprm('gii',2,gii,ld0*nspc*ndg0*nspc,ld0*nspc,ld0*nspc,ndg0*nspc)

      end
