      subroutine xxxdif(emin,emax,npts,ndos,mode,dos)
C- Differentiate Number of States to make DOS
C ----------------------------------------------------------------------
Ci Inputs
Ci   emin  :lower energy of DOS
Ci   emax  :upper enenergy of DOS
Ci   npts  :number of DOS tabulation points (input; for sampling only)
Ci   ndos  :number of energy mesh points
Ci   mode  :not used
Cio Inputs/Outputs
Cio  dos   :On input, number-of-states
Cio         :On output, density-of-states
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Sep 09 xxxdif no longer requires extra column for workspace
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npts,ndos,mode
      double precision dos(npts,ndos+1),emin,emax
C ... Local parameters
      integer ip,idos
      double precision bin,wk(npts)

      bin = 2*(emax-emin)/dble(npts-1)
      do  idos = 1, ndos
        do  ip = 2, npts-1
          wk(ip) = (dos(ip+1,idos)-dos(ip-1,idos))/bin
        enddo
        wk(1) = wk(2)
        wk(npts) = wk(npts-1)
        call dcopy(npts,wk,1,dos(1,idos),1)
      enddo
      end
