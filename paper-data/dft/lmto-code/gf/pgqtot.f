      subroutine pgqtot(npl,nsp,glist,pgplp,dosi,ipc,z,qc,zval,qtot,sumev)
C- Evaluate and printout net charges and single-particle energies
C ----------------------------------------------------------------------
Ci Inputs
Ci   npl   :number of principal layers (pgfset.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   glist :a table of PL for which diagonal G.F. are made, with
Cl          the structure defined in pgglst.
Ci   pgplp :index and dimensioning information for crystal subblocks.
Ci          The meaning of pgplp depends on the context; see susite.f
Ci   dosi  :dos and integrated dos, resolved by layer.  For each layer:
Ci          1: dos at this e : -wtkp/pi Im G
Ci          2: nos at this e : -wtkp/pi Im (wz*G)
Ci          3: 1st energy mom : -wtkp/pi Im (wz*G z)
Ci          4: 2nd energy mom : -wtkp/pi Im (wz*G z**2)
Ci          5: projection of dos to fermi level (not calculated here)
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   z     :nuclear charge
Ci   qc    :core charges
Co Outputs
Co   qtot  :total charge
Co   sumev :sum of eigenvalues
Co   zval  :total core+nuclear charge
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npl,nsp,glist(-1:npl),pgplp(6,-1:npl),ipc(*)
      double precision dosi(5,2),qtot,sumev,z(*),qc(*),zval
C ... Local parameters
      integer ipr,lgunit,stdo,kpl,ipl,ib1,ib2,i,j,scrwid
      double precision wk(10)
      parameter(scrwid=80)
      character*6 strn

      call getpr(ipr)
      stdo = lgunit(1)
      if (ipr >= 25) write(stdo,370)
  370 format('     PL      D(Ef)      N(Ef)       E_band      2nd mom      Q-Z')

      strn = ' '
      if (nsp == 2) strn = 'spin 1'
      qtot = 0
      sumev = 0
      do  70  kpl = 1, glist(-1)
        ipl = glist(-1+kpl)
        call gtibpl(ipl,npl,pgplp,ib1,ib2)
        call gtqval(ib1,ib2,z,qc,ipc,zval)
        j = 5*nsp*(kpl-1)
        call dpscop(dosi,wk,5,1+j,1,1d0)
        if (wk(5) /= 0) wk(1) = wk(5)
        if (ipr >= 25) then
          if (glist(-1) /= 1 .or. ipl /= 1) then
            write(stdo,371) ipl, (wk(i), i=1,4), wk(2)-zval/nsp
          else
            write(stdo,372) strn, (wk(i), i=1,4), wk(2)-zval/nsp
          endif
        endif
        qtot = qtot + wk(2) - zval
        sumev = sumev + wk(3)
        if (nsp == 2) then
          call dpscop(dosi,wk,5,6+j,6,1d0)
          if (wk(10) /= 0) wk(6) = wk(10)
          qtot = qtot + wk(7)
          sumev = sumev + wk(8)
          if (ipr >= 25) then
            write(stdo,372) 'spin 2', (wk(i+5),i=1,4),wk(7)-zval/nsp
            write(stdo,372) ' total', (wk(i)+wk(i+5), i=1,4),wk(2)+wk(7)-zval
            call info2(10,0,0,' N(up)-N(dn)%7f %;11,6D',wk(2)-wk(7),2)
          endif
        endif
C   ... Update since corrected DOS might be available
        call dpscop(wk,dosi,5*nsp,1,1+j,1d0)
  371   format(i7,5f12.6)
  372   format(1x,a,5f12.6)
   70 continue
      if (ipr >= 25)
     .  call awrit1('%9fdeviation from charge neutrality: %d',' ',scrwid,stdo,qtot)

      end
