      subroutine shftpq(lrel,nc,nrclas,nsp,nl,lmx,rmax,avw,pp,amom,
     .           amomr,idmod,swmod,pmin,pnu,qnu,qnur,sumevm)
C- Shift p and q according to idmod
C ----------------------------------------------------------------
Ci Inputs
Ci   lrel  :2 for Dirac moments qnur
Ci   nc    :number of classes
Ci   nrclas:nrclas(i) = number of atoms in the ith class
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   pp    :potential parameters (atomsr.f)
Ci   amom  :input values of moments qnu
Ci   amomr :input values of moments qnur
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
Ci   swmod :see Remarks
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   qnur  :relativistic moments of the sphere charges
Ci   nc,nsp,lmx,nl,rmax,pp,amom,idmod,swmod
Co Outputs
Co   pnu,qnu:pnu possibly floated, and qnu computed from amom but
Co          :1st and second energy moments shifted according to
Co          :shifted pnu
Co   sumevm :
Cv Verbosity
Cv   >=30: print out moments and shift in enu
Cr Remarks
Cr   when swmod=.true. prompt for each idmod=2 for shift in
Cr   amom and qnu can point to the same address space
Cr   command argument -enu=val sets each idmod=2 and each enu to val
Cu Updates
Cu   01 Aug 16 Some changes for Dirac case
Cu   15 Apr 15 (Vishina) New lrel=2
Cu   07 Apr 15 New idmod=4
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical swmod
      integer nl,nsp,nc,idmod(0:nl-1,nc),lmx(nc),nrclas(nc),lrel
      double precision avw,sumevm
      double precision amom(3,0:nl-1,nsp,nc),pp(6,nl,nsp,nc),rmax(nc)
      double precision pmin(0:*),pnu(0:nl-1,nsp,nc),qnu(3,0:nl-1,nsp,nc)
      double precision qnur(4,0:nl-1,2*nl,2,2,nc),amomr(4,0:nl-1,2*nl,2,2,nc)
C ... Local parameters
      integer ic,isp,l,i,j,k,idm,imu,ii
      double precision pold,pminl,ql(3),eb,pli(2),ebi(2),xx(1),mu
      double precision qnurs(4,0:nl-1,2*nl,2,2)
      procedure(integer) :: iprint,lgunit
      procedure(logical) :: cmdopt,a2bin
      character*72 outs

      if (iprint() >= 32) then
        do  j = 1, 2
          write (lgunit(j),1)
    1     format(/' CLASS L    Q0',9x,'Q1',9x,'Q2',9x,'EB',9x,'POLD',7x,'PNU')

        enddo
      endif

      do  ic = 1, nc
      do  isp = 1, nsp
      if (lrel == 2 .and. isp == 1)
     .    call dcopy(size(qnurs),amomr(1,0,1,1,1,ic),1,qnurs,1)

      do  l = 0, lmx(ic)
        idm = idmod(l,ic)
        if (idm == 4 .and. nsp == 1) idm = 1
        if (idm == 4) cycle
        if (iprint() >= 32) then
          do  j = 1, 2
C           write(lgunit(j),331) ic,l,(amom(i,l,isp,ic), i=1,3)
            call awrit4('%,4i%,4i%;10,6D%2;11,6D',' ',80,lgunit(j),ic,l,amom(1,l,isp,ic),amom(2,l,isp,ic))
          enddo
        endif
        eb = 0
        if (cmdopt('-enu=',5,0,outs)) then
          j = 5
          if (.not. a2bin(outs,eb,4,0,' ',j,72)) call rx('shftpq: bad value -enu')
          eb = eb - pp(1,l+1,isp,ic)
          idm = 2
          idmod(l,ic) = 2 + (idmod(l,ic)-mod(idmod(l,ic),10))
        endif
        if (swmod .and. idm == 2) call query('eb=',4,eb)
        pold = pnu(l,isp,ic)
        pminl = pmin(l)
C       print *, 'call enutod',l,ic,isp

        if (lrel == 2 .and. isp == 1 .and. idm /= 1 .and. idm /= 2) then
          call enutod(.true.,nl,l,ic,rmax(ic),avw,pp(1,l+1,isp,ic),amom(1,l,isp,ic),
     .      amomr(1,0,1,1,1,ic),0d0,idm,pminl,pnu(l,isp,ic),qnu(1,l,isp,ic),qnur(1,0,1,1,1,ic),eb)
        else
          call enutod(.false.,nl,l,ic,rmax(ic),avw,pp(1,l+1,isp,ic),amom(1,l,isp,ic),
     .    xx,0d0,idm,pminl,pnu(l,isp,ic),qnu(1,l,isp,ic),xx,eb)
        endif
C       (I think ...)
        if (eb /= 0 .and. idm == 2) pp(1,l+1,isp,ic) = pp(1,l+1,isp,ic) + eb
        if (iprint() >= 32 .and. eb /= 0) then
          do  j = 1, 2
            if (pminl == pnu(l,isp,ic)) then
C             write(lgunit(j),332) (qnu(i,l,isp,ic), i=1,3), eb, pold, pnu(l,isp,ic), ' *'
              call awrit5('%8f%;10,6D%2;11,6D%;11,6D%;11,6D%;11,6D *',' ',80,lgunit(j),
     .          amom(1,l,isp,ic),amom(2,l,isp,ic),eb,pold,pnu(l,isp,ic))
            else
C             write(lgunit(j),332) (qnu(i,l,isp,ic), i=1,3), eb, pold, pnu(l,isp,ic)
              call awrit5('%8f%;10,6D%2;11,6D%;11,6D%;11,6D%;11,6D',' ',80,lgunit(j),
     .          amom(1,l,isp,ic),amom(2,l,isp,ic),eb,pold,pnu(l,isp,ic))
            endif
          enddo
        endif
      enddo
      enddo

      if (lrel == 2 .and. iprint() >= 32) then
        j = 1
        write (lgunit(j),2)
    2   format('  l  mu   ms1 ms2',5x,'q0',11x,'q1',11x,'q2',11x,'enu')
        do  l = 0, lmx(ic)
          do  imu = 1, 2*(l+1)
            mu = dble(imu-l) - 1.5d0
            do  i = 1, 2
              do  k = 1, 2
                if (qnurs(1,l,imu,i,k) /= 0) then
                write(lgunit(j),23) l,mu,i,k,(qnurs(ii,l,imu,i,k), ii=1,4)
                write(lgunit(j),23) l,mu,i,k,(qnur(ii,l,imu,i,k,ic), ii=1,4)
   23           format(i3,f5.1,2i4,4f13.8)
                endif
              enddo
            enddo
          enddo
          write(lgunit(j),24) l,(sum(qnur(ii,l,:,1,1,ic)+qnur(ii,l,:,2,2,ic)), ii=1,3)
   24     format(i3,' l sum mu',4x,3f13.8)
        enddo
        write(lgunit(j),25)  (sum(qnur(ii,:,:,1,1,ic)+qnur(ii,:,:,2,2,ic)), ii=1,3)
   25   format(4x,'sum l+mu',4x,3f13.8)
        write(lgunit(j),26)  (sum(amom(ii,:,1:nsp,ic)), ii=1,3)
   26   format(4x,'scalar rel Q',3f13.8/)
      endif

      enddo

C ... Special treatment for idmod=4
      do  ic = 1, nc
      do  l = 0, lmx(ic)
        idm = idmod(l,ic)
        if (idm /= 4 .or. nsp == 1) cycle
        if (iprint() >= 32) then
          do  j = 1, 2
            write(lgunit(j),331) ic,l,(amom(i,l,1,ic), i=1,3)
            write(lgunit(j),331) ic,l,(amom(i,l,2,ic), i=1,3)
          enddo
        endif
        pminl = pmin(l)
        ebi = 0
        pli(1) = pnu(l,1,ic)
C       Unconstrained float spin1
        call enutod(.false.,nl,l,ic,rmax(ic),avw,pp(1,l+1,1,ic),amom(1,l,1,ic),
     .    xx,0d0,idm-3,pminl,pli(1),ql,xx,ebi(1))
C       Unconstrained float spin2
        pli(2) = pnu(l,2,ic)
        call enutod(.false.,nl,l,ic,rmax(ic),avw,pp(1,l+1,2,ic),amom(1,l,2,ic),
     .    xx,0d0,idm-3,pminl,pli(2),ql,xx,ebi(2))
        ebi(1) = ebi(1) + (pp(1,l+1,2,ic)-pp(1,l+1,1,ic))/2
C       Constrain float spin1 to put enu near average of spin 1, spin2
        pli(1) = pnu(l,1,ic)
        call enutod(.false.,nl,l,ic,rmax(ic),avw,pp(1,l+1,1,ic),amom(1,l,1,ic),
     .    xx,0d0,2,pminl,pnu(l,1,ic),qnu(1,l,1,ic),xx,ebi(1))
C       Constrain float spin2 to put enu near average of spin 1, spin2
        ebi(2) = ebi(2) - (pp(1,l+1,2,ic)-pp(1,l+1,1,ic))/2
        pli(2) = pnu(l,2,ic)
        call enutod(.false.,nl,l,ic,rmax(ic),avw,pp(1,l+1,2,ic),amom(1,l,2,ic),
     .    xx,0d0,2,pminl,pnu(l,2,ic),qnu(1,l,2,ic),xx,ebi(2))

        if (iprint() >= 32) then
          do  isp = 1, 2
            pold =pli(isp)
            eb = ebi(isp)
            do  j = 1, 2
            if (pminl == pnu(l,isp,ic)) then
              write(lgunit(j),332) (qnu(i,l,isp,ic), i=1,3), eb,
     .                              pold, pnu(l,isp,ic), ' *'
            else
              write(lgunit(j),332) (qnu(i,l,isp,ic), i=1,3), eb,
     .                              pold, pnu(l,isp,ic)
            endif
            enddo
          enddo
        endif
      enddo
      enddo

      sumevm = 0d0
      do  ic = 1, nc
        do  isp = 1, nsp
          do  l = 0, lmx(ic)
            sumevm = sumevm + (qnu(2,l,isp,ic) +
     .                         qnu(1,l,isp,ic)*pp(1,l+1,isp,ic))*nrclas(ic)
          enddo
        enddo
      enddo
  331 format(2i4,f10.6,5f11.6)
  332 format(8x,f10.6,5f11.6:a)
      end
