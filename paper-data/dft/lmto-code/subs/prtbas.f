      subroutine prtbas(nbas,ipclas,clabnw,ncell,nclnw,mapc,smap,iwk,
     .  tau,z,avw,rmax,ves,switch,nl,pl,ql)
C- Output basis vectors of supercell
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas,ipclas,clabnw,tau
Ci   nbas:   number of atoms in (old) basis
Ci   ipclas: table of class pointers in (old) basis
Ci   mapc:   table of mappings from old to new class pointers
Ci   iopt:   0, write classes as found in original
Ci           1, make classes unique
Ci   ncell:  number of cells in supercell
Ci   nclnw:  number of (new) classes
Ci   clabnw: class labels in supercell
Co Outputs
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nsp
      parameter (nsp=1)
      integer nbas,ipclas(nbas),ncell,nclnw,mapc(ncell,4),icold,nl,
     .  smap(ncell*nbas),iwk(8)
      logical switch(14)
      character*(8) clabnw(1)
      double precision tau(3,nbas,ncell), z(*), rmax(*), avw
      double precision pl(nl,nsp,1),ql(3,nl,nsp,1),ves(*)
C Local parameters
      integer i,k,fopn,ifi,lgunit,ic,ii,m,icell,isp,l,ix,jc,ifi2
      character*13 strn

      call dfclos(fopn('LOG'))
      ifi = fopn('LOG')

      ifi = lgunit(1)
      ifi2 = lgunit(2)
      if (smap(1) >= 0) call iinit(iwk,ncell*nbas)

C --- Output one line of STRUC in ctrl format ---
      call awrit3('STRUC   NBAS=%i  NCLASS=%i  NL=%i',
     .  ' ',80,ifi,nbas*ncell,nclnw,nl)
      call awrit3('STRUC   NBAS=%i  NCLASS=%i  NL=%i',
     .  ' ',80,ifi2,nbas*ncell,nclnw,nl)
C      write(ifi,12) nbas*ncell, nclnw, nl
C  12  format('STRUC   NBAS=',i3,'  NCLASS=',i2,'  NL=',i1)

C --- Output CLASS in ctrl format ---
      m = 0
      strn = 'CLASS   ATOM='
      do  jc = 1, nclnw
        ic = jc
C ... Get the next unwritten map sorted by tau
        if (smap(1) >= 0) then
    2     continue
          m = m+1
          ix = smap(m)+1
          i = mod(ix-1,nbas)+1
          ii = (ix-1)/nbas+1
          ic = mapc(ii,ipclas(i))
C         print *, 'm,ix,i,ii,ic=',m,ix,i,ii,ic
          if (iwk(ic) /= 0) goto 2
          iwk(ic) = 1
        endif
C ... Find icold for which ic = mapc(icell,icold)
        icold = 0
    5   continue
        icold = icold+1
        do  icell = 1, ncell
          if (mapc(icell,icold) == ic) goto 10
        enddo
        if (icold >= 1000) stop 'prtbas: bad mapc'
        goto 5
   10   continue
C          write(ifi,21) strn,clabnw(ic),nint(z(icold)),rmax(icold),
C     .      (idxdn(k,icold), k=1,nl)
C   21     format(A13,A4,'  Z=',i2,'  R/W=',f7.5,'  IDXDN=',5(i2:))
        call awrit3(strn//clabnw(ic)//' Z=%d  %29pR=%1;6d*%1;6d',
     .    ' ',80,ifi,z(icold),avw,rmax(icold))
        call awrit3(strn//clabnw(ic)//' Z=%d  %29pR=%1;6d*%1;6d',
     .    ' ',80,ifi2,z(icold),avw,rmax(icold))
        strn = '        ATOM='
      enddo

C --- Output SITE in ctrl format ---
      do  m = 1, ncell*nbas
        ix = m
        if (smap(1) >= 0) ix = smap(ix)+1
        i = mod(ix-1,nbas)+1
        ii = (ix-1)/nbas+1
        ic = mapc(ii,ipclas(i))
        if (ix == 1) then
          write (ifi,3) clabnw(ic),(tau(k,i,ii),k=1,3)
          write (ifi2,3) clabnw(ic),(tau(k,i,ii),k=1,3)
        endif
        if (ix /= 1) then
          write (ifi,4) clabnw(ic),(tau(k,i,ii),k=1,3)
          write (ifi2,4) clabnw(ic),(tau(k,i,ii),k=1,3)
        endif
      enddo

C --- Output START P,Q in ctrl format ---
      do  ic = 1, nclnw
C old class is icold for which ic = mapc(icell,icold)
        icold = 0
   15   continue
        icold = icold+1
        do  icell = 1, ncell
          if (mapc(icell,icold) == ic) goto 20
        enddo
        if (icold >= 1000) stop 'prtbas: bad mapc'
        goto 15
   20   continue

c          print *, '**', ic, icold, ncell
        write (ifi,6) clabnw(ic),((pl(l,isp,icold),l=1,nl),isp=1,nsp)
        write (ifi2,6) clabnw(ic),((pl(l,isp,icold),l=1,nl),isp=1,nsp)
        write (ifi,7) (((ql(i,l,isp,icold),i=1,3),l=1,nl),isp=1,nsp)
        write (ifi2,7) (((ql(i,l,isp,icold),i=1,3),l=1,nl),isp=1,nsp)
        if (switch(14)) then
          write (ifi,8) ves(icold)
          write (ifi2,8) ves(icold)
        endif
      enddo

C --- Output in lc format ---
      m = 0
      do  i = 1, nbas
        do  ii = 1, ncell
          m = m+1
          ic = mapc(ii,ipclas(i))
c         write(ifi,120) m, (tau(k,i,ii), k=1,3), clabnw(ic)
    1     format(' ATOM ',i3,2x,3F12.7,'   <',a4,'>')
        enddo
      enddo
    3 format('SITE    ATOM=',a4,'  POS=',3F12.7)
    4 format('        ATOM=',a4,'  POS=',3F12.7)
    6 format('          ATOM=',a4,'  P=',5(3F11.7:,/23x))
    7 format(19x,'  Q=',5(3F11.7:,/23x))
    8 format(19x,'  V=',f11.7)
      end
