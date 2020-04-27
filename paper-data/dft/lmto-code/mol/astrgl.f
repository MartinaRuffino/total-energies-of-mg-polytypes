      subroutine astrgl(nphi,lphi,nel,el,nlma,jb,pos,n0,nbas,ips,alat,
     .  cg,jcg,indxcg,cy,nhs,nhl,nstate,ioffb,ioffl,ixp,dxp,z,
     .  ceadd,bz,gbz,gdz,bzl,gbzl,gdzl)
C- Grad strux * evec for LMTO forces (vectorizes, linked basis)
C  All strx are expanded about site jb;
C  grad strx decomposed into ib for explicit decomposition of force.
      implicit none
      integer n0,nbas,nel,nhs,nhl,nphi(1),lphi(n0,1),jcg(1),indxcg(1),
     .  ips(1),ioffb(nbas,nel),ioffl(nbas+1),ixp(1),nlma,nstate
      double precision el(1),pos(3,1),dr(3),cg(1),cy(1),alat,e,
     .  gbz(nlma,nstate,nbas,3,nel),z(nhs,nstate),dxp(1),ceadd(25,5,1),
     .  gdz(nlma,nstate,3,nel),bz(nlma,nstate,nel),
     .  gbzl(nlma,nstate,nbas,3,nel),
     .  gdzl(nlma,nstate,3,nel),bzl(nlma,nstate,nel)
      integer ib,ie,je,is,jb,m,nlmb,nlmbx,nlmp,npow,nlmj,ih,nlmh,ll,
     .  iofz,ocf,oip,oikl,ojkl(0:20),ojj,ohl,ohd,os,osd,ogs,ogd,owk,
     .  nlk,lmxaj
      real w(1)
      common /w/ w
C     character *80 strn
C     integer i

C --- Setup ---
      nlk = 0
      if (nhl > 0) nlk = 1
      lmxaj = ll(nlma)
      call strxsu(nlma,nphi,lphi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,  nlma*nlmbx)
      call defrr(osd, nlma*nlmbx)
      call defrr(ogs, nlma*nlmbx*3)
      call defrr(ogd, nlma*nlmbx*3)
      call defrr(owk, nlma*nlmbx)
      call dpzero(bz, nlma*nstate*nel)
      call dpzero(gdz,nlma*nstate*nel*3)

C --- ib loop over function centers; expand around jb ---
      do  10  ib = 1, nbas
      is = ips(ib)
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,ib)-pos(m,jb))
C  ...  Strx for all energies, this pair
        nlmh = (ll(nlmp)+2)**2
        call defrr(ohl,nlmh*(nphi(is)+nlk))
        call defrr(ohd,nlmh*(nphi(is)+nlk))
        call rstr0(nphi(is)+nlk,lphi(1,is),el,nlmh,1,dr(1),dr(2),dr(3),
     .    lmxaj+1,0,w(ohl),w(ohd))
        do  12  ie = 1, nel
          nlmb = ioffb(ib+1,ie)-ioffb(ib,ie)
          if (nlmb <= 0) goto 12
          if (ixp(1) == 0 .and. ib == jb) then
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
            nlmj = ioffb(jb+1,ie)-ioffb(jb,ie)
            iofz = ioffb(ib,ie)+1
            call axtrg0(nlma,nlmb,nlmj,nhs,nbas,nstate,z(iofz,1),
     .        w(os),w(ogs),w(ogd),bz(1,1,ie),gbz(1,1,ib,1,ie),
     .        gdz(1,1,1,ie))

C       print *, 'z for ie,ib=', ie,ib,nlmb
C       call prmx(' ',z(iofz,1),nhs,nlmb,nstate)
C       strn = ' '
C       call awrit2('cp out.dat xx/z.%i%i',strn,len(strn),0,ie,ib)
C       call fsystm(strn,i)
C       print *, 'strux for jb,ie,ib=', jb,ie,ib
C       call prmx(' ',w(os),nlma,nlma,nlmb)
C       print *, 'grad strux for jb,ie,ib=', jb,ie,ib
C       call prmx(' ',w(ogs),nlma,nlma,nlmb)
C       print *, 'dot strux for jb,ie,ib=', jb,ie,ib
C       call prmx(' ',w(osd),nlma,nlma,nlmb)
C       print *, 'grad dot strux for jb,ie,ib=', jb,ie,ib
C       call prmx(' ',w(ogd),nlma,nlma,nlmb)
C       strn = ' '
C       call awrit3('cp out.dat xx/b.%i%i%i',strn,len(strn),0,jb,ie,ib)
C       call fsystm(strn,i)
C       print *, 'bz for jb,ie,ib=', jb,ie,ib
C       call prmx(' ',bz(1,1,ie),nlma,nlma,nstate)
C       strn = ' '
C       call awrit3('cp out.dat xx/bz.%i%i%i ',strn,len(strn),0,jb,ie,ib)
C       call fsystm(strn,i)
C       print *, 'gbz for ib,ie=', ib,ie
C       call prmx(' ',gbz(1,1,ib,1,ie),nlma,nlma,nstate)
C       strn = ' '
C       call awrit3('cp out.dat xx/xgbz.%i%i%i ',strn,80,0,jb,ie,ib)
C       call fsystm(strn,i)
C       print *, 'gdz for ib,ie=', ib,ie
C       call prmx(' ',gdz(1,1,1,ie),nlma,nlmj,nstate)
C       strn = ' '
C       call awrit3('cp out.dat xx/xgdz.%i%i%i ',strn,80,0,jb,ie,ib)
C       call fsystm(strn,i)

          endif
   12   continue

C ---   b*z, grad b*z, bd*z grad bd*z for linked basis ---
        if (nhl > 0) then
          ie = nel+1
          nlmbx = ioffl(ib+1)-ioffl(ib)
          if (nlmbx <= 0) goto 16
          if (ixp(1) == 0 .and. ib == jb) then
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmbx,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
            do  15  ie = 1, nel
              nlmb = ioffb(ib+1,ie)-ioffb(ib,ie)
              if (nlmb <= 0) goto 15
              if (ixp(1) == 0 .and. ib == jb) then
              else
                nlmj = 0
                do  17  je = 1, nel
   17           nlmj = max(nlmj,ioffb(jb+1,je)-ioffb(jb,je))
                iofz = ioffb(ib,ie)+1
                call axtrgl(nlma,nlmbx,nlmb,nlmj,nhs,nbas,nstate,
     .            z(iofz,1),w(os),w(ogs),w(ogd),ceadd(1,ie,is),w(owk),
     .            bzl(1,1,ie),gbzl(1,1,ib,1,ie),gdzl(1,1,1,ie))

C      print *, 'linked dot strux for jb,ie,ib=', jb,ie,ib
C      call prmx(' ',w(osd),nlma,nlma,nlmb)
C      print *, 'linked grad dot strux for jb,ie,ib=', jb,ie,ib
C      call prmx(' ',w(ogd),nlma,nlma,nlmb)
C      print *, 'linked strux for jb,ie,ib=', jb,ie,ib
C      call prmx(' ',w(os),nlma,nlma,nlmb)
C      print *, 'linked grad strux for jb,ie,ib=', jb,ie,ib
C      call prmx(' ',w(ogs),nlma,nlma,nlmb)
C      strn = ' '
C      call awrit3('cp out.dat xx/xgbl.%i%i%i',strn,len(strn),0,jb,ie,ib)
C      call fsystm(strn,i)
C      print *, 'bzl for jb,ie,ib=', jb,ie,ib
C      call prmx(' ',bzl(1,1,ie),nlma,nlma,nstate)
C      strn = ' '
C      call awrit3('cp out.dat xx/bzl.%i%i%i ',strn,len(strn),0,jb,ie,ib)
C      call fsystm(strn,i)
C      print *, 'xgbzl for jb,ie,ib=', jb,ie,ib
C      call prmx(' ',gbzl(1,1,ib,1,ie),nlma,nlma,nstate)
C      strn = ' '
C      call awrit3('cp out.dat xx/xgbzl.%i%i%i',strn,len(strn),0,jb,ie,ib)

C      call fsystm(strn,i)

              endif
   15       continue
          endif
   16     continue
        endif

        call rlse(ohl)
   10 continue


      call rlse(ocf)

      end
      subroutine axtrg0(nlma,nlmb,nlmj,nhs,
     .  nbas,nstate,z,s,gs,gd,bz,gbz,gdz)
C- Grad strux * evec for one block
      implicit none
      integer nlma,nlmb,nhs,nlmj,ix,nbas,nstate
      double precision s(nlma,nlmb),gs(nlma,nlmb,3),z(nhs,nstate),
     .  bz(nlma,nstate),gbz(nlma,nstate,nbas,3),
     .  gd(nlma,nlmb,3),gdz(nlma,nstate,3)

C --- Accumulate b * z from this block ---
      call dgemm('N','N',nlma,nstate,nlmb,1d0,s,nlma,z,nhs,1d0,bz,nlma)

C --- Grad b * z,  Grad b-dot z ---
      do  10  ix= 1, 3
        call dgemm('N','N',nlma,nstate,nlmb,1d0,gs(1,1,ix),nlma,
     .    z,nhs,0d0,gbz(1,1,1,ix),nlma)
        call dgemm('N','N',nlmj,nstate,nlmb,1d0,gd(1,1,ix),nlma,
     .    z,nhs,1d0,gdz(1,1,ix),nlma)
   10 continue

      end
      subroutine axtrgl(nlma,nlmbx,nlmb,nlmj,nhs,
     .  nbas,nstate,z,s,gs,gd,ceadd,wk,bz,gbz,gdz)
C- Grad strux * evec for one block, linking basis
      implicit none
      integer nlma,nlmbx,nlmb,nhs,nlmj,nbas,nstate,ix,i,j
      double precision s(nlma,nlmbx),gs(nlma,nlmbx,3),z(nhs,nstate),
     .  bz(nlma,nstate),gbz(nlma,nstate,nbas,3),ceadd(nlmb),
     .  gd(nlma,nlmbx,3),gdz(nlma,nstate,3),wk(nlma,nlmb)

C --- b * ceadd * z ---
      do  10  i = 1, nlma
      do  10  j = 1, nlmb
   10 wk(i,j) = s(i,j)*ceadd(j)
      call dgemm('N','N',nlma,nstate,nlmb,1d0,wk,nlma,z,nhs,1d0,bz,nlma)

C --- Grad b * ceadd * z,  Grad b-dot * ceadd * z ---
      do  20  ix= 1, 3
        do  24  i = 1, nlma
        do  24  j = 1, nlmb
   24   wk(i,j) = gs(i,j,ix)*ceadd(j)
        call dgemm('N','N',nlma,nstate,nlmb,1d0,wk,nlma,
     .    z,nhs,0d0,gbz(1,1,1,ix),nlma)
        do  28  i = 1, nlma
        do  28  j = 1, nlmb
   28   wk(i,j) = gd(i,j,ix)*ceadd(j)
        call dgemm('N','N',nlmj,nstate,nlmb,1d0,wk,nlma,
     .    z,nhs,1d0,gdz(1,1,ix),nlma)
   20 continue

      end
