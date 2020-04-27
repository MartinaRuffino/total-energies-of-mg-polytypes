      subroutine mkstrx(ldip,nbas,nlmq1,nlmq,bas,ipc,lmxl,awld,alat,
     .                  vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy,pv,MOL,
     .                  strx,dstrx)
C- Make lookup table of strux and radial derivatives
C ----------------------------------------------------------------------
Ci Inputs:
Ci  ldip  : 3 include 'spherical' dipole correction to Ewald
Ci        : 2 include 'slab' dipole correction to Ewald
Ci        : any other number - skip the dipole correction
Ci   nbas,bas,awld,alat,vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy
Ci   nlmq1: leading dimension of strx, nlmq1=(ll(nlmq)+1)**2
Ci   nlmq : max L-cutoff for multipoles, leading dimensions of dstrx,
Ci           second dimension of strx
Ci   lmxl : l-cutoff for multipoles for each class of species
Ci   pv   : true if radial derivatives required
Ci   MOL  : true for molecule (cluster) branch
Co Outputs:
Co  strx  : coefficients B_LL'(R'-R) (structure constants)
Co  dstrx : derivatives of structure constants (x dB_LL'(x)/dx) at x = |R'-R|
Co          if pv = F dstrx is not touched (see hstra.f)
Cr Remarks
Cr   Instead of making B on each sc iteration, the program prepares the whole B
Cr   in the beginning of the self-consistency cycle to be used as a lookup table.
Cr   This results in faster calculation at an expence of memory.
Cr
Cr   Calling mkstrx is set as the default option. To call instead hstr on each
Cr   iteration (and therefore to reduce memory but increase computational time),
Cr   use switch --sfly
Cr
Cr   Efficiency issues: there are two symmetry properties of zero energy
Cr   structure constants that are used in mkstrx:
Cr     B_LL'(R'-R) = B_L'L(R-R')                (1)
Cr     B_LL'(R'-R) = (-1)^(l+l') B_L'L(R'-R)    (2)
Cr   same properties hold for the radial derivative of B.
Cr
Cr   According to (1), mkstrx calculates strx(:,:,ib,jb) only for the
Cr   lower triangle {ib = 1:nbas, jb = 1:ib}
Cr   property (2) is used in hstra in the same way, ie strx(ilm,jlm,:,:)
Cr   is sought only for the triangle {ilm = 1:Lmax, jlm=1:ilm}.
Cr   hstra should therefore be called with Lmax >= Lmax', otherwise it stops.
Cr
Cr   mkstrx fills only lower triangles of strx and dstrx in the above sense;
Cr   properties (1) and (2) are used later, when strx and dstrx are copied into
Cr   'local' arrays fstrx and fdstrx in tbesel.f We decided to keep it this way
Cr   because in future we plan to reduce the size of strx and dstrx by cutting of
Cr   the unused parts of arrays strx and dstrx.
Cr
Cr   hstra vs hstr: a call to hstr (a routine that actually makes strux) is
Cr   replaced with a call to hstra, an 'optimised' version of hstr. The difference
Cr   between hstr and hstra is that hstra computes only the lower triangle (see above)
Cr   of strux but is restricted to Lmax >= Lmax'. hstr does not have any restrictions
Cr   and is called if the --sfly switch is on.
Cr
Cb Bugs
Cb   mkstrx is not written in the most efficient way and could be further refined wrt to
Cb   amount of calculated elements of strux. At the moment mkstrx is a result of certain
Cb   trade off between performance and clarity of the program.
Cb
Cu Updates
Cu    05 Mar 2010 (SL)  optimization (make strux up to Lmax and for triangle
Cu                      jb<=ib only)
Cu    19 Jan 2010 (ATP) first written
C ----------------------------------------------------------------------
      use mod_tbfuns
      implicit none
C Passed Parameters
      logical, intent(in) :: pv,MOL
      integer, intent(in) :: nbas,nlmq,nlmq1,ipc(nbas),lmxl(*)
      integer, intent(in) :: ldip,indxcg(*),jcg(*),nkd,nkg
      double precision, intent(in)  :: bas(3,nbas),alat,awld,vol,
     .                                 dlat(*),glat(*),cy(*),cg(*)
      double precision, intent(out) :: strx(nlmq1,nlmq,nbas,nbas),
     .                                 dstrx(nlmq,nlmq,nbas,nbas)
C Local Variables
      integer ib,jb,lmax,lmxf,lmxst,nlx,nlf,i1mach,iprint
      integer li,li1,lj
      double precision tau(3),taua(3),hl(100),bl(100)

      call tcn('mkstrx')

      do  ib = 1, nbas
        li = lmxl(ipc(ib))
c       nlmi = (li+1)**2
C...  nlmq1 == nlmq if no force or pressure are required.
C     Otherwise need to go up to lmax+1 in the first index of B_LL'
        if (nlmq1 > nlmq) then
          li1 = li + 1
        else
          li1 = li
        endif
        do  jb = 1, ib
          lj = lmxl(ipc(jb))
C...  The lines below are because we shall be filling the jb > ib triangle
C     using symmetry properties of B (see tbesel.f) AND because we need l+1
C     for forces. No complications if forces are not required.
          if (nlmq1 > nlmq) then
            lmax = max(li,lj)
            lmxf = lmax+1
            nlx = (lmax+1)**2
            nlf = (lmax+2)**2
            lmxst = 2*lmax+1
          else
            nlx = (lj+1)**2
            nlf = (li+1)**2
            nlf = max(nlf,nlx)
            lmxst = li + lj
          endif
c         call dmadd(bas(1,jb),3,1,1d0,bas(1,ib),3,1,-1d0,tau,3,1,3,1)
          tau(1:3) = bas(1:3,jb)-bas(1:3,ib)
          if (.not. MOL) then
            call shortn(tau,tau,dlat,nkd)
            call rcnsl0(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl)
c           le = 0
          else
            taua(1:3) = tau(1:3) * alat
            call soldhj(taua,0d0,0,lmxst,hl,bl,cy)
c           le = 1
          endif
          if (iprint() > 60)
     .      call awrit2('tau=%3:1d, rstrx=%100:1d',' ',1028,i1mach(2),
     .                   tau,hl)
          call hstra(MOL,pv,ldip,strx(1,1,ib,jb),dstrx(1,1,ib,jb),
     .               nlf,nlx,nlmq1,nlmq,hl,cg,indxcg,jcg,vol)
          call tbshfl(0,nlmq1,nlf,nlx,strx(1,1,ib,jb))
          call strfac(0,nlmq1,nlf,nlx,strx(1,1,ib,jb))
          if (pv) then
            call tbshfl(0,nlmq,nlx,nlx,dstrx(1,1,ib,jb))
            call strfac(0,nlmq,nlx,nlx,dstrx(1,1,ib,jb))
          endif
        enddo
      enddo

      call tcx('mkstrx')
      end subroutine mkstrx
