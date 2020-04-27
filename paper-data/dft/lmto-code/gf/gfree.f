      subroutine gfree(mode,nlmi,nlmj,ndim,avw,v,rmaxi,rmaxj,z,indxcg,
     .  jcg,cg,cy,dpfunf,dpfuns,ddpfun,gij,gi2)
C- ASA Greens's function for free electrons
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0 return gij in complex*16 format
Ci          1 return gij with real and imaginary parts separated
Ci   nlmi  :dimension of G inside sphere Ri containing field point
Ci   nlmj  :dimension of G inside sphere Rj containing source point
Ci   ndim  :leading dimension of gij
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   v     :connecting vector Rj-Ri, in a.u.
Ci   rmaxi :augmentation radius for site Ri, in a.u.
Ci   rmaxj :augmentation radius for site Rj, in a.u.
Ci   z     :complex energy
Ci   indxcg,jcg:indices for Clebsch Gordan coefficients
Ci   cg    :Clebsch Gordan coefficients
Ci   cy    :Normalization constants for spherical harmonics
Ci   dpfunf,dpfuns:
Ci          scaling to convert gfree to some screened representation.
Ci          usu. sqrt(P-dot) or ratio of sqrt(P-dot), viz:
Ci          to Gamma, make dpfun by calling mkpotj with iopt=130
Ci          to alpha, make dpfun by calling mkpotj with iopt=40
Ci          g^alpha_ij = dpfun_i g^bare_ij dpfun_j for i<>j
Ci          dpfunf,dpfuns = P for field and source points, respectively
Ci          To skip scaling and return bare g, set dpfun(1) = 0
Ci   ddpfun:scaling to convert diagonal gfree to a screened repsn.
Ci          usu sqrt P-dot-dot, bare rep a difference in two, viz:
Ci          to Gamma, make dpfun by calling mkpotj with iopt=120
Ci          to alpha, make dpfun by calling mkpotj with iopt=50
Co Outputs
Co   gij   :Free electron Green's function.
Co          For bare g, see Appendix, see PRB 27, 7144, 1983, Eqn. B3.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      integer mode,nlmi,nlmj,ndim,indxcg(*),jcg(*)
      double precision v(3)
      double precision cg(*),cy(*),avw,rmaxi,rmaxj
      double precision gij(ndim,*),gi2(ndim,*)
      double complex z,dpfunf(ndim),dpfuns(ndim),ddpfun(ndim)
C Local variables
      logical lscale
      integer ll,ilm,jlm,li,lj,l,m,l1,l2,m1,m2,nds,mode0
      parameter (nds=49)
      double precision vv,ddot
      double complex s(nds,nds),gl,g
      double complex wkki(0:20),wkji(0:20),wjki(0:20),wjji(0:20)
      double complex wkkj(0:20),wkjj(0:20),wjkj(0:20),wjjj(0:20)

      mode0 = mod(mode,10)
      if (nlmi == 0 .or. nlmj == 0) return
      lscale = abs(dpfunf(1)) /= 0d0
      vv = dsqrt(ddot(3,v,1,v,1))
      li = ll(nlmi)
      lj = ll(nlmj)

C --- On-site g ---
      if (vv == 0d0) then
        call wrnhjz(z,(0d0,0d0),rmaxi,max(li,lj),avw,wkki,wkji,wjki,
     .    wjji,10)
        if (mode0 == 0) then
          do  20  jlm = 1, nlmj
          do  20  ilm = 1, nlmi
            gij(2*ilm-1,2*jlm-1) = 0
            gij(2*ilm,2*jlm-1) = 0
   20     continue
        else
          do  21  jlm = 1, nlmj
          do  21  ilm = 1, nlmi
            gij(ilm,jlm) = 0
            gi2(ilm,jlm) = 0
   21     continue
        endif

        ilm = 0
        do  22  l = 0, min(li,lj)
          gl = -4/avw**2*wjji(l)*wkji(l)
          g = gl
          do  22  m = -l, l
            ilm = ilm+1

            if (lscale) g =  dpfunf(ilm)*gl*dpfuns(ilm) + ddpfun(ilm)

            if (mode0 == 0) then
              gij(2*ilm-1,2*ilm-1) = dble(g)
              gij(2*ilm,2*ilm-1) = dimag(g)
            else
              gij(ilm,ilm) = dble(g)
              gi2(ilm,ilm) = dimag(g)
            endif

   22   continue

C --- Off=site g ---
      else

        call wrnhjz(z,(0d0,0d0),rmaxi,li,avw,wkki,wkji,wjki,wjji,10)
        call wrnhjz(z,(0d0,0d0),rmaxj,lj,avw,wkkj,wkjj,wjkj,wjjj,10)
        call mstrxz(z*avw**2,v,nlmi,nlmj,nds,cg,indxcg,jcg,cy,11,s,s)

        jlm = 0
        do  30  l2 = 0, lj
        do  30  m2 = -l2, l2
          jlm = jlm+1

          ilm = 0
          do  32  l1 = 0, li
          do  32  m1 = -l1, l1
            ilm = ilm+1

            g = -4/avw**2*wjji(l1)*s(ilm,jlm)*wjjj(l2)
            if (lscale) g = dpfunf(ilm)*g*dpfuns(jlm)

            if (mode0 == 0) then
              gij(2*ilm-1,2*jlm-1) = dble(g)
              gij(2*ilm,2*jlm-1)   = dimag(g)
            else
              gij(ilm,jlm) = dble(g)
              gi2(ilm,jlm) = dimag(g)
            endif
   32     continue
   30   continue

      endif

C      call yprm('gfree',3,gij,ndim*ndim,ndim,nlmi,nlmj)
C      call yprm('gfree-real',1,gij,ndim*ndim,ndim,nlmi,nlmj)
C      call yprm('gfree-imag',1,gi2,ndim*ndim,ndim,nlmi,nlmj)

      end
C      subroutine chkgfj(n1,n2,n3,nsp,ib,jb,idim,jdim,offi,offj,gij,alat,
C     .  avw,rmaxi,rmaxj,z,gv,kv,ng,tau,npf,dpfun,ddpfun,
C     .  indxcg,jcg,cg,cy)
CC- Check gfree
CC ----------------------------------------------------------------------
CCi Inputs
CCi   z         complex energy
CCi   n1,n2,n3  QP mesh
CCi   idim,jdim dimensions of G(i,j) subblock
CCi   ldiag     T, gji and gij are identical
CCi   dpfun     sqrt P-dot, bare rep: call to mkpotj with 130
CCi   ddpfun    sqrt P-dot-dot, bare rep: call to mkpotj with 120
CC ----------------------------------------------------------------------
C      implicit none
CC Passed variables
C      integer idim,jdim,offi,offj,ib,jb,n1,n2,n3,nsp,ng,npf,indxcg(*),
C     .  jcg(*),kv(ng,3)
C      double precision cg(*),cy(*),alat,avw,rmaxi,rmaxj
C      double precision tau(3),gv(ng,3)
C      double complex gij(n1,n2,n3,idim,jdim,nsp),z
C      double complex dpfun(npf),ddpfun(npf)
CC Local variables
C      integer id,jd,isp,iset,ig,ndim,lmax,ll,ilm,jlm,l,m,i1,i2,i3,i0,k
C      integer ifi,fopna
C      parameter (ndim=49)
C      double precision v(3)
C      double complex g,gf(ndim,ndim),gasa(ndim,ndim)
C      character*20 strn
C
C      print *, 'ib,jb,energy=',ib,jb,z
C
CC --- Loop over orbital pairs ---
C      iset= 0
C      call fftz3s(n1,n2,n3,n1,n2,n3,iset)
C      do  10  isp = 1, nsp
C      do  10  id = 1, idim
C      do  10  jd = 1, jdim
C
Cc        call zprm3(' gij %i',id+100+jd,gij(1,1,1,id,jd,isp),n1,n2,n3)
C
C        call fftz3(gij(1,1,1,id,jd,isp),n1,n2,n3,n1,n2,n3,1,iset,-1)
C
Cc        call zprm3(' gij %i',id+100+jd,gij(1,1,1,id,jd,isp),n1,n2,n3)
C
C   10 continue
C
CC ... Check on-site G
C      if (ib == jb) then
C      lmax = ll(idim)
C      call dpzero(v,3)
C      call gfree(0,idim,idim,ndim,avw,v,rmaxi,rmaxj,z,indxcg,jcg,cg,
C     .  cy,dpfun(1+offi),dpfun(1+offi),ddpfun(1+offi),gf,gf)
C      ilm = 0
C      do  20  l = 0, lmax
C        do  22  m = -l, l
C          ilm = ilm+1
C          k = ilm+offi
CC     ... Convert g(bare) to G(gamma)
Cc         g = dpfun(k)*gf(ilm,ilm)*dpfun(k) + ddpfun(k)
C          g = gf(ilm,ilm)
C
C          print 333, ilm,g,gij(1,1,1,ilm,ilm,1),gij(1,1,1,ilm,ilm,1)/g
C  333     format(i4,6f12.6)
C   22   continue
C   20 continue
C       call yprm('diag g',3,gf,ndim*ndim,ndim,idim,jdim)
C      endif
C
CC ... Loop through list of R vectors
C      i0 = 1
C      if (ib == jb) i0 = 2
C      do  50  ig = i0, ng, 1
C        i1 = kv(ig,1)
C        i2 = kv(ig,2)
C        i3 = kv(ig,3)
C        v(1) = -alat/avw*(gv(ig,1) - tau(1))
C        v(2) = -alat/avw*(gv(ig,2) - tau(2))
C        v(3) = -alat/avw*(gv(ig,3) - tau(3))
C        print 335, v(1)*avw/alat,v(2)*avw/alat,v(3)*avw/alat
C  335   format(' connecting vector',3f12.6)
C
C        call gfree(0,idim,jdim,ndim,avw,v,rmaxi,rmaxj,z,indxcg,jcg,cg,
C     .    cy,dpfun(1+offi),dpfun(1+offj),ddpfun(1+offi),gf,gf)
C
C        do  52  ilm = 1, idim
C        do  52  jlm = 1, jdim
CC     ... Convert g(bare) to G(gamma)
Cc         g = dpfun(ilm+offi)*gf(ilm,jlm)*dpfun(jlm+offj)
C          g = gf(ilm,jlm)
C          if (abs(g) > 1d-6) then
C          print 334, ilm,jlm, g, gij(i1,i2,i3,ilm,jlm,1),
C     .        (g/gij(i1,i2,i3,ilm,jlm,1))
C          gasa(ilm,jlm) = gij(i1,i2,i3,ilm,jlm,1)
C          endif
C  334     format(2i4,6f14.8)
C   52   continue
C
C        ifi = fopna('out1',-1,0)
C        call ywrm(0,' ',3,ifi,'(9f15.10)',gasa,ndim*ndim,ndim,idim,jdim)
C        call fclose(ifi)
C        ifi = fopna('out2',-1,0)
C        call ywrm(0,' ',3,ifi,'(9f15.10)',gf,ndim*ndim,ndim,idim,jdim)
C        call fclose(ifi)
C        call cwrite(' quit (q) ',0,9,0)
C        read(*,'(a20)',err=99,end=99) strn
C        if (strn == 'x') stop
C        if (strn == 'q') goto 99
C   50 continue
C
C   99 print *, 'done ib,jb=',ib,jb
C      end
