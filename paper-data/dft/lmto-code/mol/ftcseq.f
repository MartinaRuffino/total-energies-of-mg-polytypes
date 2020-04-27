      subroutine ftcseq(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,d0,ndist,adec,
     .  nb1,nb2,npgo,ri1,ri2,jfi)
C- Makes 2CF sequence of distances separately, writes on file jfi
      implicit real*8 (a-h,p-z), integer (o)
      dimension lx1(8),lx2(8),ex1(8),ex2(8),dist(120),err(100)
      real w(1)
      common /w/ w
      call getpr(ipr)
      if (ndist > 120) call rx('ftcseq: ndist too big')
      call ftcdis(d0,ndist,adec,dist)
CL      write(71,710) ndist,adec,d0,dist(ndist)
  710 format(' ndist',i4,'   adec',f8.3,'   d0',f9.4,'   dn',f9.4)

C ------ loop over distances ---
      do 50 idist=1,ndist
      d=dist(idist)
      if((ipr >= 30.and.(idist == 1.or.idist == ndist))
     .   .or.ipr >= 40) write(6,550) idist,ndist,d
  550 format(/' ftcseq:  distance',i4,' of',i4,'   d=',f9.6)
      lchk=0
      if(ipr >= 60) lchk=1
      jpr=ipr-20
      if(idist == 1) jpr=ipr
      call pshpr(jpr)

C ------ define the mesh ------
      npmx=2*nb1*nb2+600
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
C|    call bisinl(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
C|   .  zc1,zc2)
      call holint(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
     .  zc1,zc2)

c ------ make fit at this distance -----------
      lpr=0
      if(idist == 1) lpr=1
      call ftcgen(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,rmt1,rmt2,
     .  rsmp1,rsmp2,rsm1,rsm2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npgo,
     .  zc1,zc2,ri1,ri2,jfi,lchk,err,lpr)

      call rlse(oxp)
      call poppr
   50 continue
      end
