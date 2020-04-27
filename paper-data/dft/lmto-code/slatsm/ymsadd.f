      subroutine ymsadd(nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .  alpha,beta,src,offsi,dest,offdi)
C- Add a subblock block of a complex matrix to a destination
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma : row subblock size
Ci   nlmb : col subblock size
Ci   ndas : row dimension of source matrix
Ci   ndad : row dimension of destination matrix
Ci   ofas : offset to first row of source matrix
Ci   ofbs : offset to first col of source matrix
Ci   ofad : offset to first row of destination matrix
Ci   ofbd : offset to first col of destination matrix
Ci   alpha: complex scalar alpha scaling src; see Remarks
Ci   beta : complex scalar beta scaling dest; see Remarks
Ci   src  : source matrix
Ci   offsi: separation between real,imaginary parts of src
Ci   offdi: separation between real,imaginary parts of dest
Co Outputs
Co   dest : is overwritten by alpha * src + beta * dest
Cr Remarks
Cr   ymsadd adds a subblock alpha * src into a subblock of dest:
Cr   dest is overwritten by
Cr       dest <- alpha * src + beta * dest
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,offsi,offdi
      double precision src(ndas,1),dest(ndad,1),alpha(2),beta(2)
C ... Local parameters
      integer ia,ib,ofsi,ofdi

      if (nlma <= 0 .or. nlmb <= 0) return

      ofsi = ofas + offsi
      ofdi = ofad + offdi

      if (beta(1) == 0 .and. beta(2) == 0) then
C       dest = src
        if (alpha(1) == 1 .and. alpha(2) == 0) then
          do  10  ib = 1, nlmb
          do  10  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = src(ia+ofsi,ib+ofbs)
   10     continue
C       dest = -src
        elseif (alpha(1) == -1 .and. alpha(2) == 0) then
          do  11  ib = 1, nlmb
          do  11  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = -src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = -src(ia+ofsi,ib+ofbs)
   11     continue
C       dest = alpha * src
        else
          do  12  ib = 1, nlmb
          do  12  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = alpha(1)*src(ia+ofas,ib+ofbs) -
     .                              alpha(2)*src(ia+ofsi,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = alpha(1)*src(ia+ofsi,ib+ofbs) +
     .                              alpha(2)*src(ia+ofas,ib+ofbs)
   12     continue
        endif
      elseif (beta(1) == 1 .and. beta(2) == 0) then
C       dest = dest + src
        if (alpha(1) == 1 .and. alpha(2) == 0) then
          do  20  ib = 1, nlmb
          do  20  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = dest(ia+ofad,ib+ofbd) +
     .                              src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = dest(ia+ofdi,ib+ofbd) +
     .                              src(ia+ofsi,ib+ofbs)
   20     continue
C       dest = dest - src
        elseif (alpha(1) == -1 .and. alpha(2) == 0) then
          do  21  ib = 1, nlmb
          do  21  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = dest(ia+ofad,ib+ofbd) -
     .                              src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = dest(ia+ofdi,ib+ofbd) -
     .                              src(ia+ofsi,ib+ofbs)
   21     continue
C       dest = dest + alpha * src
        else
          do  22  ib = 1, nlmb
          do  22  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = dest(ia+ofad,ib+ofbd) +
     .                              alpha(1)*src(ia+ofas,ib+ofbs) -
     .                              alpha(2)*src(ia+ofsi,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = dest(ia+ofdi,ib+ofbd) +
     .                              alpha(1)*src(ia+ofsi,ib+ofbs) +
     .                              alpha(2)*src(ia+ofas,ib+ofbs)
   22     continue
        endif
C     dest = beta * dest + alpha * src
      else
        do  30  ib = 1, nlmb
        do  30  ia = 1, nlma
          dest(ia+ofad,ib+ofbd) = beta(1)*dest(ia+ofad,ib+ofbd) -
     .                            beta(2)*dest(ia+ofdi,ib+ofbd) +
     .                            alpha(1)*src(ia+ofas,ib+ofbs) -
     .                            alpha(2)*src(ia+ofsi,ib+ofbs)
          dest(ia+ofdi,ib+ofbd) = beta(1)*dest(ia+ofdi,ib+ofbd) +
     .                            beta(2)*dest(ia+ofad,ib+ofbd) +
     .                            alpha(1)*src(ia+ofsi,ib+ofbs) +
     .                            alpha(2)*src(ia+ofas,ib+ofbs)
   30   continue
      endif
      end
C#ifdefC TEST
C      subroutine fmain
CC      implicit none
C      integer lda,ldb,nr,nc,ir,ic,i1mach,ld1,ld2,ldr,ldr1,oi,oi1,kcplxi,
C     .  kcplxf,lerr,so,lda1,ldb1
C      parameter (lda=5,ldb=6,nr=4,nc=3,lda1=7,ldb1=8)
C      double precision alpha(2), beta(2)
C      double precision yy(lda,ldb,2),y(lda,2,ldb),z(2,lda,ldb)
C      double precision yy0(lda,ldb,2),diff
C      double precision yy1(lda1,ldb1,2),yy10(lda1,ldb1,2)
C      character*8 fm
C      equivalence (yy,y)
C      equivalence (yy,z)
C
C      so = i1mach(2)
C      fm = '(8f8.2)'
C
C      lerr = 0
C
C      do  10  ir = 1, lda
C      do  10  ic = 1, ldb
C        yy(ir,ic,1) = 100*ir + ic
C        yy(ir,ic,2) = -(100*ir + ic)
C
C        yy0(ir,ic,1) = 100*ir + ic
C        yy0(ir,ic,2) = -(100*ir + ic)
C
C   10 continue
C
C      do  20  ir = 1, lda1
C      do  20  ic = 1, ldb1
C        yy1(ir,ic,1) = 100*ir + ic
C        yy1(ir,ic,2) = -(100*ir + ic)
C
C        yy10(ir,ic,1) = 100*ir + ic
C        yy10(ir,ic,2) = -(100*ir + ic)
C
C   20 continue
C
C      do  30  kcplxi = 0, 2, 2
C      call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C      call cplxdm(kcplxi,lda1,ldb1,ld1,ld2,ldr1,oi1)
C      call ztoyy(yy,lda,ldb,lda,ldb,0,kcplxi)
C      call ztoyy(yy1,lda1,ldb1,lda1,ldb1,0,kcplxi)
C      call ztoyy(yy10,lda1,ldb1,lda1,ldb1,0,kcplxi)
C
CC ... straight copy, with an offset
C      alpha(1) = 1
C      alpha(2) = 0
C      beta(1) = 0
C      beta(2) = 0
C      call dvset(yy1,1,lda1*ldb1*2,-1d0)
C      call ymscop(4,2,3,ldr,ldr1,2,1,1,2,yy,oi,yy1,oi1)
C      call ywrm(0,'real copy yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
CC     call ymsadd(2,3,ldr,ldr1,2,1,1,2,alpha,beta,yy,oi,yy1,oi1)
C      call ymscop(0,2,3,ldr,ldr1,2,1,1,2,yy,oi,yy1,oi1)
C      call ywrm(0,'complex copy yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C
CC ... add with alpha=1-i, beta=1
C      alpha(1) = 1
C      alpha(2) = -1
C      beta(1) = 1
C      beta(2) = 0
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
C      call ywrm(0,'yy1+(1-i)yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C
CC ... with alpha=-1, beta=0
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      alpha(1) = -1
C      alpha(2) = 0
C      beta(1) = 0
C      beta(2) = 0
C      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
C      call ywrm(0,'- yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C
CC ... with alpha=-1, beta=1
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      alpha(1) = -1
C      alpha(2) = 0
C      beta(1) = 1
C      beta(2) = 0
C      call ymsadd(2,3,ldr,ldr1,2,1,1,2,alpha,beta,yy,oi,yy1,oi1)
C      call ywrm(0,'yy1 - yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
CC ... ditto, for real matrix using ymscop mode=5
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      call dscal(lda*ldb*2,-1d0,yy,1)
C      call ymscop(5,2,3,ldr,ldr1,2,1,1,2,yy,oi,yy1,oi1)
C      call ywrm(0,'real yy1 - yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C      call dscal(lda*ldb*2,-1d0,yy,1)
C
CC ... with beta=0
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      alpha(1) = 0
C      alpha(2) = 1
C      beta(1) = 0
C      beta(2) = 0
C      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
C      call ywrm(0,'i * yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C
CC ... add with alpha=beta=1
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      alpha(1) = 1
C      alpha(2) = 0
C      beta(1) = 1
C      beta(2) = 0
C      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
C      call ywrm(0,'y + yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C
CC ... add with general alpha,beta
C      alpha(1) = 1
C      alpha(2) = -1
C      beta(1) = -1d0/2
C      beta(2) = 1d0/2
C      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
C      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
C      call ywrm(0,'(1-i)(yy-yy1/2)',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C
C   30 continue
C
C      end
C#endif
