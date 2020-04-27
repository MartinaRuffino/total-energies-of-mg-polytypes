      subroutine gfqdys(n1,n2,n3,idim,npfun,pfun,pref,gij,gref,wij,gdy)
C- g = gref + g dP gref by FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1,n2,n3 QP mesh
Ci   idim     dimensions of g
Ci   pfun,npfun  potential functions, and leading dimension
Ci   pref     potential functions, reference G
Ci   gij      Green's function, on FFT mesh
Ci   gref     reference Green's function, on FFT mesh
Ci   wij      work array of same dimensions as gij,gref
Co  Outputs
Co    gdy     Dyson product g + g dP gref
Cr  Remarks
Cr    k-space:
Cr    mc -qr gref pref pfun -- -v2dia -x -qr g -x gref -+
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      integer idim,npfun,n1,n2,n3
      double complex pfun(npfun),pref(npfun)
      double complex gij(n1,n2,n3,idim,idim),wij(idim,n1,n2,n3,idim)
      double complex gref(n1,n2,n3,idim,idim),gdy(idim,idim)
C Local variables
      integer id,jd,i1,i2,i3,iset,lprod,jdim,i1m,i2m,i3m

      jdim = idim

C --- R.S. gf by inverse FFT, one orbital pair at a time ---
      iset = 0
      call fftz3s(n1,n2,n3,n1,n2,n3,iset)
      do  50  jd = 1, jdim
      do  50  id = 1, idim
        call fftz3(gij(1,1,1,id,jd),n1,n2,n3,n1,n2,n3,1,iset,-1)
        call fftz3(gref(1,1,1,id,jd),n1,n2,n3,n1,n2,n3,1,iset,-1)
   50 continue

C --- Copy gij dPj into format suitable for zgemm ---
      do  40  id = 1, idim
      do  40  jd = 1, jdim
      do  10  i3 = 1, n3
      i3m = n3+2-i3
      if (i3 == 1) i3m = 1

        do  20  i2 = 1, n2
        i2m = n2+2-i2
        if (i2 == 1) i2m = 1

          do  30  i1 = 1, n1
          i1m = n1+2-i1
          if (i1 == 1) i1m = 1

          wij(id,i1m,i2m,i3m,jd) = gij(i1,i2,i3,id,jd)*
     .                             (pref(jd)-pfun(jd))

   30     continue
   20   continue
   10 continue

   40 continue

C --- Calculate g_ik delP gref_kj ---
      lprod = n1*n2*n3*jdim
      call zgemm('N','N',idim,jdim,lprod,(1d0,0d0),wij,idim,gref,lprod,
     .  (0d0,0d0),gdy,idim)

      call yprm('delta g-dyson',3,gdy,idim*idim,idim,idim,idim)

C --- Add gref_ij to gdy ---
      do  60  jd = 1, jdim
      do  60  id = 1, idim
        gdy(id,jd) = gref(1,1,1,id,jd) + gdy(id,jd)
        if (abs(gdy(id,jd)-gij(1,1,1,id,jd)) > 1d-12) then
          print 333, id,jd,gdy(id,jd),gij(1,1,1,id,jd)
  333     format(2i4,' oops',4f14.8)
        endif
   60 continue

      call yprm('g-dyson',3,gdy,idim*idim,idim,idim,idim)

      end
