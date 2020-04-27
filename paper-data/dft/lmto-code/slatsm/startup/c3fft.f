      subroutine c2fft(c,id,nm,nn,wm,wn,isig,iord,iwork,ierr)
****purpose:
*       this routine performs a 2-dimensional complex fourier transform,
*       of order nm*nn.
****usage:
*       the user is expected to provide the data in a 2-dimensional
*       complex array c, dimensioned in the calling program c(id,nn);
*       id can be different from nm, and it is recommended that it is
*       chosen equal to nm+1 if nm is even or to nm if nm is odd.
*       the elements c(k,*),nm<k<=id must be zeroed. the routine is
*       intended for repeated usage, thus separate set-up and
*       operating calls are available : the user should always perform
*       a set-up call ( isig=0 ) passing the chosen parameters, before
*       performing the actual transform ( isig= +1 or -1 ); the user can
*       choose whether to obtain the results of the direct transform
*       in natural order (isig=-1,iord=1) or leave them in the
*      'bit-reversed' order( isig=-1,iord=0); this choice saves
*       some computer time, and it is recommended in cases discussed
*       in the long write-up. analogously,the inverse transform accepts
*       input ( please note! ) data in natural order ( isig=1,iord=1),
*       or data already subjected to a bit-reversal permutation( isig=1
*       iord=0 ).
****parameters :
*       input :
*       c : array to be transformed; declared complex c(id,nn) in the
*           calling program;
*       id : first dimension of c in the calling program
*       isig : option flag : isig=0 : set-up run, c not used
*                            isig=-1: direct transform
*                            isig=+1: inverse transform
*       wm,wn : integer arrays , used to host tables for the transform;
*               dimensioned in the calling program at least (4*nm+14)
*               and (4*nn+14) respectively; if isig /= 0 ,they are
*               assumed to have been set by a previous call with isig=0
*               (all other arguments unchanged), and never have been
*               modified since then;
*               if nm=nn, wm and wn do not need to be distinct;
*       nm : order of the transform along the columns of c
*       nn : order of the transform along the rows of c
*       iord : option flag : =1 : output in natural order (isig=-1)
*                                 input in natural order  (isig=+1)
*                            =0 : output in bit-reversed order(isig=-1)
*                                 input in bit-reversed order(isig=+1)
*       iwork : integer array, used as work area for reordering if iord=
*                1, it must be at least max(nm,nn) words long.
*                it is unused if iord=0.
*
****output :
*       c : transformed array
*       wm, wn : only if isig=0, wm and wn are filled with the
*               appropriate  tables (twiddle factors).
*       iwork : undefined
*       ierr  : error code : =0 : successful
*                            =1 : data dimensions are not correct
*                            =2 : prime factors different from 2,3,5
*                                 are present in data dimensions
*                            =3 : tables not correctly initialized
*
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision wm(-14:*),wn(-14:*)
      integer iwork(*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      if(id < nm)then
         ierr=iderr
         return
      endif
      ierr=0
      if(isig == 0)then
         call mfftp(nm,wm(0),wm,0,ierr)
         if(ierr /= 0)return
         if(nn /= nm)then
         call mfftp(nn,wn(0),wn,0,ierr)
         if(ierr /= 0)return
         else
         call mfftz0(wm,1,2*nm+7,wn,1)
         endif
         return

      else if(isig > 0)then

         if(iord /= 0)then
         call mfftov(c,1,id,nm,nn,wm(3*nm),iwork)
         call mfftov(c,id,1,nn,nm,wn(3*nn),iwork)
         endif

         call mfftiv(c,1,id,nm,nn,wm,wm(0),ierr)
         if(ierr /= 0)return
         call mfftiv(c,id,1,nn,nm,wn,wn(0),ierr)
         if(ierr /= 0)return

      else
         call mfftdv(c,1,id,nm,nn,wm,wm(0),ierr)
         if(ierr /= 0)return
         call mfftdv(c,id,1,nn,nm,wn,wn(0),ierr)
         if(ierr /= 0)return

         if(iord /= 0)then
         call mfftov(c,1,id,nm,nn,wm(2*nm),iwork)
         call mfftov(c,id,1,nn,nm,wn(2*nn),iwork)
         endif
      endif

      end

      subroutine c3fft(c,id,nl,nm,nn,wl,wm,wn,iopt,isig,iord,iwork,ierr)

****purpose:
*       this routine performs a 3-dimensional complex fourier transform,
*       of order nl*nm*nn  .
*       it is a user interface for a package (mfft), specially designed
*       for high performance on cray x-mp machines.
****usage:
*       the user is expected to provide the data in a 3-dimensional
*       complex array c, dimensioned in the calling program c(id,nm,nn);
*       id can be different from nl, and it is recommended that it is
*       chosen odd for maximum performance. the routine is
*       intended for repeated usage, thus separate set-up and
*       operating calls are available : the user should always perform
*       a set-up call ( isig=0 ) passing the chosen parameters, before
*       performing the actual transform ( isig= +1 or -1 ); the user can
*       choose whether to obtain the results of the direct transform
*       in natural order (isig=-1,iord=1) or leave them in the
*       bit-reversed  order( isig=-1,iord=0); this choice saves
*       some computer time, and it is recommended in cases discussed
*       in the long write-up. analogously, the inverse transform accepts
*       input ( please note! ) data in natural order ( isig=1,iord=1),
*       or data already subjected to a bit-reversal permutation( isig=1
*       iord=0).
*       a special treatment is available to speed up the transform of
*       small matrices. this treatment is activated by the flag iopt. in
*       this case the tables for the second dimension ( wm ) are larger,
*       but the increase in performance is substantial when nm<32.
****arguments :
*       input :
*       c : array to be transformed; declared complex c(id,nm,nn) in the
*           calling program;
*       id : first dimension of c in the calling program
*       isig : option flag : isig=0 : set-up run, c not used
*                            isig=-1: direct transform
*                            isig=+1: inverse transform
*       wl,wm,wn : integer arrays,used to host tables for the transforms
*               dimensioned in the calling program at least (4*nl+14)
*               (4*nm+14) and (4*nn+14) respectively; if iopt=1
*               wm must be dimensioned at least (4*nm*(id+1)+14)
*               if isig /= 0, they are assumed to have been set by a
*               previous call with isig=0 and other arguments equal, and
*               never have been modified ;
*               when the corresponding orders are equal, they do not
*               need to be distinct.
*       nl : order of the transform along the columns of c
*       nm : order of the transform along the rows of c
*       nn : order of the transform along the third dimension of c
*       iopt : option flag : =0 : normal treatment
*                            =1 : special treatment for improving
*                                 vectorization on matrices with
*                                 small nl; requires long wm(see);if
*                                 requested, must be present in both
*                                 the set-up and transform calls;
*       iord : option flag : =1 : output in natural order (isig=-1)
*                                 input in natural order  (isig=+1)
*                            =0 : output in bit-reversed order(isig=-1)
*                                 input in bit-reversed order(isig=+1)
*       iwork : integer array, used as work area for reordering if
*               iord=1; it must be at least max(nl,nm,nn) words long.
*
*        output :
*       c : transformed array
*       wl, wm, wn : only if isig=0, wl,wm and wn are filled with the
*                     appropriate tables
*       iwork : undefined
*       ierr  : error code : =0 : successful
*                          : =1 : data dimensions are not correct
*                            =2 : prime factors different from 2,3,5
*                                 are present in data dimensions
*                            =3 : tables not correctly initialized

      implicit double precision (a-h,o-z)
        double complex c(0:*)
        double precision wm(-14:*),wn(-14:*),wl(-14:*)
        integer iwork(*)
        integer iderr,facerr,tberr
        parameter (iderr=1,facerr=2,tberr=3)

        if(id < nl)then
        ierr=iderr
        return
        endif
        ierr=0
        if(isig == 0)then
        call mfftp(nm,wm(0),wm,id*iopt,ierr)
        if(ierr /= 0)return
        if(nm /= nn)then
        call mfftp(nn,wn(0),wn,0,ierr)
        if(ierr /= 0)return
        else
        call mfftz0(wm,1,2*nm+7,wn,1)
        endif
        if(nl == nm)then
        call mfftz0(wm,1,2*nm+7,wl,1)
        else if(nl == nn)then
        call mfftz0(wn,1,2*nn+7,wl,1)
        else
        call mfftp(nl,wl(0),wl,0,ierr)
        if(ierr /= 0)return
        endif
        return

        else   if(isig > 0)then

        if(iord /= 0)then
        call mfftov(c,1,id,nl,nm*nn,wl(3*nl),iwork)
        call mfftov(c,id*nm,1,nn,id*nm,wn(3*nn),iwork)
        call mfftom(c,id,id*nm,1,nm,nn,nl,wm(3*nm),
     $                iwork)
        endif
        call mfftiv(c,1,id,nl,nm*nn,wl,wl(0),ierr)
        if(ierr /= 0)return
        if(iopt == 0)then
        call mfftim(c,id,id*nm,1,nm,nn,nl,wm,wm(0),ierr)
        if(ierr /= 0)return
        else
        call mfftis(c,id,id*nm,1,nm,nn,nl,wm,wm(0),ierr)
        if(ierr /= 0)return
        endif
        call mfftiv(c,id*nm,1,nn,id*nm,wn,wn(0),ierr)
        if(ierr /= 0)return

      else

        call mfftdv(c,1,id,nl,nm*nn,wl,wl(0),ierr)
        if(ierr /= 0)return
        if(iopt == 0)then
        call mfftdm(c,id,id*nm,1,nm,nn,nl,wm,wm(0),ierr)
        if(ierr /= 0)return
        else
        call mfftds(c,id,id*nm,1,nm,nn,nl,wm,wm(0),ierr)
        if(ierr /= 0)return
        endif
        call mfftdv(c,id*nm,1,nn,id*nm,wn,wn(0),ierr)
        if(ierr /= 0)return
        if(iord /= 0)then
        call mfftov(c,1,id,nl,nm*nn,wl(2*nl),iwork)
        call mfftov(c,id*nm,1,nn,id*nm,wn(2*nn),iwork)
        call mfftom(c,id,id*nm,1,nm,nn,nl,wm(2*nm),
     $                iwork)
        endif

      endif

      end

      subroutine r2fft(c,id,nm,nn,wm,wn,isig,iord,iwork,ierr)
*
****purpose:
*       this routine performs a 2-dimensional real fourier transform,
*       of order nm*nn.
****usage:
*       the user is expected to provide the data in a 2-dimensional
*       real array c, dimensioned in the calling program c(id,nn)
*       and equivalenced to a complex array ccom(id/2,nn).
*       when a direct trasnform is performed, half of the output
*       coefficients are contained in ccom(i,j) as follows
*        x(i,j)= ccom(i,j) with i=1,(nm/2+1) and j=1,nn
*       the remaining ones are obtained exploiting the conjugated
*       symmetry relation
*        x(i,j) = conjg(ccom(nm-i+2,1+mod(nn-j+1,nn)))
*       with i=nm/2+2,nm  and  j=1,nn
*       nm (actual first dimension of data) must be an even number;
*       id (declared first dimension) must be gr. or eq. to nm+2.
*       the elements c(k,*), nm<k<=id must be zeroed. the routine is
*       intended for repeated usage, thus separate set-up and
*       operating calls are available : the user should always perform
*       a set-up call ( isig=0 ) passing the chosen parameters, before
*       performing the actual transform ( isig= +1 or -1 ); the user can
*       choose whether to obtain the results of the direct transform
*       in natural order (isig=-1,iord=1) or leave them in the
*       bit-reversed  order( isig=-1,iord=0); this choice saves
*       some computer time, and it is recommended in cases discussed
*       in the long write-up. analogously, the inverse transform accepts
*       input ( please note! ) data in natural order ( isig=1,iord=1),
*       or data already subjected to a bit-reversal permutation( isig=1
*       iord=0 ).
****parameters :
*       input :
*       c : array to be transformed; declared real c(id,nn) in the
*           calling program;
*       id : first dimension of c in the calling program
*            it must be an even number equal to nm+2.
*       isig : option flag : isig=0 : set-up run, c not used
*                            isig=-1: direct transform
*                            isig=+1: inverse transform
*       wm,wn : integer arrays , used to host tables for the transform;
*               dimensioned in the calling program in this way :
*               wm : at least (6*nm+14)
*               wn : at least (4*nn+14)
*               if isig /= 0 ,they are assumed to have been set
*                by a previous call with isig=0 and
*               all the other arguments unchanged, and never have been
*               modified since then;
*               if nm=nn, wm and wn do not need to be distinct;
*       nm : order of the transform along the columns of c
*       nn : order of the transform along the rows of c
*       iord : option flag : =1 : output in natural order (isig=-1)
*                                 input in natural order  (isig=+1)
*                            =0 : output in bit-reversed order(isig=-1)
*                                 input in bit-reversed order(isig=+1)
*                                 warning:the first dimension is ordered
*                                 in any case because of post-processing
*                                 (for direct) and pre-processing for
*                                 inverse real transforms
*       iwork : integer array, used as work area for reordering if iord=
*               1, must be at least max(nm,nn) words long.
*
****output :
*       c : transformed array
*       wm, wn : only if isig=0, wm and wn are filled with the
*               appropriate  tables
*       iwork : undefined
*       ierr  : error code : =0 : successful
*                            =1 : wrong id parameter
*                            =2 : prime factors different from 2,3,5
*                                 are present in data dimensions
*                            =3 : tables not correctly initialized
*                            =4 : first dimension is an odd number
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision wm(-14:*),wn(-14:*)
      integer iwork(*)
*
      integer iderr,facerr,tberr,odderr
      parameter (iderr=1,facerr=2,tberr=3,odderr=4)
*
      if(id < nm+2)then
        ierr=iderr
        return
      endif
      nm1=nm/2
      if(nm1*2 /= nm)then
        ierr=odderr
        return
      endif
      ierr=0
*
      if(isig == 0) then
        call mfftrp(nm,wm(4*nm))
        call mfftp(nm1,wm(0),wm,0,ierr)
        if(ierr /= 0)return
        if(nn /= nm1) then
          call mfftp(nn,wn(0),wn,0,ierr)
          if(ierr /= 0)return
        else
          call mfftz0(wm,1,2*nm1+7,wn,1)
        endif
        return
      else  if(isig > 0) then

        if(iord /= 0) then
          call mfftov(c,id/2,1,nn,nm1+1,wn(3*nn),iwork)
        endif
        call mfftiv(c,id/2,1,nn,nm1+1,wn,wn(0),ierr)
        if(ierr /= 0)return
        call mfftri(c,1,id/2,nm1,nn,wm(4*nm))
        call mfftov(c,1,id/2,nm1,nn,wm(3*nm1),iwork)
        call mfftiv(c,1,id/2,nm1,nn,wm,wm(0),ierr)
        if(ierr /= 0)return

      else
        call mfftdv(c,1,id/2,nm1,nn,wm,wm(0),ierr)
        if(ierr /= 0)return
        call mfftov(c,1,id/2,nm1,nn,wm(nm1*2),iwork)
        call mfftrd(c,1,id/2,nm1,nn,wm(4*nm))
        call mfftdv(c,id/2,1,nn,nm1+1,wn,wn(0),ierr)
        if(ierr /= 0)return

        if(iord /= 0) then
          call mfftov(c,id/2,1,nn,nm1+1,wn(nn*2),iwork)
        endif
*
      endif
*
      end
c     ###########    fft1 ends here     ###################

      subroutine r3fft(c,id,nl,nm,nn,wl,wm,wn,iopt,isig,iord,iwork,ierr)
*
****purpose:
*       this routine performs a 3-dimensional real fourier transform,
*       of order nl*nm*nn  .
****usage:
*       the user is expected to provide the data in a 3-dimensional
*       real array c, dimensioned in the calling program c(id,nm,nn);
*       id has to be an even integer, equal to nl+2.
*       for output data arrengement see notes to r2fft here above.
*       this routine is
*       intended for repeated usage, thus separate set-up and
*       and operating calls are available: the user should in any case
*       perform a set-up call (isig=0) passing the parameters before
*       performing an actual transform ( isig= +1 or -1 ); the user can
*       choose whether to obtain the results of the direct transform
*       in natural order (isig=-1,iord=1) or leave them in the
*       bit-reversed  order( isig=-1,iord=0); this choice saves
*       some computer time, and it is recommended in cases discussed
*       in the long write-up. analogously, the inverse transform accepts
*       input ( please note! ) data in natural order ( isig=1,iord=1),
*       or data already subjected to a bit-reversal permutation( isig=1
*       iord=0).
*       a special treatment is available to speed up the transform of
*       small matrices. this treatment is activated by the flag iopt. in
*       this case the tables for the second dimension ( wm ) are larger,
*       but the increase in performance is substantial when nm<32.
****arguments :
*       input :
*       c : array to be transformed; declared complex c(id,nm,nn) in the
*           calling program;
*       id : first dimension of c in the calling program
*            it has to be an even integer >= nl+2.
*       isig : option flag : isig=0 : set-up run, c not used
*                            isig=-1: direct transform
*                            isig=+1: inverse transform
*       wl,wm,wn : integer arrays,used to host tables for the transforms
*               dimensioned in the calling program at least (6*nl+14)
*               (4*nm+14) and (4*nn+14) respectively; if
*               iopt=1, wm must be dimensioned at least 4*nm*(id/2+1)+14
*               if isig /= 0, they are assumed to have been set by a
*               previous call with isig=0 and other arguments equal, and
*               never have been modified ;
*               when the corresponding orders are equal, they do not
*               need to be distinct
*       nl : order of the transform along the columns of c
*            it has to be an even integer.
*       nm : order of the transform along the rows of c
*       nn : order of the transform along the third dimension of c
*       iopt : option flag : =0 : normal treatment
*                            =1 : special treatment for improving
*                                 vectorization on matrices with
*                                 small nl; requires long wm(see);if
*                                 requested, must be present in both
*                                 the set-up and transform calls;
*       iord : option flag : =1 : output in natural order (isig=-1)
*                                 input in natural order  (isig=+1)
*                            =0 : output in bit-reversed order(isig=-1)
*                                 input in bit-reversed order(isig=+1)
*       iwork : integer array, used as work area for reordering if
*               iord=1; it must be at least max(nl,nm,nn) words long.
*
*        output :
*       c : transformed array
*       wl, wm, wn : only if isig=0, wl,wm and wn are filled with the
*                     appropriate tables
*       iwork : undefined
*       ierr  : error code : =0 : successful
*                            =1 : wrong id parameter
*                            =2 : prime factors different from 2,3,5
*                                 are present in data dimensions
*                            =3 : tables not correctly initialized
*                            =4 : first dimension is an odd number
*
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision wl(-14:*),wm(-14:*),wn(-14:*)
      integer  iwork(*)
*
*
      integer iderr,facerr,tberr,odderr
      parameter (iderr=1,facerr=2,tberr=3,odderr=4)
*
      if(id < nl+2)then
        ierr=iderr
        return
      endif
      nl1=nl/2
      if(nl1*2 /= nl)then
        ierr=odderr
        return
      endif
      ierr=0

      nmpn=nm*nn
*
*
      if(isig == 0) then
        call mfftp(nm,wm(0),wm,id/2*iopt,ierr)
        if(ierr /= 0)return

        call mfftrp(nl,wl(4*nl))
        if(nl1 /= nm) then
          call mfftp(nl1,wl(0),wl,0,ierr)
          if(ierr /= 0)return
        else
          call mfftz0(wm,1,2*nm+7,wl,1)
        endif
*
        if(nm == nn) then
          call mfftz0(wm,1,2*nm+7,wn,1)
        else if(nn == nl1) then
          call mfftz0(wl,1,2*nl1+7,wn,1)
        else
          call mfftp(nn,wn(0),wn,0,ierr)
          if(ierr /= 0)return
        endif
        return
*
      else   if(isig > 0) then
*
        if(iord /= 0) then
          call mfftom(c,id/2,id/2*nm,1,nm,nn,nl1+1,wm(nm*3),iwork)
          call mfftov(c,id/2*nm,1,nn,id/2*nm,wn(nn*3),iwork)
        endif
*
        call mfftiv(c,id/2*nm,1,nn,id/2*nm,wn,wn(0),ierr)
        if(ierr /= 0)return
*
        if(iopt == 0) then
          call mfftim(c,id/2,id/2*nm,1,nm,nn,nl1+1,wm,wm(0),ierr)
          if(ierr /= 0)return
        else
          call mfftis(c,id/2,id/2*nm,1,nm,nn,nl1+1,wm,wm(0),ierr)
          if(ierr /= 0)return
        endif
*
        call mfftri(c,1,id/2,nl1,nmpn,wl(4*nl))
        call mfftov(c,1,id/2,nl1,nmpn,wl(nl1*3),iwork)
        call mfftiv(c,1,id/2,nl1,nmpn,wl,wl(0),ierr)
        if(ierr /= 0)return
*
*
      else
*
*
        call mfftdv(c,1,id/2,nl1,nmpn,wl,wl(0),ierr)
        if(ierr /= 0)return
        call mfftov(c,1,id/2,nl1,nmpn,wl(nl1*2),iwork)
        call mfftrd(c,1,id/2,nl1,nmpn,wl(4*nl))
*
*
        if(iopt == 0) then
          call mfftdm(c,id/2,id/2*nm,1,nm,nn,nl1+1,wm,wm(0),ierr)
          if(ierr /= 0)return
        else
          call mfftds(c,id/2,id/2*nm,1,nm,nn,nl1+1,wm,wm(0),ierr)
          if(ierr /= 0)return
        endif
*
        call mfftdv(c,id/2*nm,1,nn,id/2*nm,wn,wn(0),ierr)
        if(ierr /= 0)return
*
        if(iord /= 0) then
          call mfftov(c,id/2*nm,1,nn,id/2*nm,wn(nn*2),iwork)
          call mfftom(c,id/2,id/2*nm,1,nm,nn,nl1+1,wm(nm*2),iwork)
        endif
*
      endif
*
      end
      subroutine mfftdm(c,imsx,ivsx,iesx,nmx,nvx,nex,tables,w,ierr)
*   purpose:
*       this subroutine performs a direct fourier transform along
*       the second dimension of a 3-dimensional matrix, using the
*       gentleman-sande algorithm.
*       the sequence to be transformed is c[imsx,nmx], whose components
*       are the 2-vectors c(m)[ivsx,nvx [iesx,nex]].
*       see ref.[1] for notations.
*  example:
*       let c be a 3-d matrix c(n1,n2,n3), declared via
*              dimension c(nn1,n2,n3)
*       with nn1 >= n1
*       then its dft along the second dimension is obtained by
*              call mfftdm(c,nn1,nn1*n2,1,n2,n3,n1,tables,ierr)
*  implementation:
*       the transformation is implemented through repeated calls to the
*       "butterfly" transformation mfft?6; parameters of the "butterfly"
*       are communicated through the common block mfftpa.
*  arguments:
*    input  :
*       c : array to be transformed.
*       imsx,ivsx,iesx,nmx,nvx,nex: these arguments define the structure
*           of c according to the definitions above. they are unchanged
*           on output
*       tables : array prepared by mfftp. it  is not changed on output.
*                it should be declared integer tables(4*nm+14);
*                it must be initialized by mfftp before usage.
*    output:
*       c : transform of the original array; "bit reversed" order
*       ierr : error code : =0 : successful
*                         : =3 :  'tables' not correctly initialized
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision w(0:*)
      integer tables(-14:*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
*  loading the common block : constants
      ims=imsx
      ivs=ivsx
      ies=iesx
      nm=nmx
      nv=nvx
      ne=nex

      ivlim=(nv-1)*ivs
      ilim=(ne-1)*ies
      mstep=ims
*  loading the common block : iteration-dependent quantities: initializa
      lx=1
      mx=nm
*  select the highest factor of nm
      ifac=tables(-1)
      goto(200,300,500)ifac
      ierr=tberr
      return
*...  radix 5 loop
500   continue
      do 510 im=1,tables(-12)
        mx=mx/5
        nustep=mx*mstep
        mlim=nustep-mstep
        lstep=nustep*5
        llim=nm*mstep-lstep
        call mfftc6(c,w)
        lx=lx*5
510   continue
*..  radix 3 loop
300   continue
      do 310 im=1,tables(-13)
        mx=mx/3
        nustep=mx*mstep
        mlim=nustep-mstep
        lstep=nustep*3
        llim=nm*mstep-lstep
        call mfftb6(c,w)
        lx=lx*3
310   continue
*..  radix 2 loop
200   continue
      do 210 im=1,tables(-14)
        mx=mx/2
        nustep=mx*mstep
        mlim=nustep-mstep
        md2lim=nustep/2-mstep
        lstep=nustep*2
        llim=nm*mstep-lstep
        call mffta6(c,w)
        lx=lx+lx
210   continue
      end
      subroutine mfftds(c,imsx,ivsx,iesx,nmx,nvx,nex,tables,w,ierr)

*   purpose:
*       the same as mfftdm.    it is a variant of mfftdm
*       aimed at optimal performance on matrices having the
*       first dimension smaller than 64. it requires that mfftp has
*       been called with id /= 0
*  implementation:
*       the transformation is implemented through repeated calls to the
*       "butterfly" transformation mfft?8; parameters of the "butterfly"
*       are communicated through the common block mfftpa.
*  arguments:
*    input  :
*       c : array to be transformed.
*       imsx,ivsx,iesx,nmx,nvx,nex: these arguments define the structure
*           of c according to the definitions above. they are unchanged
*           on output
*       tables : array prepared by mfftp. it  is not changed on output.
*                it should be declared integer tables(4*nm+14);
*                it must be initialized by mfftp before usage.
*    output:
*       c : transform of the original array; "bit reversed" order
*       ierr : error code : =0 : successful
*                         : =3 :  'tables' not correctly initialized
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision w(0:*)
      integer tables(-14:*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      save ibase
*  loading the common block : constants
      ims=imsx
      ivs=ivsx
      nm=nmx
      nv=nvx
      ne=nex
      ivlim=(nv-1)*ivs
      mstep=ims
*  loading the common block : iteration-dependent quantities: initializa
      lx=1
      mx=nm

      ibase=2*nm
*  select the highest factor of nm
      ifac=tables(-1)
      goto(200,300,500)ifac
      ierr=tberr
      return
*...  radix 5 loop
500   continue
      do 510 im=1,tables(-12)
        mx=mx/5
        nustep=mx*mstep
        ilim=nustep-1
        lstep=nustep*5
        llim=nm*mstep-lstep
        call mfftc8(c,w(2*ibase))
        lx=lx*5
        ibase=ibase+nustep*4
510   continue
*..  radix 3 loop
300   continue
      do 310 im=1,tables(-13)
        mx=mx/3
        nustep=mx*mstep
        ilim=nustep-1
        lstep=nustep*3
        llim=nm*mstep-lstep
        call mfftb8(c,w(2*ibase))
        lx=lx*3
        ibase=ibase+nustep*2
310   continue
*..  radix 2 loop
200   continue
      do 210 im=1,tables(-14)
        mx=mx/2
        nustep=mx*mstep
        ilim=nustep-1
        md2lim=nustep/2-mstep
        lstep=nustep*2
        llim=nm*mstep-lstep
        call mffta8(c,w(2*ibase))
        lx=lx+lx
        ibase=ibase+nustep
210   continue

      end
      subroutine mfftdv(c,ivsx,iesx,nvx,nex,tables,w,ierr)
*   purpose:
*       this subroutine performs a direct fourier transform along
*       one dimension of a 2-dimensional matrix, using the
*       gentleman-sande algorithm; no reordering is performed.
*       the sequence to be transformed is c[ivsx,nvx], whose components
*       are the vectors c(m)[iesx,nex].
*       see ref.[1] for notations.
*  example:
*       let c be a 2-d matrix c(n1,n2) declared via
*                 dimension c(id,n2)
*       with id >= n1.
*       then the dft along the first dimension is obtained by
*                 call mfftdv(c,1,id,n1,n2,tables,ierr)
*       the dft along the second dimension is obtained by
*                 call mfftdv(c,id,1,n2,n1,tables,ierr)
*  implementation:
*       the transformation is implemented through repeated calls to the
*       "butterfly" transformation mfft?4; parameters of the "butterfly"
*       are communicated through the common block mfftpa.
*  arguments:
*       c : array to be transformed.
*       ivsx,iesx,nvx,nex: these arguments define the structure of
*           c according to the definitions above. they are unchanged on
*           output;
*       tables : array prepared by mfftp. it  is not changed on output.
*                it should be declared integer tables(4*nm+14);
*                it must be initialized by mfftp before usage.
*    output:
*       c : transform of the original array; "bit reversed" order
*       ierr : error code : =0 : successful
*                         : =3 :  'tables' not correctly initialized
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision w(0:*)
      integer tables(-14:*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
*  loading the common block : constants
      ivs=ivsx
      ies=iesx
      nv=nvx
      ne=nex
      ilim=(ne-1)*ies
      mstep=ivs
*  loading the common block : iteration-dependent quantities: initializa
      lx=1
      mx=nv
*  select the highest factor of nv
      ifac=tables(-1)
      goto(200,300,500)ifac
      ierr=tberr
      return
*...  radix 5 loop
500   continue
      do 510 im=1,tables(-12)
        mx=mx/5
        nustep=mx*mstep
        mlim=nustep-mstep
        lstep=nustep*5
        llim=nv*mstep-lstep
        call mfftc4(c,w)
        lx=lx*5
510   continue
*..  radix 3 loop
300   continue
      do 310 im=1,tables(-13)
        mx=mx/3
        nustep=mx*mstep
        mlim=nustep-mstep
        lstep=nustep*3
        llim=nv*mstep-lstep
        call mfftb4(c,w)
        lx=lx*3
310   continue
*..  radix 2 loop
200   continue
      do 210 im=1,tables(-14)
        mx=mx/2
        nustep=mx*mstep
        mlim=nustep-mstep
        md2lim=nustep/2-mstep
        lstep=nustep*2
        llim=nv*mstep-lstep
        call mffta4(c,w)
        lx=lx+lx
210   continue
      end

      subroutine mfftim(c,imsx,ivsx,iesx,nmx,nvx,nex,tables,w,ierr)

*   purpose:
*       this subroutine performs an inverse fourier transform along
*       the second dimension of a 3-dimensional matrix, using the
*       cooley-tukey algorithm.
*       the input sequence is assumed to have been subjected to a
*       "bit reversal" permutation, through a call to mfftom, or
*       because it is the output of mfftdm.
*       the sequence to be transformed is c[imsx,nmx], whose components
*       are the 2-vectors c(m)[ivsx,nvx [iesx,nex]].
*       see ref.[1] for notations.
*  example:
*       let c be a 3-d matrix c(n1,n2,n3), declared via
*              dimension c(nn1,n2,n3)
*       with nn1 >= n1
*       then its idft along the second dimension is obtained by
*              call mfftim(c,nn1,nn1*n2,1,n2,n3,n1,tables,ierr)
*  implementation:
*       the transformation is implemented through repeated calls to the
*       "butterfly" transformation mfft?7; parameters of the "butterfly"
*       are communicated through the common block mfftpa.
*  arguments:
*    input  :
*       c : array to be transformed.
*       imsx,ivsx,iesx,nmx,nvx,nex: these arguments define the structure
*           of c according to the definitions above. they are unchanged
*           on output
*       tables : array prepared by mfftp. it  is not changed on output.
*                it should be declared integer tables(4*nm+14);
*                it must be initialized by mfftp before usage.
*    output:
*       c : transform of the original array;
*       ierr : error code : =0 : successful
*                         : =3 :  'tables' not correctly initialized
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision w(0:*)
      integer tables(-14:*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
*  loading the common block : constants
      ims=imsx
      ivs=ivsx
      ies=iesx
      nm=nmx
      nv=nvx
      ne=nex
      ivlim=(nv-1)*ivs
      ilim=(ne-1)*ies
      lstep=ims
*  loading the common block : iteration-dependent quantities: initializa
      lx=1
      mx=nm
      ifac=tables(-1)
      if(ifac > 3)then
        ierr=tberr
        return
      endif
*..   radix 2 loop
200   continue
      do 210 im=1,tables(-14)
        mx=mx/2
        nustep=lx*lstep
        llim=nustep-lstep
        ld2lim=nustep/2-lstep
        mstep=nustep*2
        mlim=nm*lstep-mstep
        call mffta7(c,w)
        lx=lx+lx
210   continue
      if(ifac == 1)return
*..   radix 3 loop
300   continue
      do 310 im=1,tables(-13)
        mx=mx/3
        nustep=lx*lstep
        llim=nustep-lstep
        mstep=nustep*3
        mlim=nm*lstep-mstep
        call mfftb7(c,w)
        lx=lx*3
310   continue
      if(ifac == 2)return
*..   radix 5 loop
500   continue
      do 510 im=1,tables(-12)
        mx=mx/5
        nustep=lx*lstep
        llim=nustep-lstep
        mstep=nustep*5
        mlim=nm*lstep-mstep
        call mfftc7(c,w)
        lx=lx*5
510   continue
      end

      subroutine mfftis(c,imsx,ivsx,iesx,nmx,nvx,nex,tables,w,ierr)
*   purpose:
*       the same as mfftim.  it is a variant of mfftim, optimized for
*       maximum performance on matrices whose first dimension is
*       less than 64.
*       it requires that tables has been prepared by a call to mfftp
*       with id /= 0 .
*  implementation:
*       the transformation is implemented through repeated calls to the
*       "butterfly" transformation mfft?9; parameters of the "butterfly"
*       are communicated through the common block mfftpa.
*  arguments:
*    input:
*       c : array to be transformed.
*       imsx,ivsx,iesx,nmx,nvx,nex: these arguments define the structure
*           of c according to the definitions above. they are unchanged
*           on output
*       tables : array prepared by mfftp. it  is not changed on output.
*                it should be declared integer tables(4*nm+14);
*                it must be initialized by mfftp before usage.
*    output:
*       c : transform of the original array; "bit reversed" order
*       ierr : error code : =0 : successful
*                         : =3 :  'tables' not correctly initialized
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision w(0:*)
      integer tables(-14:*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      save ibase
*  loading the common block : constants
      ims=imsx
      ivs=ivsx
      ies=1
      nm=nmx
      nv=nvx
      ne=nex

      ivlim=(nv-1)*ivs
      lstep=ims
*  loading the common block : iteration-dependent quantities: initializa
      lx=1
      mx=nm

      ibase=(2+ims)*nm
      ifac=tables(-1)
      if(ifac > 3)then
        ierr=tberr
        return
      endif
*..   radix 2 loop
200   continue
      do 210 im=1,tables(-14)
        mx=mx/2
        nustep=lx*lstep
        ilim=nustep-1
        mstep=nustep*2
        mlim=nm*lstep-mstep
        call mffta9(c,w(2*ibase))
        lx=lx+lx
        ibase=ibase+nustep
210   continue
      if(ifac == 1)return
*..   radix 3 loop
300   continue
      do 310 im=1,tables(-13)
        mx=mx/3
        nustep=lx*lstep
        ilim=nustep-1
        mstep=nustep*3
        mlim=nm*lstep-mstep
        call mfftb9(c,w(2*ibase))
        lx=lx*3
        ibase=ibase+nustep*2
310   continue
      if(ifac == 2)return
*..   radix 5 loop
500   continue
      do 510 im=1,tables(-12)
        mx=mx/5
        nustep=lx*lstep
        ilim=nustep-1
        mstep=nustep*5
        mlim=nm*lstep-mstep
        call mfftc9(c,w(2*ibase))
        lx=lx*5
        ibase=ibase+nustep*4
510   continue
      end
      subroutine mfftiv(c,ivsx,iesx,nvx,nex,tables,w,ierr)
*
*   purpose:
*       this subroutine performs an inverse fourier transform along
*       one dimension of a 2-dimensional matrix, using the
*       cooley-tukey algorithm.
*       the input matrix is assumed to have been subjected
*        to a "bit-reversal" reordering ( through mfftov or
*       because it is the output of mffdv and has not been reordered)
*       the sequence to be transformed is c[ivsx,nvx], whose components
*       are the vectors c(m)[iesx,nex].
*       see ref.[1] for notations.
*  example:
*       let c be a 2-d matrix c(n1,n2) declared via
*                 dimension c(id,n2)
*       with id >= n1.
*       then the idft along the first dimension is obtained by
*                 call mfftiv(c,1,id,n1,n2,tables,ierr)
*       the idft along the second dimension is obtained by
*                 call mfftiv(c,id,1,n2,n1,tables,ierr)
*  implementation:
*       the transformation is implemented through repeated calls to the
*       "butterfly" transformation mfft?5; parameters of the "butterfly"
*       are communicated through the common block mfftpa.
*  arguments:
*       c : array to be transformed.
*       ivsx,iesx,nvx,nex: these arguments define the structure of
*           c according to the definitions above. they are unchanged on
*           output;
*       tables : array prepared by mfftp. it  is not changed on output.
*                it should be declared integer tables(4*nm+14);
*                it must be initialized by mfftp before usage.
*  output:
*       c : transform of the original array; "bit reversed" order
*       ierr : error code : =0 : successful
*                         : =3 :  'tables' not correctly initialized
      implicit double precision (a-h,o-z)
      double complex c(*)
      double precision w(0:*)
      integer tables(-14:*)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
*  loading the common block : constants
      ivs=ivsx
      ies=iesx
      nv=nvx
      ne=nex
      ilim=(ne-1)*ies
      lstep=ivs
*  loading the common block : iteration-dependent quantities: initializa
      lx=1
      mx=nv
*  select the highest factor of nv
      ifac=tables(-1)
      if(ifac > 3)then
        ierr=tberr
        return
      endif
*..  radix 2 loop
200   continue
      do 210 im=1,tables(-14)
        mx=mx/2
        nustep=lx*lstep
        llim=nustep-lstep
        ld2lim=nustep/2-lstep
        mstep=nustep*2
        mlim=nv*lstep-mstep
        call mffta5(c,w)
        lx=lx+lx
210   continue
      if(ifac == 1)return
*..  radix 3 loop
300   continue
      do 310 im=1,tables(-13)
        mx=mx/3
        nustep=lx*lstep
        llim=nustep-lstep
        mstep=nustep*3
        mlim=nv*lstep-mstep
        call mfftb5(c,w)
        lx=lx*3
310   continue
      if(ifac == 2)return
*..   radix 5 loop
500   continue
      do 510 im=1,tables(-12)
        mx=mx/5
        nustep=lx*lstep
        llim=nustep-lstep
        mstep=nustep*5
        mlim=nv*lstep-mstep
        call mfftc5(c,w)
        lx=lx*5
510   continue
      end
      subroutine mfftom(c,ims,ivs,ies,nm,nv,ne,index,itemp)
*
*   purpose :
*     this routine performs a reordering of a vector-of-2 vectors
*     of complex c[ims,nm [ivs,nv [ies,ne]]], according to a
*     permutation index "index".
*     see ref.[1] for notations, and comments to mfftdm above.
*  arguments
*     c : vector-of-2 vectors to be reordered
*     ims,ivs,ies,nm,nv,ne: these arguments describe the structure of
*     c, according to the above definition;
*     index: integer array , containing the permutation index  ;
*            it is nm elements long ; prepared by mfftp.
*     iwork: integer array, of length at least nm, used as workspace;
*
      implicit double precision (a-h,o-z)
      integer index(0:*),itemp(0:*)
      double complex c(0:ims-1,0:*),t

      i3lim=(nv-1)*ivs
      jlim=(ne-1)*ies
      do 1 i=0,nm-1
1     itemp(i)=index(i)

      do 4 i=1,nm-3
2       if(itemp(i) /= i)then
          idest=itemp(i)
          do 3 i3=0,i3lim,ivs
          do 3 j=i3,i3+jlim,ies
            t=c(j,i)
            c(j,i)=c(j,idest)
            c(j,idest)=t
3         continue
          itemp(i)=itemp(idest)
          itemp(idest)=idest
        goto 2
        endif
4     continue
      end
      subroutine mfftov(c,ivs,ies,nv,ne,index,itemp)
*
*   purpose :
*     this routine performs a reordering of a vector-of-vectors
*     of complex c[ivs,nv [ies,ne]], according to a
*     permutation index "index".
*     see ref.[1] for notations, and comments to mfftdv.
*  arguments
*     c : vector-of-vectors to be reordered
*     ivs,ies,nv,ne: these arguments describe the structure of
*     c, according to the above definition;
*     index: integer array , containing the permutation index ;
*            it is nv elements long; prepared by mfftp.
*     iwork: integer array, of length at least nv, used as workspace;
*
      implicit double precision (a-h,o-z)
      integer index(0:nv-1),itemp(0:nv-1)
      double complex c(ivs,0:*),t

      neies=ne*ies
      do 1 i=0,nv-1
1     itemp(i)=index(i)

      do 4 i=1,nv-3
2       if(itemp(i) /= i)then
          idest=itemp(i)
          do 3 j=1,neies,ies
            t=c(j,i)
            c(j,i)=c(j,idest)
            c(j,idest)=t
3         continue
          itemp(i)=itemp(idest)
          itemp(idest)=idest
        goto 2
        endif
4     continue
      end
      subroutine mfftp(n,w,iw,id,ierr)
*
*     this subroutine prepares tables for use by the mfft
*     routines. the parameters are
*
*     n : is the order of the transform;
*
*     w : array of length at least 4*n+14  words if id=0,
*      and 4*n*(id+1)+14 if id > 0 ; this array is filled
*      by the present routine, and should not be modified by the
*      user; it is required by all the operating routines;
*      the routines mfftis and mfftds
*      require that w has been filled by a call to mfftp with
*      id > 0 ; all the routines do not modify the contents
*      of w ;
*      warning: different portions of w are handled as complex or
*               integer variables by different routines.
*
*     id : if id == 0 the tables are set for a normal
*      transform; if id > 0 the tables are set for both
*      a normal and a "special" transform (i.e. optimized for small
*      first dimension data arrays, see ref.[2]); in this case it
*      should be equal to the first dimension of the array to be
*      transformed, as declared in the calling program.
*
*     ierr : error code : =0 : successful
*                       : =2 : factorization error
*
*
*************************************************************
*
*      reference information : layout of w
*
*  in all cases
*
* word address      type    n. of el. word length   purpose
*
*  0                integer      14      14     factorization of n
*
*  14               complex     n       2*n     exp(i*p/(2*n)*k),k=
*                                               0,n-1
*  3*n+14           integer     n        n      permutation index:kofi
*  2*n+14           integer     n        n      permutation index:iofk
*
*  only if id > 0
*
*  4*n+14           complex    n*id     2*n*id  tables for mfftds
*  4*n+2*n*id+14    complex    n*id     2*n*id  tables for mfftis
*
*
*
*******************************************************************
      implicit double precision (a-h,o-z)
      integer iw(-14:*)
      double precision w(0:1,0:*)

*...  factorization of n in iw(-14)..w(-1)
      call mfftp1(iw(-14),n,ierr)
      if(ierr /= 0)return
*
*     preparation of permutation indexes in w(n)..w(2*n-1)
*     warning : iw(0)..iw(n-1) used as a work space
      call mfftp2(w(2*n,0),w(3*n,0),w(0,0),iw(-14),n)
*
*     preparation of phase factor table in w(0)..w(2n-1)
      w(0,0)=1.0d0
      w(1,0)=0.0d0
      pi2dn=datan(1.d0)*8.d0/n
      do 1  i=1,n-1
        w(0,i) = dcos(pi2dn*i)
        w(1,i) = dsin(pi2dn*i)
1     continue
*
*     if tables for special transform are requested
      if(id > 0)then
        call mfftp4(w(0,0),w(4*n,0),iw(-14),n,id)
      endif
*
      end
      subroutine mfftrd(c,isv,ise,nv,ne,rw)
*
*   purpose:
*       this routine performs the post-processing phase for
*       real 2-dimensional dft's, according to formula (2.7)
*       in ref.[1].
*       post-processing acts after computing the complex dft
*       and eventual reordering (calls to mfftdv and mfftov).
*       it applies to a vector-of-vectors-of-complex
*               c[ivs,nv [ies,ne]].
*       see ref.[1] for notations.
*
*   arguments:
*      input:
*       c : data array, output from mfftdv, mfftov; to be declared
*                      real c(ise*2,ne)
*           in the calling program.
*     isv : separation of elements in a column of c (usually 1)
*     ise : separation of elements in a row of c, divided by 2
*      nv : no. of elements to be processed in a column of c
*      ne : no. of elements in a row of c, divided by 2.
*      rw : complex array of lenght at least nv; it must
*           be initialized by a call to mfftrp; it remains
*           unchanged in output.
*
*  output : post-processed array c
*
*
      implicit double precision (a-h,o-z)
      double complex c(0:isv-1,0:*),rw(0:*),t1,t2
*
      if (nv > 1) then
      do 200 iv=1,(nv-1)/2
       do 190 ie=0,(ne-1)*ise,ise
*
      t1=c(ie,iv)
      t2=c(ie,nv-iv)
      c(ie,iv)=((t1+dconjg(t2))+(rw(iv)*(t1-dconjg(t2))))*0.5d0
      c(ie,nv-iv)=(dconjg(t1+dconjg(t2))-dconjg(rw(iv)*(t1-dconjg(t2))))
     .  *0.5d0
*
 190   continue
 200  continue
*
*
      if(2*iv == nv) then
        do 210 ie=0,(ne-1)*ise,ise
        c(ie,iv)=dconjg(c(ie,iv))
 210  continue
*
      endif
*
       endif

      do 300 ie=0,(ne-1)*ise,ise
        t1=c(ie,0)
        c(ie,0)=(dble(t1)+dimag(t1))
        c(ie,nv)=(dble(t1)-dimag(t1))
 300  continue
*
*
      end
      subroutine mfftri(c,isv,ise,nv,ne,rw)
*
*   purpose:
*       this routine performs the pre-processing phase for
*       real 2-dimensional inverse dft's (see formula (2.7)
*       in ref.[1]).
*       pre-processing acts on sequental data before a call to
*       mfftov an before computing the idft (a call to mfftiv)
*       and eventual reordering (calls to mfftiv and mfftov).
*       it applies to a vector-of-vectors-of-complex
*               c[ivs,nv [ies,ne]].
*       see ref.[1] for notations.
*
*   arguments:
*      input:
*       c : data array, output from mfftdv, mfftov; to be declared
*                      real c(ise*2,ne)
*           in the calling program.
*     isv : separation of elements in a column of c (usually 1)
*     ise : separation of elements in a row of c, divided by 2
*      nv : no. of elements to be processed in a column of c
*      ne : no. of elements in a row of c, divided by 2.
*      rw : complex array of lenght at least nv; it must
*           be initialized by a call to mfftrp; it remains
*           unchanged in output.
*
*  output : post-processed array c
*
*
      implicit double precision (a-h,o-z)
      double complex c(0:isv-1,0:*),rw(0:*),t1,t2
*
*
      if(nv > 1) then
*
      do 200 iv=1,(nv-1)/2
        do 190 ie=0,(ne-1)*ise,ise
        t1=c(ie,iv)
        t2=c(ie,nv-iv)
*
        c(ie,iv)=(t1+dconjg(t2))+(dconjg(rw(iv))*(t1-dconjg(t2)))
        c(ie,nv-iv)=dconjg(t1+dconjg(t2))-dconjg(dconjg(rw(iv))*(t1-
     $              dconjg(t2)))


 190    continue
 200   continue
*
*
        if(2*iv == nv) then
*
        do 210 ie=0,(ne-1)*ise,ise
        c(ie,iv)=2*dconjg(c(ie,iv))
 210   continue
*
*
      endif
*
      endif
*
*
      do 220 ie=0,(ne-1)*ise,ise
         rp=dble(c(ie,0))+dble(c(ie,nv))
         rm=dble(c(ie,0))-dble(c(ie,nv))
         c(ie,0)=dcmplx(rp,rm)
 220  continue
*
*
      end

      subroutine mfftrp(nn,rw)
*
*  this routine prepares the phase factor tables to be used in the
*  post-processing phase in real transforms, i.e. r2fft and r3fft.
*
*  the definition of phase factors we use is the following
*
*             rw(iv) = -i*exp(-i*2*pi/nn*iv)
*
*  where the symbols are:
*
*  i: imaginary unit
*  nn : number of real elements
*  iv: running index from 0 to nn-1
*
*
*
*
      implicit double precision (a-h,o-z)
      double complex rw(0:*)
*
*
      pi2=8*datan(1.d0)
      phase=pi2/nn
      rw(0)=(0.d0,-1.d0)
*
*
      do 10 i=1,nn-1
      rw(i)=dcmplx(-dsin(i*phase),-dcos(i*phase))
 10   continue
*
*
      end
      subroutine mffta4(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 2 step applied to a vector-of
*       vectors-of-complex c[ivs,nv [ies,ne]]. see ref.[1] for notations
*       this routine can be used only by routine mfftdv, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*)
      double complex t0,f


        if(mx > 2*lx)then
          do 100 lam=0,llim,lstep
            muf=lx
            do 90 mu=lam+mstep,lam+md2lim,mstep
              f=dconjg(fac(muf))
              do 80 i=mu,mu+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=(t0-c(i,1))*f
80            continue
              muf=muf+lx
90          continue
            muf=muf-lx
            do 91 mu=lam+md2lim+2*mstep,lam+mlim,mstep
              f=-fac(muf)
              do 81 i=mu,mu+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=(t0-c(i,1))*f
81            continue
              muf=muf-lx
91          continue
            do 82 i=md2lim+mstep+lam,md2lim+mstep+lam+ilim,ies
               t0=c(i,0)
              c(i,0)=t0+c(i,1)
              c(i,1)=dcmplx(dimag(t0-c(i,1)),-dble(t0-c(i,1)))
82          continue
            do 83 i=lam,lam+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=t0-c(i,1)
83          continue
100       continue
        else
          if(mx == 1)goto 1000
* if mx > 1 come here
            muf=lx
            do 200 mu=mstep,md2lim,mstep
              f=dconjg(fac(muf))
              do 190 lam=mu,mu+llim,lstep
                do 180 i=lam,lam+ilim,ies
                  t0=c(i,0)
                  c(i,0)=t0+c(i,1)
                  c(i,1)=(t0-c(i,1))*f
  180           continue
  190         continue
              muf=muf+lx
200         continue
            muf=muf-lx
            do 201 mu=md2lim+2*mstep,mlim,mstep
              f=-fac(muf)
              do 191 lam=mu,mu+llim,lstep
                do 181 i=lam,lam+ilim,ies
                  t0=c(i,0)
                  c(i,0)=t0+c(i,1)
                  c(i,1)=(t0-c(i,1))*f
  181           continue
  191         continue
              muf=muf-lx
201         continue
             do 192 lam=md2lim+mstep,md2lim+mstep+llim,lstep
             do 182 i=lam,lam+ilim,ies
               t0=c(i,0)
               c(i,0)=t0+c(i,1)
               c(i,1)=dcmplx(dimag(t0-c(i,1)),-dble(t0-c(i,1)))
182           continue
192         continue
1000        do 193 lam=0,llim,lstep
              do 183 i=lam,lam+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=t0-c(i,1)
183           continue
193         continue
        endif
      end
      subroutine mffta5(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 2 step applied to a vector-of
*       vectors-of-complex c[ivs,nv [ies,ne]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftiv, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:*),t0,f
      double complex fac(0:*)

        if(2*mx >= lx)then
          if(lx == 1)goto 1000
          lamf=mx
          do 100 lam=lstep,ld2lim,lstep
            f=fac(lamf)
            do 90 mu=lam,lam+mlim,mstep
              do 80 i=mu,mu+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
80            continue
90          continue
            lamf=lamf+mx
100       continue
          lamf=lamf-mx
          do 101 lam=ld2lim+2*lstep,llim,lstep
            f=-dconjg(fac(lamf))
            do 91 mu=lam,lam+mlim,mstep
              do 81 i=mu,mu+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
81            continue
91          continue
            lamf=lamf-mx
101       continue
          do 93 mu=ld2lim+lstep,ld2lim+lstep+mlim,mstep

            do 83 i=mu,mu+ilim,ies
              t0=dcmplx(-dimag(c(i,1)),dble(c(i,1)))
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
83          continue
93        continue
1000      do 92 mu=0,mlim,mstep

            do 82 i=mu,mu+ilim,ies
              t0=c(i,1)
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
82          continue
92        continue
        else
          do 200 mu=0,mlim,mstep
            lamf=mx
            do 190 lam=mu+lstep,mu+ld2lim,lstep
              f=fac(lamf)
              do 180 i=lam,lam+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
180           continue
              lamf=lamf+mx
190         continue
            lamf=lamf-mx
            do 191 lam=mu+ld2lim+2*lstep,mu+llim,lstep
              f=-dconjg(fac(lamf))
              do 181 i=lam,lam+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
181           continue
              lamf=lamf-mx
191         continue
            do 182 i=mu,mu+ilim,ies
              t0=c(i,1)
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
182         continue
            do 183 i=mu+ld2lim+lstep,mu+ld2lim+lstep+ilim,ies
              t0=dcmplx(-dimag(c(i,1)),dble(c(i,1)))
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
183         continue
200       continue
        endif
      end
      subroutine mffta6(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 2 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdm, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*)
      double complex t0,f


        if(mx > 2*lx)then
          do 100 lam=0,llim,lstep
            muf=lx
            do 90 mu=lam+mstep,lam+md2lim,mstep
              f=dconjg(fac(muf))
              do 80  iv=mu,mu+ivlim,ivs
              do 80 i=iv,iv+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=(t0-c(i,1))*f
80            continue
              muf=muf+lx
90          continue
            muf=muf-lx
            do 91 mu=lam+md2lim+2*mstep,lam+mlim,mstep
              f=-fac(muf)
              do 81 iv=mu,mu+ivlim,ivs
              do 81 i=iv,iv+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=(t0-c(i,1))*f
81            continue
              muf=muf-lx
91          continue
            do 82 iv=md2lim+mstep+lam,md2lim+mstep+lam+ivlim,ivs
            do 82 i=iv,iv+ilim,ies
               t0=c(i,0)
              c(i,0)=t0+c(i,1)
              c(i,1)=dcmplx(dimag(t0-c(i,1)),-dble(t0-c(i,1)))
82          continue
            do 83 iv=lam,lam+ivlim,ivs
            do 83 i=iv,iv+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=t0-c(i,1)
83          continue
100       continue
        else
          if(mx == 1)goto 1000
* if mx > 1 come here
            muf=lx
            do 200 mu=mstep,md2lim,mstep
              f=dconjg(fac(muf))
              do 190 lam=mu,mu+llim,lstep
                do 180 iv=lam,lam+ivlim,ivs
                do 180 i=iv,iv+ilim,ies
                  t0=c(i,0)
                  c(i,0)=t0+c(i,1)
                  c(i,1)=(t0-c(i,1))*f
  180           continue
  190         continue
              muf=muf+lx
200         continue
            muf=muf-lx
            do 201 mu=md2lim+2*mstep,mlim,mstep
              f=-fac(muf)
              do 191 lam=mu,mu+llim,lstep
                do 181 iv=lam,lam+ivlim,ivs
                do 181 i=iv,iv+ilim,ies
                  t0=c(i,0)
                  c(i,0)=t0+c(i,1)
                  c(i,1)=(t0-c(i,1))*f
  181           continue
  191         continue
              muf=muf-lx
201         continue
             do 192 lam=md2lim+mstep,md2lim+mstep+llim,lstep
             do 182 iv=lam,lam+ivlim,ivs
             do 182 i=iv,iv+ilim,ies
               t0=c(i,0)
               c(i,0)=t0+c(i,1)
               c(i,1)=dcmplx(dimag(t0-c(i,1)),-dble(t0-c(i,1)))
182           continue
192         continue
1000        do 193 lam=0,llim,lstep
              do 183 iv=lam,lam+ivlim,ivs
              do 183 i=iv,iv+ilim,ies
                t0=c(i,0)
                c(i,0)=t0+c(i,1)
                c(i,1)=t0-c(i,1)
183           continue
193         continue
        endif
      end
      subroutine mffta7(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 2 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftim, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*)
      double complex t0,f

        if(2*mx > lx)then
          if(lx == 1)goto 1000
*   else here
          lamf=mx
          do 100 lam=lstep,ld2lim,lstep
            f=fac(lamf)
            do 90 mu=lam,lam+mlim,mstep
              do 80 iv=mu,mu+ivlim,ivs
              do 80 i=iv,iv+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
80            continue
90          continue
            lamf=lamf+mx
100       continue
          lamf=lamf-mx
          do 101 lam=ld2lim+2*lstep,llim,lstep
            f=-dconjg(fac(lamf))
            do 91 mu=lam,lam+mlim,mstep
              do 81 iv=mu,mu+ivlim,ivs
              do 81 i=iv,iv+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
81            continue
91          continue
            lamf=lamf-mx
101       continue
          do 93 mu=ld2lim+lstep,ld2lim+lstep+mlim,mstep

            do 83 iv=mu,mu+ivlim,ivs
            do 83 i=iv,iv+ilim,ies
              t0=dcmplx(-dimag(c(i,1)),dble(c(i,1)))
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
83          continue
93        continue
1000      do 92 mu=0,mlim,mstep

            do 82 iv=mu,mu+ivlim,ivs
            do 82 i=iv,iv+ilim,ies
              t0=c(i,1)
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
82          continue
92        continue
        else
          do 200 mu=0,mlim,mstep
            lamf=mx
            do 190 lam=mu+lstep,mu+ld2lim,lstep
              f=fac(lamf)
              do 180 iv=lam,lam+ivlim,ivs
              do 180 i=iv,iv+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
180           continue
              lamf=lamf+mx
190         continue
            lamf=lamf-mx
            do 191 lam=mu+ld2lim+2*lstep,mu+llim,lstep
              f=-dconjg(fac(lamf))
              do 181 iv=lam,lam+ivlim,ivs
              do 181 i=iv,iv+ilim,ies
                t0=c(i,1)*f
                c(i,1)=c(i,0)-t0
                c(i,0)=c(i,0)+t0
181           continue
              lamf=lamf-mx
191         continue
            do 182 iv=mu,mu+ivlim,ivs
            do 182 i=iv,iv+ilim,ies
              t0=c(i,1)
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
182         continue
            do 183 iv=mu+ld2lim+lstep,mu+ld2lim+lstep+ivlim,ivs
            do 183 i=iv,iv+ilim,ies
              t0=dcmplx(-dimag(c(i,1)),dble(c(i,1)))
              c(i,1)=c(i,0)-t0
              c(i,0)=c(i,0)+t0
183         continue
200       continue
        endif
      end

      subroutine mffta8(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix-2 step applied to a vector-of-
*       2-vectors-of-complex [ims,nm [ivs,nv [ies,ne]]], optimized for
*       small ne matrices.
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdm, which controls
*       its operation through common mfftpa.
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*),t0
*
      if (mx /= 1) then
*
      do 200 lam=0,llim,lstep
        do 150 iv=lam,lam+ivlim,ivs
          imuf=0
          do 100 imu=iv,iv+ilim
            t0=c(imu,0)
            c(imu,0)=t0+c(imu,1)
            c(imu,1)=(t0-c(imu,1))*fac(imuf)
            imuf=imuf+1
100       continue
150      continue
200   continue
*
      else
        do 400 lam=0,llim,lstep
          do 350 iv=lam,lam+ivlim,ivs
            do 300 imu=iv,iv+ilim
            t0=c(imu,0)
            c(imu,0)=t0+c(imu,1)
            c(imu,1)=t0-c(imu,1)
 300        continue
 350      continue
 400    continue
      endif
*
      end
      subroutine mffta9(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix-2 step applied to a vector-of-
*       2-vectors-of-complex [ims,nm [ivs,nv [ies,ne]]], optimized for
*       small ne matrices.
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftim, which controls
*       its operation through common mfftpa.
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*),t0
*
      if(lx /= 1) then
*
      do 200 mu=0,mlim,mstep
        do 150 iv=mu,mu+ivlim,ivs
          ilamf=0
          do 100 ilam=iv,iv+ilim
            t0=c(ilam,1)*fac(ilamf)
            c(ilam,1)=c(ilam,0)-t0
            c(ilam,0)=c(ilam,0)+t0
            ilamf=ilamf+1
100       continue
150     continue
200   continue
*
      else
      do 400 mu=0,mlim,mstep
        do 350 iv=mu,mu+ivlim,ivs
          do 300 ilam=iv,iv+ilim
          t0=c(ilam,1)
          c(ilam,1)=c(ilam,0)-t0
          c(ilam,0)=c(ilam,0)+t0
 300      continue
 350    continue
 400  continue
      endif
*
      end
      subroutine mfftb4(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 3 step applied to a vector-of
*       vectors-of-complex c[ivs,nv [ies,ne]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdv, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:2),fac(0:*)
      double complex t0,t1,t2,f1,f2

      double precision sin60
      parameter ( sin60 =  8.6602540378443864d-1)

*..  mu > 1
            muf=lx
            do 200 mu=mstep,mlim,mstep
              f1=dconjg(fac(muf))
              f2=dconjg(fac(2*muf))
              do 190 lam=mu,mu+llim,lstep
                do 180 i=lam,lam+ilim,ies
                t0=c(i,1)+c(i,2)
                t1=c(i,0)-0.5d0*t0
                t2=(c(i,1)-c(i,2))*sin60
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1-dcmplx(-dimag(t2),dble(t2)))*f1
                c(i,2)=(t1+dcmplx(-dimag(t2),dble(t2)))*f2
  180           continue
  190         continue
              muf=muf+lx
200         continue

*..  mu=0
1000        do 193 lam=0,llim,lstep
              do 183 i=lam,lam+ilim,ies
                t0=c(i,1)+c(i,2)
                t1=c(i,0)-0.5d0*t0
                t2=(c(i,1)-c(i,2))*sin60
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1-dcmplx(-dimag(t2),dble(t2)))
                c(i,2)=(t1+dcmplx(-dimag(t2),dble(t2)))
183           continue
193         continue
      end
      subroutine mfftb5(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 3 step applied to a vector-of
*       vectors-of-complex c[ivs,nv [ies,ne]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftiv, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:2),fac(0:*)
      double complex t0,t1,t2,f1,f2
      double precision sin60
      parameter ( sin60 =  8.6602540378443864d-1)

*..  lam > 0
          lamf=mx
          do 100 lam=lstep,llim,lstep
            f1=fac(lamf)
            f2=fac(2*lamf)
            do 90 mu=lam,lam+mlim,mstep
              do 80 i=mu,mu+ilim,ies
                t0=c(i,1)*f1+c(i,2)*f2
                t2=(c(i,1)*f1-c(i,2)*f2)*sin60
                t1=c(i,0)-0.5d0*t0
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1+dcmplx(-dimag(t2),dble(t2)))
                c(i,2)=(t1-dcmplx(-dimag(t2),dble(t2)))
80            continue
90          continue
            lamf=lamf+mx
100       continue

*.. lam=0
          do 92 mu=0,mlim,mstep
            do 82 i=mu,mu+ilim,ies
                t0=c(i,1)+c(i,2)
                t2=(c(i,1)-c(i,2))*sin60
                t1=c(i,0)-0.5d0*t0
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1+dcmplx(-dimag(t2),dble(t2)))
                c(i,2)=(t1-dcmplx(-dimag(t2),dble(t2)))
82          continue
92        continue
      end
      subroutine mfftb6(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 3 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdm, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:2),fac(0:*)
      double complex t0,t1,t2,f1,f2
      double precision sin60
      parameter ( sin60 =  8.6602540378443864d-1)


*..  mu > 0
            muf=lx
            do 200 mu=mstep,mlim,mstep
              f1=dconjg(fac(muf))
              f2=dconjg(fac(muf*2))
              do 190 lam=mu,mu+llim,lstep
                do 180 iv=lam,lam+ivlim,ivs
                do 180 i=iv,iv+ilim,ies
                  t0=c(i,1)+c(i,2)
                  t2=(c(i,1)-c(i,2))*sin60
                  t1=c(i,0)-0.5d0*t0
                  c(i,0)=c(i,0)+t0
                  c(i,1)=(t1-dcmplx(-dimag(t2),dble(t2)))*f1
                  c(i,2)=(t1+dcmplx(-dimag(t2),dble(t2)))*f2
  180           continue
  190         continue
              muf=muf+lx
200         continue
            do 193 lam=0,llim,lstep
              do 183 iv=lam,lam+ivlim,ivs
              do 183 i=iv,iv+ilim,ies
                t0=c(i,1)+c(i,2)
                t2=(c(i,1)-c(i,2))*sin60
                t1=c(i,0)-0.5d0*t0
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1-dcmplx(-dimag(t2),dble(t2)))
                c(i,2)=(t1+dcmplx(-dimag(t2),dble(t2)))
183           continue
193         continue
      end
      subroutine mfftb7(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 3 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftim, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:2),fac(0:*)
      double complex t0,t1,t2,f1,f2
      double precision sin60
      parameter ( sin60 =  8.6602540378443864d-1)
*..  lam>0
          lamf=mx
          do 100 lam=lstep,llim,lstep
            f1=fac(lamf)
            f2=fac(lamf*2)
            do 90 mu=lam,lam+mlim,mstep
              do 80 iv=mu,mu+ivlim,ivs
              do 80 i=iv,iv+ilim,ies
                t0=c(i,1)*f1+c(i,2)*f2
                t2=(c(i,1)*f1-c(i,2)*f2)*sin60
                t1=c(i,0)-0.5d0*t0
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1+dcmplx(-dimag(t2),dble(t2)))
                c(i,2)=(t1-dcmplx(-dimag(t2),dble(t2)))
80            continue
90          continue
            lamf=lamf+mx
100       continue
          do 92 mu=0,mlim,mstep
            do 82 iv=mu,mu+ivlim,ivs
            do 82 i=iv,iv+ilim,ies
                t0=c(i,1)+c(i,2)
                t2=(c(i,1)-c(i,2))*sin60
                t1=c(i,0)-0.5d0*t0
                c(i,0)=c(i,0)+t0
                c(i,1)=(t1+dcmplx(-dimag(t2),dble(t2)))
                c(i,2)=(t1-dcmplx(-dimag(t2),dble(t2)))
82          continue
92        continue
      end
      subroutine mfftb8(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 3 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdm, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*)
      double complex t0,t1,t2,f1,f2
      double precision sin60
      parameter ( sin60 =  8.6602540378443864d-1)


      if(mx /= 1)then
            do 200 lam=0,llim,lstep
              do 190 iv=lam,lam+ivlim,ivs
                imuf=0
                do 180 imu=iv,iv+ilim

                  t0=c(imu,1)+c(imu,2)
                  t2=(c(imu,1)-c(imu,2))*sin60
                  t1=c(imu,0)-0.5d0*t0
                  c(imu,0)=c(imu,0)+t0
                  c(imu,1)=(t1-dcmplx(-dimag(t2),dble(t2)))*fac(imuf)
                  c(imu,2)=(t1+dcmplx(-dimag(t2),dble(t2)))*fac(imuf+
     $                     nustep)
                  imuf=imuf+1
  180           continue
  190         continue
200         continue
        else
            do 400 lam=0,llim,lstep
              do 390 iv=lam,lam+ivlim,ivs
                do 380 imu=iv,iv+ilim
                  t0=c(imu,1)+c(imu,2)
                  t2=(c(imu,1)-c(imu,2))*sin60
                  t1=c(imu,0)-0.5d0*t0
                  c(imu,0)=c(imu,0)+t0
                  c(imu,1)=(t1-dcmplx(-dimag(t2),dble(t2)))
                  c(imu,2)=(t1+dcmplx(-dimag(t2),dble(t2)))
380             continue
390           continue
400         continue
        endif
      end
      subroutine mfftb9(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 3 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftim, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:1),fac(0:*)
      double complex t0,t1,t2,f1,f2
      double precision sin60
      parameter ( sin60 =  8.6602540378443864d-1)

      if(lx /= 1)then
          do 200 mu=0,mlim,mstep
            do 150 iv=mu,mu+ivlim,ivs
              ilamf=0
              do 100 ilam=iv,iv+ilim
                t0=c(ilam,1)*fac(ilamf)+c(ilam,2)*fac(ilamf+nustep)
                t2=(c(ilam,1)*fac(ilamf)-c(ilam,2)*fac(ilamf+nustep))*
     $             sin60
                t1=c(ilam,0)-0.5d0*t0
                c(ilam,0)=c(ilam,0)+t0
                c(ilam,1)=(t1+dcmplx(-dimag(t2),dble(t2)))
                c(ilam,2)=(t1-dcmplx(-dimag(t2),dble(t2)))
                ilamf=ilamf+1
100           continue
150         continue
200       continue
      else
          do 400 mu=0,mlim,mstep
            do 350 iv=mu,mu+ivlim,ivs
              do 300 ilam=iv,iv+ilim
                t0=c(ilam,1)+c(ilam,2)
                t2=(c(ilam,1)-c(ilam,2))*sin60
                t1=c(ilam,0)-0.5d0*t0
                c(ilam,0)=c(ilam,0)+t0
                c(ilam,1)=(t1+dcmplx(-dimag(t2),dble(t2)))
                c(ilam,2)=(t1-dcmplx(-dimag(t2),dble(t2)))
300           continue
350         continue
400       continue
      endif
      end
      subroutine mfftc4(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 5 step applied to a vector-of
*       vectors-of-complex c[ivs,nve[ies,ne]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdv, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:4),fac(0:*)
      double complex t1,t2,t3,t4,t5,f1,f2,f3,f4
      double precision sin72,rad5d4,s36d72
      parameter (
     $ sin72 =  9.51056516295153572116439333d-1,
     $ rad5d4 =  5.59016994374947424102293417d-1,
     $ s36d72 =  6.18033988749894848204586834d-1 )
*..  mu > 1
            muf=lx
            do 200 mu=mstep,mlim,mstep
              f1=dconjg(fac(muf))
              f2=dconjg(fac(2*muf))
              f3=dconjg(fac(3*muf))
              f4=dconjg(fac(4*muf))
              do 190 lam=mu,mu+llim,lstep
                do 180 i=lam,lam+ilim,ies
                  t1=c(i,1)+c(i,4)
                  t2=c(i,2)+c(i,3)
                  t3=(c(i,1)-c(i,4))*sin72
                  t4=(c(i,2)-c(i,3))*sin72
                  t5=t1+t2
                  t1=rad5d4*(t1-t2)
                  t2=c(i,0)-0.25d0*t5
                  c(i,0)=c(i,0)+t5
                  t5=t2+t1
                  t2=t2-t1
                  t1=t3+s36d72*t4
                  t3=s36d72*t3-t4
                  c(i,1)=(t5-dcmplx(-dimag(t1),dble(t1)))*f1
                  c(i,4)=(t5+dcmplx(-dimag(t1),dble(t1)))*f4
                  c(i,2)=(t2-dcmplx(-dimag(t3),dble(t3)))*f2
                  c(i,3)=(t2+dcmplx(-dimag(t3),dble(t3)))*f3
  180           continue
  190         continue
              muf=muf+lx
200         continue
*..  mu=0
1000        do 193 lam=0,llim,lstep
              do 183 i=lam,lam+ilim,ies
                t1=c(i,1)+c(i,4)
                t2=c(i,2)+c(i,3)
                t3=(c(i,1)-c(i,4))*sin72
                t4=(c(i,2)-c(i,3))*sin72
                t5=t1+t2
                t1=rad5d4*(t1-t2)
                t2=c(i,0)-0.25d0*t5
                c(i,0)=c(i,0)+t5
                t5=t2+t1
                t2=t2-t1
                t1=t3+s36d72*t4
                t3=s36d72*t3-t4
                c(i,1)=t5-dcmplx(-dimag(t1),dble(t1))
                c(i,4)=t5+dcmplx(-dimag(t1),dble(t1))
                c(i,2)=t2-dcmplx(-dimag(t3),dble(t3))
                c(i,3)=t2+dcmplx(-dimag(t3),dble(t3))
183           continue
193         continue
      end
      subroutine mfftc5(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 5 step applied to a vector-of
*       vectors-of-complex c[ivs,nv [ies,ne]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftiv, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:4),fac(0:*)
      double complex t1,t2,t3,t4,t5,f1,f2,f3,f4
      double precision sint2,rad5d4,s36d72
      parameter (
     $ sin72 =  9.51056516295153572116439333d-1,
     $ rad5d4 =  5.59016994374947424102293417d-1,
     $ s36d72 =  6.18033988749894848204586834d-1 )
*..  lam > 0
          lamf=mx
          do 100 lam=lstep,llim,lstep
            f1=fac(lamf)
            f2=fac(2*lamf)
            f3=fac(3*lamf)
            f4=fac(4*lamf)
            do 90 mu=lam,lam+mlim,mstep
              do 80 i=mu,mu+ilim,ies
                t1=c(i,1)*f1+c(i,4)*f4
                t2=c(i,2)*f2+c(i,3)*f3
                t3=(c(i,1)*f1-c(i,4)*f4)*sin72
                t4=(c(i,2)*f2-c(i,3)*f3)*sin72
                t5=t1+t2
                t1=rad5d4*(t1-t2)
                t2=c(i,0)-0.25d0*t5
                c(i,0)=c(i,0)+t5
                t5=t2+t1
                t2=t2-t1
                t1=t3+s36d72*t4
                t3=s36d72*t3-t4
                c(i,1)=t5+dcmplx(-dimag(t1),dble(t1))
                c(i,4)=t5-dcmplx(-dimag(t1),dble(t1))
                c(i,2)=t2+dcmplx(-dimag(t3),dble(t3))
                c(i,3)=t2-dcmplx(-dimag(t3),dble(t3))
80            continue
90          continue
            lamf=lamf+mx
100       continue

*.. lam=0
          do 92 mu=0,mlim,mstep
            do 82 i=mu,mu+ilim,ies
              t1=c(i,1)+c(i,4)
              t2=c(i,2)+c(i,3)
              t3=(c(i,1)-c(i,4))*sin72
              t4=(c(i,2)-c(i,3))*sin72
              t5=t1+t2
              t1=rad5d4*(t1-t2)
              t2=c(i,0)-0.25d0*t5
              c(i,0)=c(i,0)+t5
              t5=t2+t1
              t2=t2-t1
              t1=t3+s36d72*t4
              t3=s36d72*t3-t4
              c(i,1)=t5+dcmplx(-dimag(t1),dble(t1))
              c(i,4)=t5-dcmplx(-dimag(t1),dble(t1))
              c(i,2)=t2+dcmplx(-dimag(t3),dble(t3))
              c(i,3)=t2-dcmplx(-dimag(t3),dble(t3))
82          continue
92        continue
      end
c     ###################    fft5 ends here    ##################

      subroutine mfftc6(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 5 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftdm, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:4),fac(0:*)
      double complex t1,t2,t3,t4,t5,f1,f2,f3,f4
      double precision sin72,rad5d4,s36d72
      parameter (
     $ sin72 =  9.51056516295153572116439333d-1,
     $ rad5d4 =  5.59016994374947424102293417d-1,
     $ s36d72 =  6.18033988749894848204586834d-1 )
*..  mu > 0
            muf=lx
            do 200 mu=mstep,mlim,mstep
              f1=dconjg(fac(muf))
              f2=dconjg(fac(muf*2))
              f3=dconjg(fac(muf*3))
              f4=dconjg(fac(muf*4))
              do 190 lam=mu,mu+llim,lstep
                do 180 iv=lam,lam+ivlim,ivs
                do 180 i=iv,iv+ilim,ies
                  t1=c(i,1)+c(i,4)
                  t2=c(i,2)+c(i,3)
                  t3=(c(i,1)-c(i,4))*sin72
                  t4=(c(i,2)-c(i,3))*sin72
                  t5=t1+t2
                  t1=rad5d4*(t1-t2)
                  t2=c(i,0)-0.25d0*t5
                  c(i,0)=c(i,0)+t5
                  t5=t2+t1
                  t2=t2-t1
                  t1=t3+s36d72*t4
                  t3=s36d72*t3-t4
                  c(i,1)=(t5-dcmplx(-dimag(t1),dble(t1)))*f1
                  c(i,4)=(t5+dcmplx(-dimag(t1),dble(t1)))*f4
                  c(i,2)=(t2-dcmplx(-dimag(t3),dble(t3)))*f2
                  c(i,3)=(t2+dcmplx(-dimag(t3),dble(t3)))*f3
  180           continue
  190         continue
              muf=muf+lx
200         continue
            do 193 lam=0,llim,lstep
              do 183 iv=lam,lam+ivlim,ivs
              do 183 i=iv,iv+ilim,ies
                t1=c(i,1)+c(i,4)
                t2=c(i,2)+c(i,3)
                t3=(c(i,1)-c(i,4))*sin72
                t4=(c(i,2)-c(i,3))*sin72
                t5=t1+t2
                t1=rad5d4*(t1-t2)
                t2=c(i,0)-0.25d0*t5
                c(i,0)=c(i,0)+t5
                t5=t2+t1
                t2=t2-t1
                t1=t3+s36d72*t4
                t3=s36d72*t3-t4
                c(i,1)=t5-dcmplx(-dimag(t1),dble(t1))
                c(i,4)=t5+dcmplx(-dimag(t1),dble(t1))
                c(i,2)=t2-dcmplx(-dimag(t3),dble(t3))
                c(i,3)=t2+dcmplx(-dimag(t3),dble(t3))
183           continue
193         continue
      end
      subroutine mfftc7(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 5 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftim, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:4),fac(0:*)
      double complex t1,t2,t3,t4,t5,f1,f2,f3,f4
      double precision sin72,rad5d4,s36d72
      parameter (
     $ sin72 =  9.51056516295153572116439333d-1,
     $ rad5d4 =  5.59016994374947424102293417d-1,
     $ s36d72 =  6.18033988749894848204586834d-1 )
*..  lam>0
          lamf=mx
          do 100 lam=lstep,llim,lstep
            f1=fac(lamf)
            f2=fac(lamf*2)
            f3=fac(lamf*3)
            f4=fac(lamf*4)
            do 90 mu=lam,lam+mlim,mstep
              do 80 iv=mu,mu+ivlim,ivs
              do 80 i=iv,iv+ilim,ies
                t1=c(i,1)*f1+c(i,4)*f4
                t2=c(i,2)*f2+c(i,3)*f3
                t3=(c(i,1)*f1-c(i,4)*f4)*sin72
                t4=(c(i,2)*f2-c(i,3)*f3)*sin72
                t5=t1+t2
                t1=rad5d4*(t1-t2)
                t2=c(i,0)-0.25d0*t5
                c(i,0)=c(i,0)+t5
                t5=t2+t1
                t2=t2-t1
                t1=t3+s36d72*t4
                t3=s36d72*t3-t4
                c(i,1)=t5+dcmplx(-dimag(t1),dble(t1))
                c(i,4)=t5-dcmplx(-dimag(t1),dble(t1))
                c(i,2)=t2+dcmplx(-dimag(t3),dble(t3))
                c(i,3)=t2-dcmplx(-dimag(t3),dble(t3))
80            continue
90          continue
            lamf=lamf+mx
100       continue
          do 92 mu=0,mlim,mstep
            do 82 iv=mu,mu+ivlim,ivs
            do 82 i=iv,iv+ilim,ies
              t1=c(i,1)+c(i,4)
              t2=c(i,2)+c(i,3)
              t3=(c(i,1)-c(i,4))*sin72
              t4=(c(i,2)-c(i,3))*sin72
              t5=t1+t2
              t1=rad5d4*(t1-t2)
              t2=c(i,0)-0.25d0*t5
              c(i,0)=c(i,0)+t5
              t5=t2+t1
              t2=t2-t1
              t1=t3+s36d72*t4
              t3=s36d72*t3-t4
              c(i,1)=t5+dcmplx(-dimag(t1),dble(t1))
              c(i,4)=t5-dcmplx(-dimag(t1),dble(t1))
              c(i,2)=t2+dcmplx(-dimag(t3),dble(t3))
              c(i,3)=t2-dcmplx(-dimag(t3),dble(t3))
82          continue
92        continue
      end
      subroutine mfftc8(c,fac)
*
*   purpose:
*       elementary gentleman-sande radix 5 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftds, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:4),fac(0:*)
      double complex t1,t2,t3,t4,t5
      double precision sin72,rad5d4,s32d72
      parameter (
     $ sin72 =  9.51056516295153572116439333d-1,
     $ rad5d4 =  5.59016994374947424102293417d-1,
     $ s36d72 =  6.18033988749894848204586834d-1 )
      if(mx /= 1)then
            do 200 lam=0,llim,lstep
              do 190 iv=lam,lam+ivlim,ivs
                imuf=0
                imuf2=nustep
                imuf3=2*nustep
                imuf4=3*nustep
                do 180 imu=iv,iv+ilim
                  t1=c(imu,1)+c(imu,4)
                  t2=c(imu,2)+c(imu,3)
                  t3=(c(imu,1)-c(imu,4))*sin72
                  t4=(c(imu,2)-c(imu,3))*sin72
                  t5=t1+t2
                  t1=rad5d4*(t1-t2)
                  t2=c(imu,0)-0.25d0*t5
                  c(imu,0)=c(imu,0)+t5
                  t5=t2+t1
                  t2=t2-t1
                  t1=t3+s36d72*t4
                  t3=s36d72*t3-t4
                  c(imu,1)=(t5-dcmplx(-dimag(t1),dble(t1)))*fac(imuf )
                  c(imu,4)=(t5+dcmplx(-dimag(t1),dble(t1)))*fac(imuf4)
                  c(imu,2)=(t2-dcmplx(-dimag(t3),dble(t3)))*fac(imuf2)
                  c(imu,3)=(t2+dcmplx(-dimag(t3),dble(t3)))*fac(imuf3)
                  imuf=imuf+1
                  imuf2=imuf2+1
                  imuf3=imuf3+1
                  imuf4=imuf4+1
  180           continue
  190         continue
200         continue
        else
            do 400 lam=0,llim,lstep
              do 390 iv=lam,lam+ivlim,ivs
                do 380 imu=iv,iv+ilim
                  t1=c(imu,1)+c(imu,4)
                  t2=c(imu,2)+c(imu,3)
                  t3=(c(imu,1)-c(imu,4))*sin72
                  t4=(c(imu,2)-c(imu,3))*sin72
                  t5=t1+t2
                  t1=rad5d4*(t1-t2)
                  t2=c(imu,0)-0.25d0*t5
                  c(imu,0)=c(imu,0)+t5
                  t5=t2+t1
                  t2=t2-t1
                  t1=t3+s36d72*t4
                  t3=s36d72*t3-t4
                  c(imu,1)=t5-dcmplx(-dimag(t1),dble(t1))
                  c(imu,4)=t5+dcmplx(-dimag(t1),dble(t1))
                  c(imu,2)=t2-dcmplx(-dimag(t3),dble(t3))
                  c(imu,3)=t2+dcmplx(-dimag(t3),dble(t3))
380             continue
390           continue
400         continue
        endif
      end
      subroutine mfftc9(c,fac)
*
*   purpose:
*       elementary cooley-tukey radix 5 step applied to a vector-of
*       2-vectors-of-complex c[ims,nm [ivs,nv [ies,ne]]].
*       see ref.[1] for notations.
*       this routine can be used only by routine mfftis, which controls
*       its operation through the mfftpa common
*
*   dummy arguments :
*
*   c   array being fourier  transformed
*   fac phase factors, prepared by mfftp; not modified in output
*
      implicit double precision (a-h,o-z)
      common /mfftpa/  ims,ivs,ies,nm,nv,ne,mx,lx,mlim,mstep,llim,lstep,
     $ nustep,ivlim,ilim,md2lim,ld2lim
      double complex c(0:nustep-1,0:4),fac(0:*)
      double complex t1,t2,t3,t4,t5
      double precision sin72,rad5d4,s36d72
      parameter (
     $ sin72 =  9.51056516295153572116439333d-1,
     $ rad5d4 =  5.59016994374947424102293417d-1,
     $ s36d72 =  6.18033988749894848204586834d-1 )

      if(lx /= 1)then
          do 200 mu=0,mlim,mstep
            do 150 iv=mu,mu+ivlim,ivs
              ilamf=0
              ilamf2=nustep
              ilamf3=2*nustep
              ilamf4=3*nustep
              do 100 ilam=iv,iv+ilim
                t1=c(ilam,1)*fac(ilamf)+c(ilam,4)*fac(ilamf4)
                t2=c(ilam,2)*fac(ilamf2)+c(ilam,3)*fac(ilamf3)
                t3=(c(ilam,1)*fac(ilamf)-c(ilam,4)*fac(ilamf4))*sin72
                t4=(c(ilam,2)*fac(ilamf2)-c(ilam,3)*fac(ilamf3))*sin72
                t5=t1+t2
                t1=rad5d4*(t1-t2)
                t2=c(ilam,0)-0.25d0*t5
                c(ilam,0)=c(ilam,0)+t5
                t5=t2+t1
                t2=t2-t1
                t1=t3+s36d72*t4
                t3=s36d72*t3-t4
                c(ilam,1)=t5+dcmplx(-dimag(t1),dble(t1))
                c(ilam,4)=t5-dcmplx(-dimag(t1),dble(t1))
                c(ilam,2)=t2+dcmplx(-dimag(t3),dble(t3))
                c(ilam,3)=t2-dcmplx(-dimag(t3),dble(t3))
                ilamf=ilamf+1
                ilamf2=ilamf2+1
                ilamf3=ilamf3+1
                ilamf4=ilamf4+1
100           continue
150         continue
200       continue
      else
          do 400 mu=0,mlim,mstep
            do 350 iv=mu,mu+ivlim,ivs
              do 300 ilam=iv,iv+ilim
                t1=c(ilam,1)+c(ilam,4)
                t2=c(ilam,2)+c(ilam,3)
                t3=(c(ilam,1)-c(ilam,4))*sin72
                t4=(c(ilam,2)-c(ilam,3))*sin72
                t5=t1+t2
                t1=rad5d4*(t1-t2)
                t2=c(ilam,0)-0.25d0*t5
                c(ilam,0)=c(ilam,0)+t5
                t5=t2+t1
                t2=t2-t1
                t1=t3+s36d72*t4
                t3=s36d72*t3-t4
                c(ilam,1)=t5+dcmplx(-dimag(t1),dble(t1))
                c(ilam,4)=t5-dcmplx(-dimag(t1),dble(t1))
                c(ilam,2)=t2+dcmplx(-dimag(t3),dble(t3))
                c(ilam,3)=t2-dcmplx(-dimag(t3),dble(t3))
300           continue
350         continue
400       continue
      endif
      end
      subroutine mfftp1(w,nm,ierr)
*
*   purpose :
*     factorization of nm storing the powers of each factor in w
*     ( nm = 2**w(1)*3**w(2)*... )
*     the maximum factor found  is stored in w(14)
*
      implicit double precision (a-h,o-z)
      parameter (maxfac=3)
      integer w(14),factor(maxfac)
      integer iderr,facerr,tberr
      parameter (iderr=1,facerr=2,tberr=3)
      data (factor(i),i=1,maxfac)/2,3,5/
      n=nm
      do 100 i=1,maxfac
      w(i)=0
 10   if(mod(n,factor(i)) == 0)then
      w(i)=w(i)+1
      n=n/factor(i)
      goto 10
      endif
      if(n == 1)goto 200
100   continue
      ierr=facerr
      return
200   w(14)=i
      end
      subroutine mfftp2(indx,indxi,i2,iw,nm)
*
*     this routine computes two index tables for the
*     permutation due to representation inversion ("bit reversal")
*     the index tables are stored one after the other, and are
*     reciprocal.
*
*     warning: by "bit reversal" we mean the shuffling of indexes
*              as required by in-place fft algorithms, regardless
*              of what is their radix (i.e. 2,3,5).
*              the shuffling is effectively an 'integer bit reversal'
*              only in case of radix-2 algorithms.
*
      implicit double precision (a-h,o-z)
      integer indx(0:*),indxi(0:*),i2(0:*),iw(14)
      parameter (maxfac=3)
      integer factor(maxfac)
      data (factor(i),i=1,maxfac)/2,3,5/
      do 1 i=0,nm-1
        indx(i)=i
1     continue
      ifac=iw(14)
      lx=1
      mx=nm
      do 30 ifac=1,ifac
        do 10 i=1,iw(ifac)-1,2
          mx=mx/factor(ifac)
          call mfftp3(indx,i2,mx,factor(ifac),lx)
          lx=lx*factor(ifac)
          mx=mx/factor(ifac)
          call mfftp3(i2,indx,mx,factor(ifac),lx)
          lx=lx*factor(ifac)
10      continue

*...  if w(ifac)odd,then
        if (i == iw(ifac))then
          mx=mx/factor(ifac)
          call mfftp3(indx,i2,mx,factor(ifac),lx)
          lx=lx*factor(ifac)
          do 20 i=0,nm-1
            indx(i)=i2(i)
20        continue
        endif
30    continue
*...     inverse permutation
         do 40 i=0,nm-1
           i2(i)=indx(i)
           indxi(i)=i
40       continue
         do 59 i=1,nm-3
51         if(i2(i) /= i)then
             idest=i2(i)
             it=indxi(i)
             indxi(i)=indxi(idest)
             indxi(idest)=it
             i2(i)=i2(idest)
             i2(idest)=idest
           goto 51
           endif
59       continue

      end
      subroutine mfftp3(indx,i1,mx,nx,lx)
*
*     this subroutine performs a "bit reversal" permutation
*
      implicit double precision (a-h,o-z)
      integer indx(mx,nx,lx),i1(mx,lx,nx)
      do 1 nu=1,nx
      do 1 mu=1,mx
      do 1 lam=1,lx
        i1(mu,lam,nu)=indx(mu,nu,lam)
1     continue
      end
      subroutine mfftp4(exptab,spetab,factab,n,n1)
*
*     this subroutine builds the twiddle factor tables for use
*     of mfft?s routines (special optimization for small
*     data matrices); it must be used by mfftp only.
*
*     parameters:
*     exptab: twiddle factor table
*     spetab: special twiddle factor table
*     factab: factorization of n
*     n: order of the transform
*     n1: first dimension of the array to be transformed ( >= n)
*
      implicit double precision (a-h,o-z)
      integer maxfac
      parameter(maxfac=3)
      integer factab(14)
      double complex exptab(0:*),spetab(0:*)
      integer factor(maxfac)
      data factor/2,3,5/

      mx=n
      lx=1
      j=0
      do 50 ifact=factab(14),1,-1
         nx=factor(ifact)
         do 40 ipow=1,factab(ifact)
          mx=mx/nx
          do 30 nu=1,nx-1
            do 20 mu=0,mx-1
              do 10 i1=0,n1-1
                spetab(j)=dconjg(exptab(mu*lx*nu))
                j=j+1
10            continue
20          continue
30        continue
          lx=lx*nx
40      continue
50    continue
      mx=n
      lx=1
      j=n*n1
      do 100 ifact=1,factab(14)
        nx=factor(ifact)
        do 90 ipow=1,factab(ifact)
          mx=mx/nx
          do 80 nu=1,nx-1
            do 70 lambda=0,lx-1
              do 60 i1=0,n1-1
                spetab(j)=exptab(mx*lambda*nu)
                j=j+1
60            continue
70          continue
80        continue
          lx=lx*nx
90      continue
100   continue
      end
         subroutine mfftz0(w1,ise1,ne,w2,ise2)
**** purpose : vector copy routine
*    parameters :
*        w1 : vector to be copied
*        ise1: stride of w1
*        ne : number of elements
*        w2 : destination vector
*        ise2: stride of w2
      implicit double precision (a-h,o-z)
         double complex w1(0:*),w2(0:*)
         j=0
         do 1 i=0,(ne-1)*ise1,ise1
         w2(j)=w1(i)
         j=j+ise2
1     continue
         end
