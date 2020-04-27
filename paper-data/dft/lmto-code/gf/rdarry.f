      integer function rdarry(l1,l2,s,array,n,ldble,errh,ifi)
C- File read array with option to skip over record or match to array
C ----------------------------------------------------------------------
Ci Inputs
Ci   l1    :0 do not read array
Ci         :1 read next record into 'array'
Ci         :2 read next record into temporaray space and require that it
Ci            match 'array'
Ci   l2    :F record is missing from file
Ci         :T record is present in file
Ci   s     :string describing record (for error messages)
Ci   n     :size of array
Ci   ldble :F arrray is integer
Ci         :T array is double
Ci   errh  :abort with error message if
Ci   ifi   :file logical unit
Co Outputs
Ci   array :array to be (optionally) read from file unit ifi
Cr Remarks
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical l2,ldble,lddump,lidump
      character *(*) s
      integer l1,n,ifi,errh
C ... Dynamically allocated local arrays
      integer, allocatable :: iwk(:)
      real(8), allocatable :: wk(:)
C ... Local parameters
      double precision array(n)
      character ss*20
      integer k,idamax,iiamax,ix
      double precision xx

      rdarry = 0
      if (l1 /= 0 .and. .not. l2) goto 99
      if (.not. l2) return

C ... Array is present but is skipped over
      if (l1 == 0) then
        read(ifi)

C ... Read array into 'array'
      elseif (l1 == 1) then

        if (ldble) then
          if (.not. lddump(array,n,ifi)) goto 99
        else
          if (.not. lidump(array,n,ifi)) goto 99
        endif
      elseif (l1 == 2) then
        if (ldble) then
          allocate(wk(n))
          if (.not. lddump(wk,n,ifi)) goto 99
          call daxpy(n,-1d0,array,1,wk,1)
          k = idamax(n,wk,1)
          xx = wk(k)
          deallocate(wk)
          if (xx /= 0) goto 99
        else
          allocate(iwk(n))
          if (.not. lidump(iwk,n,ifi)) goto 99
          call iaxpy(n,-1,array,1,iwk,1)
          k = iiamax(n,iwk,1)
          ix = iwk(k)
          deallocate(iwk)
          if (ix /= 0) goto 99
        endif
      else
        call sanrg(.true.,l1,0,2,'rdarry:','l1')
      endif
C     Read array into work area and match to array

      return

C --- Error handling when failed to read or match array ---
   99 continue
      if (errh == 0) call rxs('RDARRY: file mismatch, array ',s)
      if (errh > 0) then
        ss = s
        print *, 'RDARRY (warning) failed to read array '//ss
      endif
      rdarry = -1

      end
