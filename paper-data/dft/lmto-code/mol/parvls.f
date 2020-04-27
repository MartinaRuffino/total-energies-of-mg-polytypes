      subroutine parvls(sopts,vals,fn)
C- Parse for plotting vals from the command line
C ----------------------------------------------------------------------
Ci Inputs:  sopts (command line string following --plot)
Ci
Co Outputs: vals, fn
Co
Cr Remarks: Adapted from lmf ioden for parsing plot arguments:
Cr          vals(2:4) are three atoms to define the plane of the plot
Cr          vals(5) is the index to the atom at the origin
Cr          vals(1) is the permutation for the axes (123 by default)
Cr          vals(6:7) is the range in x
Cr          vals(8:9) is the range in y
Cr          vals(10) is the number of divisions in x
Cr          vals(11) is the number of divisions in y
Cr          vals(12) is the displacement of the plane in z
Cr          All vals are double, converted to int later.
Cr          fn is the filename ('rho' by default)
Cr Example:
Cr  lmrho --plot:atoms=1,2,3:o=1:fn=plot:x=-3,3:y=-4,4:nx=31:ny=41:dz=0
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      double precision vals(15)
      character*(*) sopts, fn
C Local Variables
      integer int,i,j1,j2,kv(4),ivec(3),length,a2vec,len
      character dc*1

      length = len(sopts)
      dc = sopts(1:1)

      if (dc /= ' ') then
C ... Return here to resume parsing for arguments
        j2 = 0
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) == dc) goto 10
        j1 = min(length,j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
          elseif (sopts(j1:j1+2) == 'fn=')  then
            if (j1+3 > j2) goto 99
            fn = sopts(j1+3:j2)
          elseif (sopts(j1:j1+1) == 'o=')  then
            if (j1+2 > j2) goto 99
            i = j1 + 1
            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,int) /= 1)
     .        goto 99
            vals(5) = int
          elseif (sopts(j1:j1+5) == 'atoms=')  then
            if (j1+6 > j2) goto 99
            i = j1+5
            if (a2vec(sopts,j2,i,2,', '//dc,3,2,3,kv,ivec) /= 3)
     .        goto 99
            vals(2) = ivec(1)
            vals(3) = ivec(2)
            vals(4) = ivec(3)
          elseif (sopts(j1:j1+1) == 'x=')  then
            if (j1+2 > j2) goto 99
            i = j1 + 1
            if (a2vec(sopts,j2,i,4,', '//dc,3,2,3,kv,vals(6)) /= 2)
     .        goto 99
          elseif (sopts(j1:j1+1) == 'y=')  then
            if (j1+2 > j2) goto 99
            i = j1 + 1
            if (a2vec(sopts,j2,i,4,', '//dc,3,2,3,kv,vals(8)) /= 2)
     .        goto 99
          elseif (sopts(j1:j1+2) == 'nx=')  then
            if (j1+3 > j2) goto 99
            i = j1 + 2
            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,int) /= 1)
     .        goto 99
            vals(10) = int
          elseif (sopts(j1:j1+2) == 'ny=')  then
            if (j1+3 > j2) goto 99
            i = j1 + 2
            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,int) /= 1)
     .        goto 99
            vals(11) = int
          elseif (sopts(j1:j1+2) == 'dz=')  then
            if (j1+3 > j2) goto 99
            i = j1 + 2
            if (a2vec(sopts,j2,i,4,' '//dc,2,1,1,kv,vals(12)) /= 1)
     .        goto 99
          elseif (sopts(j1:j1+4) == 'perm=')  then
            if (j1+5 > j2) goto 99
            i = j1 + 4
            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,int) /= 1)
     .        goto 99
            vals(1) = int
          else
            call rxs('ioden: unrecognised option ... ',sopts(j1:j2))
          endif
          goto 10
        endif
      endif
      return
   99 continue
      call rxs('ioden: failed to parse option ... ',sopts(j1:j2))
      end
