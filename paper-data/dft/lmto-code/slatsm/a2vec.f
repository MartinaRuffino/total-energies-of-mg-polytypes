      integer function a2vec(str,lstr,ip,cast,sep,nsep,itrm,nvec,ix,res)
C- Parses ascii string for a vector of binary values
C ----------------------------------------------------------------------
C Inputs:
C   str(ip:lstr-1):string to parse, from (ip=0 for first char)
C   cast:        0=logical, 2=int, 3=real, 4=double
C   sep(1:nsep): class of chars that separate arguments
C   itrm,nvec:   index to sep which indicates last arg: parse terminates
C                when sep(i>=itrm) encountered or nvec arguments parsed.
C                nvec<0: sign used to allow nam=expr as an expression
C                which has the additional effect of loading the nam
C                with expr into the variables table.
C                itrm<0:  sign used to skip over white space before
C                parsing for next token.
C Outputs:
C   res:         binary result
C   ip:          position in str on exit
C   ix:          indices to which characters terminated expression
C   a2vec:       number of values actually converted.
C                If not all values are parsed, a2vec returns -1-#
C                where # is the number of arguments successfully parsed
C Remarks
Cu Updates
Cu   04 Aug 13   Small changes (D. Pashov)
Cu   25 Jul 11   Bug fix case no terminator found
Cu   26 May 07   When last expr ends at lstr w/out term, return w/out
Cu               error only if a2bin successfully parses expr
Cu   02 Feb 01   return w/out error when last expr ends at lstr w/out term
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ip,lstr,cast,nsep,nvec,itrm,ix(*)
      character str*(*), sep(*)
      double precision res(*)
C ... Local parameters
      integer iv,ip0,iitrm
      logical a2bin
      character term

      iitrm = iabs(itrm)

C --- For each element do ---
      a2vec = 0
      do  10  iv = 1, iabs(nvec)

C   --- Determine appropriate terminator, possibly exit ---
C   ... ip0 points to start of next token
   12   if (itrm < 0) call skipbl(str,lstr,ip)
        if (ip >= lstr) return
        ip0 = ip
        call chrps2(str,sep,nsep,lstr,ip0,ix(iv))
        if (ip0 >= lstr) ix(iv) = 0
C       No terminator was found: possibly exit with error
        if (ix(iv) == 0) then
          term = sep(1)
C         Still ok if expr runs to last char in string with no term
          ip0 = ip
          if (.not.a2bin(str,res,cast,iv-1,term,ip0,lstr-1)) then
            a2vec = -iv
            return
          endif
        else
          term = sep(ix(iv))
        endif
        if (.not. a2bin(str,res,cast,iv-1,term,ip,lstr-1)) then
C ...     Try for nam=expr
          if (nvec < 0) then
            call skipbl(str,lstr,ip)
            ip0 = ip
            call parsyv(str,lstr,1,0,ip)
            if (ip == ip0) then
              a2vec = -iv
              return
            endif
C           call shosyv(0,0,0,6)
            if (ix(iv) >= iitrm) return
            goto 12
          endif
C ...     Error exit
          a2vec = -iv
          return
        endif
        a2vec = iv
        if (ix(iv) >= iitrm) return
   10 continue

      end
