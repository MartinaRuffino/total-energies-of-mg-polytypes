      subroutine shoblk(nrow,ncol,ndim,nsr,nsc,iax,s,iopt)
C- Prints a nl x nl subblock of a two-dimensional matrix to stdout
      implicit none
      integer nrow,ncol,ndim,nsr,nsc,iopt,niax
      double precision s(0:ndim-1,0:*)
      parameter (niax=10)
      integer iax(niax,*)
      integer n1,n2,i,j,l1,l2

   10 continue
      print *, nrow,ncol,' rows and columns of blocks. block i,j? '
      read (*,*) n1,n2
      if (n1+n2 == 0) return

      i = (n2-1)*nrow + n1
      write (*,1) i,(iax(j,i),j=1,5)
    1 format(' ','vector ',i3,' atom1 ',i3,' atom2 ',i3,' transl ',3I5)

      n1 = (n1-1)*nsr
      n2 = (n2-1)*nsc
      do  l1 = n1, n1+nsr-1
        if (iopt == 0) then
          write (*,2) (s(l1,l2),l2=n2,n2+nsc-1)
        else
          write (*,2) (1000*s(l1,l2),l2=n2,n2+nsc-1)
        endif
      enddo
      goto 10
    2 format(1x,9F8.4)
      end
