      subroutine dpdump(array,length,ifile)
C- Binary I/O of an array
      double precision array(length)
      if (ifile > 0) read(ifile) array
      if (ifile < 0) write(-ifile) array
      end
      subroutine dpsdmp(array,n1,n2,ifile)
C- Binary I/O of an array segment
      integer n1,n2,ifile,length
      double precision array(n2)
      length = n2-n1+1
      if (length > 0) call dpdump(array(n1),length,ifile)
      end
