      subroutine cinit(array,leng)
C- Initializes complex array to zero
      implicit none
      integer leng
      real array(2*leng)
      call ainit(array,leng+leng)
      end
      subroutine zinit(array,leng)
C- Initializes complex array to zero
      implicit none
      integer leng
      double precision array(2*leng)
      call dpzero(array,leng+leng)
      end
      subroutine iinit(array,leng)
C- Initializes integer array to zero
      implicit none
      integer leng,i
      integer array(leng)
      forall (i=1:leng) array(i) = 0
      end
      subroutine ainit(array,leng)
C- Initializes real array to zero
      implicit none
      integer leng,i
      real array(leng)
      forall (i=1:leng) array(i) = 0
      end
      subroutine dpzero(array,leng)
C- Initializes double precision array to zero
      implicit none
      integer leng,i
      double precision array(leng)
      forall (i=1:leng) array(i) = 0
      end
      real function rval(array,index)
C- Returns the real value of ARRAY(INDEX)
      implicit none
      integer index
      real array(index)
      rval = array(index)
      end
      double precision function dval(array,index)
C- Returns the double precision value of ARRAY(INDEX)
      implicit none
      integer index
      double precision array(index)
      dval = array(index)
      end
      integer function ival(array,index)
C- Returns the integer value of ARRAY(INDEX)
      implicit none
      integer index
      integer array(index)
      ival = array(index)
      end
      integer function ival2(array,nda,i1,i2)
C- Returns the value of integer array(i1,i2)
      implicit none
      integer nda,i1,i2,array(nda,1)
      ival2 = array(i1,i2)
      end
      logical function logval(array,index)
C- Returns the integer value of ARRAY(INDEX)
      implicit none
      integer index
      logical array(index)
      logval = array(index)
      end
      complex function cval(array,index)
C- Returns the complex value of ARRAY(INDEX)
      implicit none
      integer index
      complex array(index)
      cval = array(index)
      end
      subroutine dvset(array,i1,i2,val)
C- Sets some elements of double precision array to value
      implicit none
      integer i1,i2
      double precision array(i2),val
      integer i
      forall (i=i1:i2) array(i) = val
      end
      subroutine dvadd(array,i1,i2,val)
C- Adds value some elements of double precision array
      implicit none
      integer i1,i2
      double precision array(i2),val
      integer i
      do  10  i = i1, i2
   10 array(i) = array(i) + val
      end
      subroutine ivset(array,i1,i2,val)
C- Sets some elements of integer array to value
      implicit none
      integer i1,i2,array(i2),val,i
      forall (i=i1:i2) array(i) = val
      end
      subroutine lvset(array,i1,i2,val)
C- Sets some elements of logical array to value
      implicit none
      integer i1,i2,i
      logical array(i2),val
      forall (i=i1:i2) array(i) = val
      end
