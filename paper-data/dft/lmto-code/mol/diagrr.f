      subroutine diagrr(ndim,n,smat,hmat,tln,evl)
      parameter( nmax=800 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension smat(ndim,n),hmat(ndim,n),tln(ndim,n),evl(n)
      dimension ee(nmax),dl(nmax),ibl(nmax),iperm(nmax),nbl0(nmax+1)
      if(n > nmax) call rx('change dimension nmax in diagrr')
      call reduc1(ndim,n,hmat,smat,dl)
      call tred2(ndim,n,hmat,evl,ee,tln)
      call imtql2(ndim,n,evl,ee,tln,ierr)
      call rebaka(ndim,n,1,n,smat,dl,tln)
      if(iprint() >= 30) write(6,600) evl
  600 format(' evl='/(1x,8f10.5))
C|99  ihs=0
C|    write(6,*) 'enter ihs'
C|    read (5,*)  ihs
C|    if(ihs <= 0) return
C|    write(6,677) evl(ihs),(tln(i,ihs),i=1,n)
C|677 format(' e=',f12.5/(1x,9f8.4:))
C|    if(.true.) goto 99
      return
      end
