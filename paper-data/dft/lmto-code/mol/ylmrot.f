      subroutine ylmrot(l,alfa,beta,gama,a,ndim)
c  makes (nm x nm) transformation matrix describing how the spherical
c  harmonics of one l mix under rotation (cf bradley/cracknell ch 2).
c  def:  ylm(m1,rot(p))= sum a(m1,m2)*ylm(m2,p)
      parameter( lmx=20)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(ndim,ndim),fac(0:2*lmx)
      complex*16 cpp,cpm,cmp,cmm,ci,cexpm,cexpn
      isig(nn)=1-2*mod(iabs(nn),2)
      if(l > lmx) stop '*** l gt lmx in ylmrot'
      nm=2*l+1
      if(nm > ndim) stop '*** nm gt ndim in ylmrot'
      ci=(0.d0,1.d0)
      sr2=dsqrt(2.d0)
      fac(0)=1.d0
      do 1 i=1,2*l
  1   fac(i)=fac(i-1)*i
      lp1=l+1
c ------ make the matrix dnm(beta) --------
      sbet=dsin(.5d0*beta)
      cbet=dcos(.5d0*beta)
      if(dabs(sbet) < 1.d-8) then
        do 5 m=1,nm
        do 6 n=1,nm
  6     a(m,n)=0.d0
  5     a(m,m)=1.d0
      else
      ctbet=cbet/sbet
      do 2 m=-l,l
      a(m+lp1,l+lp1)=isig(l-m)*dsqrt(fac(l+l)/(fac(l+m)*fac(l-m)))
     .    *cbet**(l+m)*sbet**(l-m)
  2   continue
      do 3 n=l-1,0,-1
      do 3 m=-n,n
      rhs=(m-n-1)*ctbet*a(m+lp1,n+1+lp1)
     .    -dsqrt((l+m)*(l-m+1.d0))*a(m-1+lp1,n+1+lp1)
  3   a(m+lp1,n+lp1)=rhs/dsqrt((l+n+1.d0)*(l-n))
      do 4 n=0,l
      do 4 m=-n,n
      a(-n+lp1,-m+lp1)=a(m+lp1,n+lp1)
      a(n+lp1,m+lp1)=isig(m+n)*a(m+lp1,n+lp1)
  4   a(-m+lp1,-n+lp1)=a(n+lp1,m+lp1)
      endif
c ------- make into real matrix -----
      do 40 m=0,l
c>    cexpm=cdexp(dcmplx(0d0,-m*gama))
      cexpm=cdexp(dcmplx(0d0,-m*gama))
      do 40 n=0,l
c>    cexpn=cdexp(dcmplx(0d0,-n*alfa))
      cexpn=cdexp(dcmplx(0d0,-n*alfa))
      xxx=0.5d0
      if(m == 0) xxx=xxx/sr2
      if(n == 0) xxx=xxx/sr2
      cpp=a(m+lp1,n+lp1)*cexpm*cexpn
c>    cmp=a(-m+lp1,n+lp1)*dconjg(cexpm)*cexpn*isig(m)
      cmp=a(-m+lp1,n+lp1)*dconjg(cexpm)*cexpn*isig(m)
c>    cpm=a(m+lp1,-n+lp1)*cexpm*dconjg(cexpn)*isig(n)
      cpm=a(m+lp1,-n+lp1)*cexpm*dconjg(cexpn)*isig(n)
c>    cmm=a(-m+lp1,-n+lp1)*dconjg(cexpm)*dconjg(cexpn)*isig(m)*isig(n)
      cmm=a(-m+lp1,-n+lp1)*dconjg(cexpm)*dconjg(cexpn)*isig(m)*isig(n)
      a(m+lp1,n+lp1)=(cpp+cmp+cpm+cmm)*xxx
      if(m > 0) a(-m+lp1,n+lp1)=ci*(cpp-cmp+cpm-cmm)*xxx
      if(n > 0) a(m+lp1,-n+lp1)=-ci*(cpp+cmp-cpm-cmm)*xxx
  40  if(m > 0.and.n > 0) a(-m+lp1,-n+lp1)=(cpp-cmp-cpm+cmm)*xxx
c -------------------------------------
c|    do 46 m=1,nm
c|46  write(6,360) (a(m,n),n=1,nm)
c|360 format(1x,9f8.5)
      return
      end
