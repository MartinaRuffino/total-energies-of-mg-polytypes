      subroutine evlwgt(qval,nhs,evl,ewt,nstate)
C- Set weights for eigenvalues, detect degeneracies
      implicit real*8 (a-h,p-z), integer (o)
      dimension evl(1),ewt(1)
      tol=1d-5
      nqval=qval+0.01
      nstat0=(qval+0.01)/2
      if(2*nstat0 /= nqval) then
C|       write(6,*) 'odd electron num -- half occupation of top level'
         nstat0=nstat0+1
      endif
c ---- look for degeneracy at highest occ state --------
      do 11 i=nstat0,1,-1
  11  if(dabs(evl(i)-evl(nstat0)) < tol) is1=i
      do 12 i=nstat0,nhs
  12  if(dabs(evl(i)-evl(nstat0)) < tol) is2=i
      nd0=is2-is1+1
c  next two lines force normal occupation
C|      is1=nstat0
C|      is2=nstat0
      ndeg=is2-is1+1
      qrest=nqval-2*(is1-1)
      wfrac=0.5d0*qrest/ndeg
      nstate=is2
      do 1 is=1,nstate
  1   ewt(is)=1d0
      do 2 is=is1,is2
  2   ewt(is)=wfrac
      if(iprint() >= 20) write(6,300) nqval,nd0,is1,is2,evl(nstate)
  300 format(/' nqval=',i4,'   degen=',i3,'   (states',i4,
     . '  to',i4,')    evl=',f11.5)
      write(71,710) nqval,is1,is2,nd0,evl(nstate),evl(nstate+1)
  710 format(' mce  qv',i5,'   is1,is2',2i5,'   deg',i2,'  ev',f9.5,
     .   '  nx',f9.5)
      end
