      subroutine iorodc(nri,rhoi,nr,lmxl,a,orho,nbas,ips,ifi)
c  input/output for decomposed density
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhoi(nri),ips(1),lmxl(1),nr(1),orho(1),a(1)
      real w(1)
      common /w/ w
      jfi=iabs(ifi)
      rewind jfi
c ------- write branch --------------
      if(ifi < 0) then
        write(jfi) nri
        call dpdump(rhoi,nri,-jfi)
        do 10 ib=1,nbas
        is=ips(ib)
        nr1=nr(is)
        nlml=(lmxl(is)+1)**2
        write(jfi) nr1,nlml,a(is)
  10    call dpdump(w(orho(ib)),nr1*nlml,-jfi)
      endif
c ------- read branch: also defines arrays --------------
      if(ifi > 0) then
        read(jfi) nridum
        if(nri /= nridum) call rx('iorodc-read:  nri disagrees')
        call dpdump(rhoi,nri,jfi)
        do 11 ib=1,nbas
        is=ips(ib)
        read(jfi) nr(is),nlmdum,a(is)
        nr1=nr(is)
        nlml=(lmxl(is)+1)**2
        if(nlml /= nlmdum) call rx('iorodc-read:  nlml disagrees')
      write(6,945) ib,nr1,nlml,a(is)
  945 format(' ib=',i5,'   nr=',i5,'   nlml=',i5,'   a=',f9.3)
        call defrr(orho(ib),   nr1*nlml)
  11    call dpdump(w(orho(ib)),nr1*nlml,jfi)
      endif
      end
