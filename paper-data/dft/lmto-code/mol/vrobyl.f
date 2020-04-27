      subroutine vrobyl(nr1,nlml,vl,rl,rofi,rwgt,z,sum,tot,job)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sum(0:10),rl(nr1,1),vl(nr1,1),rofi(1),rwgt(1)
      lmaxl=ll(nlml)
      do 33 l=0,lmaxl
  33  sum(l)=0d0
      if(job == 1) then
        xx=2d0*z*dsqrt(16d0*datan(1d0))
        do 3 ir=2,nr1
  3     sum(0)=sum(0)+rl(ir,1)*(vl(ir,1)-xx/rofi(ir))*rwgt(ir)
      else
        do 4 ir=2,nr1
  4     sum(0)=sum(0)+rl(ir,1)*vl(ir,1)*rwgt(ir)
      endif
      do 30 ilm=2,nlml
      l=ll(ilm)
      do 30 ir=2,nr1
  30  sum(l)=sum(l)+ rl(ir,ilm)*vl(ir,ilm) * rwgt(ir)
      tot=0d0
      do 35 l=0,lmaxl
  35  tot=tot+sum(l)
      end
