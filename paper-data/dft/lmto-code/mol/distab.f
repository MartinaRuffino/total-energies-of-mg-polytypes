      subroutine distab(nbas,ips,spid,alat,pos,ifi)
C- Print distance table
      implicit real*8 (a-h,p-z), integer (o)
      dimension ips(1),pos(3,nbas),d(50)
      character*8 spid(1)
      n1=8
      nrep=(nbas-1)/n1+1
      do 10 irep=1,nrep
      i1=(irep-1)*n1+1
      i2=min0(irep*n1,nbas)
      do 10 jrep=irep,nrep
      j1=(jrep-1)*n1+1
      j2=min0(jrep*n1,nbas)
      if(ifi /= 71) write(6,*) ' '
      write(ifi,120) (spid(ips(j)),j=j1,j2)
  120 format(9x,8(4x,a4))
      write(ifi,121) (j,j=j1,j2)
  121 format(6x,8i8)
      if(ifi /= 71) write(ifi,122) ('--------',j=j1,j2)
  122 format(7x,'---',8a8)
      do 11 i=i1,i2
      do 12 j=j1,j2
  12  call dpdist(pos(1,i),pos(1,j),3,d(j))
      write(ifi,130) spid(ips(i)),i,(alat*d(j),j=j1,j2)
  130 format(1x,a4,i4,8f8.3)
  11  continue
  10  continue
      end
