      subroutine prerot(bin,bout,nlm,rmat,nrot)
c  pre-rotates coefficient vector for 2-center expansion
      implicit real*8 (a-h,o-z)
      dimension bin(nlm),bout(nlm),rmat(nrot,nrot)
      do 12 ilm=1,nlm
  12  bout(ilm)=0d0
      do 10 jlm=1,nlm
      do 10 ilm=1,nlm
  10  bout(ilm)=bout(ilm)+rmat(ilm,jlm)*bin(jlm)
      end
      subroutine posrot(bin,bout,nlm,rmat,nrot)
c  post-rotates output vector made with 2-center expansion
      implicit real*8 (a-h,o-z)
      dimension bin(nlm),bout(nlm),rmat(nrot,nrot)
C     do 12 ilm=1,nlm
C 12  bout(ilm)=0d0
      do 10 ilm=1,nlm
      bout(ilm)=0d0
      do 10 jlm=1,nlm
  10  bout(ilm)=bout(ilm)+rmat(jlm,ilm)*bin(jlm)
      end
