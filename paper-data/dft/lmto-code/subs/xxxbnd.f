      subroutine xxxbnd(evl,nbnd,nqx,nqy,ifib)
      implicit none
      integer nbnd,nqx,nqy,ifib,ix,iy,ib,ifi,fopna
      double precision evl(nbnd,nqy,nqx)
      character*10 strn

      do  ib = 1, nbnd
        ifi = ifib
        if (nbnd /= 1) then
          strn = ' '
          call awrit1('bnd%i',strn,10,0,ib)
          ifi = fopna(strn,30,0)
        endif
        write (ifi,1) nqx,nqy
    1   format('% rows ',i5,' cols ',i5)
        do  ix = 1, nqx
          write (ifi,2) (evl(ib,iy,ix),iy = 1,nqy)
    2     format(6F12.6)
        enddo
      if (nbnd /= 1) call fclose(ifi)
      enddo
      end
      double precision function sclp(v1,v2)
      double precision v1(3),v2(3)
      sclp  = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
      end
