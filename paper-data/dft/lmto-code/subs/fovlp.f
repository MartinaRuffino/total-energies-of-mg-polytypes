      subroutine fovlp(ib1,ib2,ntab,iax,plat,pos,ipc,alat,rmax,z,pwr,
     .  facaa,facae,facee,fmax,f,inc)
C- Analytical function of the sphere overlaps
Cu Updates
Cu   08 Nov 17  calculates ES max overlaps separately
Cu   22 Oct 02  weight ES-ES and atom-ES overlaps differently
Cu              New argument list
      implicit none
      integer ib1,ib2,ntab(*),niax,ipc(*),inc
      parameter (niax=10)
      integer iax(niax,*)
      double precision alat,plat(3,3),pos(3,*),rmax(*),pwr,fmax(3),f,
     .  facaa,facae,facee
      double precision dsqr,dr,d,sumrs,z(*),zb,zi,fac
      integer i,ib,i0,i1,ii,jj,kk,ix,icb,ici
      integer, parameter :: NULLI=-99999

      fmax = NULLI
      f = 0
      inc = 0
      do  ib = ib1, ib2
        i0 = ntab(ib)+1
        i1 = ntab(ib+1)
        do  i = i0+1, i1
          ii = iax(3,i)-iax(3,i0)
          jj = iax(4,i)-iax(4,i0)
          kk = iax(5,i)-iax(5,i0)
          dsqr = 0
          do  ix = 1, 3
            dr = pos(ix,iax(2,i))-pos(ix,iax(1,i0))+
     .        plat(ix,1)*ii+plat(ix,2)*jj+plat(ix,3)*kk
            dsqr = dsqr + dr**2
          enddo

          icb = ipc(ib)
          ici = ipc(iax(2,i))
          zb = z(icb)
          zi = z(ici)
          fac = 1
          if (zb /= 0 .and. zi /= 0) then
            fac = facaa
          elseif (zb == 0 .and. zi == 0) then
            fac = facee
          else
            fac = facae
          endif
          sumrs = rmax(icb) + rmax(ici)
          d = alat*dsqrt(dsqr)
          if (sumrs > d) then
            f = f + fac*(sumrs/d-1)**pwr
            inc = inc+1
C          else
C            goto 10
          endif
          if (zb /= 0 .and. zi /= 0) then
            fmax(1) = max(sumrs/d-1,fmax(1))
          elseif (zb == 0 .and. zi == 0) then
            fmax(3) = max(sumrs/d-1,fmax(3))
          else
            fmax(2) = max(sumrs/d-1,fmax(2))
          endif
        enddo
      enddo
      end
