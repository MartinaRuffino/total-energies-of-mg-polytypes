      subroutine swapdV(nl2,nsites,nsp,iax,dh)
C- Cludge to swap Vsp with Vps and Vpd with Vdp and Vsd with Vds
C ----------------------------------------------------------------
C Remarks
C  Force dh(r1,l1,r2,l2) = dh(r2,l2,r1,l1)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nl2,nsites,nsp,niax
      parameter (niax=10)
      integer iax(niax,nsites)
      double precision dh(nl2,nl2,nsites*nsp**2,3)
C Local parameters
      integer i,j,k,lm1,lm2,nsp2,isp,ii,jj
      double precision temp
C     .  ,d1mach

      if (nl2 < 2) return

      nsp2 = nsp**2
      do  isp = 1, nsp2
      do  i = 1, nsites
        ii = i+(isp-1)*nsites
        do  j = i, nsites
          jj = j+(isp-1)*nsites
          if (iax(1,j) == iax(2,i) .and. iax(2,j) == iax(1,i) .and.
     .      iax(3,i) == -iax(3,j) .and.
     .      iax(4,i) == -iax(4,j) .and.
     .      iax(5,i) == -iax(5,j)) then

C --- swap sp,ps ---
              do  lm1 = 2, 4
                do  k = 1, 3
C                print *, ii,jj,k,lm1
C                print 333, dh(lm1,1,ii,k), dh(lm1,1,jj,k)
C  333           format(2g17.6)
C                if (dh(lm1,1,ii,k)*dh(lm1,1,jj,k) < -d1mach(3)) then
C                  print*,' i = ',ii,'  j = ',jj,'  k = ',k,
C     .              '  lm1 = ',lm1
C                  call rx('bug in swapdV')
C                endif
                temp = dh(lm1,1,ii,k)
                dh(lm1,1,ii,k) = dh(lm1,1,jj,k)
                dh(lm1,1,jj,k) = temp
                enddo
              enddo
              if (nl2 >= 5) then

C --- swap sd,ds ---
              do  lm1 = 5, 9
              do  k = 1, 3
C                print *, ii,jj,k,lm1
C                print 333, dh(lm1,1,ii,k), dh(lm1,1,jj,k)
C                if (dh(lm1,1,ii,k)*dh(lm1,1,jj,k) > d1mach(3)) then
C                  print*,' i = ',ii,'  j = ',jj,'  lm1 = ',lm1
C                  call rx('bug in swapdV')
C                endif
                temp = dh(lm1,1,ii,k)
                dh(lm1,1,ii,k) = -dh(lm1,1,jj,k)
                dh(lm1,1,jj,k) = -temp
              enddo
              enddo

C --- swap pd,dp ---
              do  lm1 = 5, 9
              do  lm2 = 2, 4
              do  k = 1, 3
C                  print *, ii,jj,k,lm1,lm2
C                  print 333, dh(lm1,lm2,ii,k), dh(lm1,lm2,jj,k)
C                  if (dh(lm1,lm2,ii,k)*dh(lm1,lm2,jj,k) < -d1mach(3))
C     .              then
C                    print*,' i = ',ii,'  j = ',jj,'  lm1 = ',lm1,
C     .                '  lm2 = ',lm2
C                    call rx('bug in swapdV')
C                  endif
                temp = dh(lm1,lm2,ii,k)
                dh(lm1,lm2,ii,k) = dh(lm1,lm2,jj,k)
                dh(lm1,lm2,jj,k) = temp
              enddo
              enddo
              enddo
              endif
              exit
            endif
          enddo
        enddo
      enddo
      end
