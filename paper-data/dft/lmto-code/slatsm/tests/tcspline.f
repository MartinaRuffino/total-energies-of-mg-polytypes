      subroutine fmain
C     driver for routine splint, which calls spline
      implicit none
      INTEGER NP
      double precision PI,x0
      PARAMETER(NP=10)
      INTEGER i,nfunc
      double precision x,yp1,ypn,xa(NP),ya(NP),y2(NP)
      double precision f,fp,fi,y,yp,yi

      pi = 4*datan(1d0)
      do 14 nfunc=1,2
        if (nfunc == 1) then
          write(*,*) 'Interpolate sine function from 0 to pi'
          do 11 i=1,np
            xa(i)=i*pi/np
            ya(i)=dsin(xa(i))
   11     continue
          yp1=dcos(xa(1))
          ypn=dcos(xa(np))
        else if (nfunc == 2) then
          write(*,*) 'Interpolate exponential function from 0 to 1'
          do 12 i=1,np
            xa(i)=1.0d0*i/np
            ya(i)=dexp(xa(i))
   12     continue
          yp1=dexp(xa(1))
          ypn=dexp(xa(np))
        else
          stop
        endif
        call cspline(xa,ya,np,yp1,ypn,y2)
        write(*,'(1x,t10,a1,t20,a4,t28,a13,t43,a,8x,a,3x,a,4x,a)')
     .    'x','f(x)','interpolation','diff',
     .    'f''(x)','interpolation','diff'
        do 13 i=1,11
          if (nfunc == 1) then
            x=(-0.05d0+i/10.0d0)*pi
            f=dsin(x)
            fp=dcos(x)
          else if (nfunc == 2) then
            x=-0.05d0+i/10.0d0
            f=dexp(x)
            fp = f
          endif
          call csplint(xa,ya,y2,np,x,y,yp)
          write(*,'(1x,3f12.6,1pe12.2,0p,2f12.6,1pe12.2)')
     .      x,f,y,f-y,fp,yp,fp-yp
   13   continue
        write(*,*) '-----------------------------------'
        if (nfunc == 1) then
          write(*,*) 'Integrate sine function from 1 to x'
          x0 = 1
C          x0 = pi-.01*1
        else if (nfunc == 2) then
          x0 = 0
          write(*,*) 'Integrate exponential function from 0 to x'
        endif
        write(*,'(1x,t10,a1,t20,a4,t28,a13,t43,a,8x,a,3x,a,4x,a)')
     .    'x','f(x)','interpolation','diff',
     .    'int f','interpolation','diff'
       do  i = 1, 11
          if (nfunc == 1) then
            x=(-0.05d0+i/10.0d0)*pi
            f=dsin(x)
            fp=dcos(x)
            fi = -dcos(x) + dcos(x0)
          else if (nfunc == 2) then
            x=-0.05d0+i/10.0d0
            f=dexp(x)
            fp = f
            fi = f - dexp(x0)
          endif
          call csplint(xa,ya,y2,np,x,y,yp)
          call csintg(xa,ya,y2,np,x0,x,yi)
          write(*,'(1x,3f12.6,1pe12.2,0p,2f12.6,1pe12.2)')
     .      x,f,y,f-y,fi,yi,fi-yi
        enddo

C       write(*,*) 'Press RETURN'
C       read(*,*)
   14 continue





      end
C      subroutine rx(string)
C      character *(*) string
C
C      print *, string
C      stop
C      end
