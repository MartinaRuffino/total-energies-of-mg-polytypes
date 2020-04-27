      subroutine dfftqp(plat,jr,pos,n1,n2,n3,qp,nq,jq)
C- Fourier transform of real function at a specified list of q-points
      implicit none
      integer n1,n2,n3,nq
      double precision jr(n1,n2,n3),plat(3,3),pos(3),qp(3,nq)
      double complex jq(nq)
C Local variables
      integer iq,j1,j2,j3
      double precision twopi,sp,p(3)
      double complex phase

      twopi = 8d0*datan(1d0)
      call dpzero(jq,2*nq)
      do  iq = 1, nq

      do  j1 = 0, n1-1
      do  j2 = 0, n2-1
      do  j3 = 0, n3-1
C   ... multiply pos*0 for FFT convention
C   ... multiply pos*1 for real jq
        p(1) = j1*plat(1,1)+j2*plat(1,2)+j3*plat(1,3) ! - pos(1)
        p(2) = j1*plat(2,1)+j2*plat(2,2)+j3*plat(2,3) ! - pos(2)
        p(3) = j1*plat(3,1)+j2*plat(3,2)+j3*plat(3,3) ! - pos(3)
        sp = twopi*(p(1)*qp(1,iq) + p(2)*qp(2,iq) + p(3)*qp(3,iq))
        phase = dcmplx(dcos(sp),dsin(sp))
        jq(iq) = jq(iq) + phase*jr(j1+1,j2+1,j3+1)
      enddo
      enddo
      enddo
      enddo

      end
