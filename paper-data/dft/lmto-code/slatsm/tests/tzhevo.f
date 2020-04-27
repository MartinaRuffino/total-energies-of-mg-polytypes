      subroutine fmain
      implicit none
      integer ldh,lz,n
      parameter (ldh=3,n=2,lz=2)
      double complex h(ldh,ldh),s(ldh,ldh),evec(lz,lz)
      integer nev,nevl,nevl0
      double precision eps1,eps2,del,emx,evl(n),eo(2),epsovl

      emx = 9d9
      eps1 = 2
      eps2 = 2 + 1d-6
      del =  1d-4 * 1d-3
      epsovl = 1d-8

      h(1,1) = eps1
      h(2,2) = eps2
      h(1,2) = eps1*.999d0
      h(2,1) = eps1*.999d0

      s(1,1) = 1
      s(2,2) = 1
      s(1,2) = 1+del
      s(2,1) = 1+del

C     call zprm('h',12,h,ldh,n,n)
C     call zprm('s',12,s,ldh,n,n)

C     call zhev(n,h,s,.false.,.true.,n,emx,nev,wk,.true.,0,evl,evec)
      nevl = -1
      nevl = 2
      nevl0 = nevl
      eo(1) = 99999d0
      call zhevo(n,ldh,h,s,n,emx,epsovl,0,nevl,nev,evl,eo,n,evec)
      print *, 'nevl(in), nevl(out), nev',nevl0,nevl,nev
      print *, 'eo',eo
      if (nevl0 == 1) call cexit(1,1)
      print *, 'ev',evl
      call prmx('evl',evl,n,nev,1)
      call zprm('z',2,evec,lz,n,nev)
      end

