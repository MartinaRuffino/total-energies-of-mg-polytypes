      subroutine lgstar(mode,ng,n,gv,ng0,ips0,cg)
C- Compresses F.T. of a real function, using fact it is hermitian
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0, count number of inequivalent points ng0,
Ci             and make ips0.  cg is not used.
Ci         :1  same as mode 0, but also compress cg
Ci         :2  use ips0 to undo compression of cg
Ci   ng    :number of G-vectors
Ci   n     :cg array holds n functions;
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   ips0  :(mode=2) permutation array.
Ci   cg    :list of g-vectors for each of n functions
Co Outputs
Co   ng0   :(mode=0,1) number of inequivalent points
Co   ips0  :(mode=0,1) array of permutation indices.  A negative value
Co         :signifies the point's hermitian counterpart falls earlier in
Co         :the list, and points to that element.
Cr Remarks
Cr   Hacked from svgsym, using only one symmetry operation
Cr   to reduce (G,-G) pairs to a single element.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer mode,ng,n,ips0(ng),ng0,iprint
      double precision gv(ng,3)
      double complex cg(ng,n)
      integer i,i0,i00,irep,k,j,j0,lwarn,m
      double precision v(3),df

      if (mode < 0 .or. mode > 2) call rxi('lgstar, bad mode',mode)
      lwarn = 0
      if (mode == 2) goto 200

C --- mode = 0,1 ---
      ng0 = 0
      do  i = 1, ng
        ips0(i) = 0
      enddo

C --- Main loop: look for next unclassified vector ---
      i00 = 1
      do  irep = 1, ng+1
        i0 = 0
        do  i = i00, ng
          i0 = i
          if (ips0(i) == 0) goto 80
        enddo
        goto 81
   80   continue

C   ... Apply all point ops, find in list, add to phase sum
        ng0 = ng0 + 1
        ips0(i0) = ng0
        if (mode == 1)  then
          do  m = 1, n
            cg(ng0,m) = cg(i0,m)
          enddo
        endif
        do  k = 1, 1
          v(1) = gv(i0,1)
          v(2) = gv(i0,2)
          v(3) = gv(i0,3)
          do  j = i0+1, ng
            df = (v(1)+gv(j,1))**2+(v(2)+gv(j,2))**2+(v(3)+gv(j,3))**2
            j0 = j
            if (df < 1d-8) goto 70
          enddo
C     ... No matching vector here ... should only happen for G=0
          i00 = i0
          goto 10
   70     continue
          ips0(j0) = -i0
          if (mode == 1) then
            if (abs(cg(i0,1)-dconjg(cg(j0,1))) > 1d-9) lwarn = lwarn+1
          endif
        enddo
        i00 = i0
   10   continue
      enddo
      call rxi('bug in lgstar, irep=',irep)
   81 continue
      if (lwarn > 1 .and. iprint() >= 10) print 1,lwarn
    1 format(' lgstar (warning):',i6,' points not hermitian')
      return

C --- mode = 2 ---
  200 continue

C ... Unpack original points first
      do  m = 1, n
        do  i = ng, 1,-1
          k = ips0(i)
          if (k > 0) cg(i,m) = cg(k,m)
        enddo
      enddo

C ... Unpack hermitian points
      do  m = 1, n
        do  i = 1, ng
          k = -ips0(i)
          if (k > 0) cg(i,m) = dconjg(cg(k,m))
        enddo
      enddo
      end
