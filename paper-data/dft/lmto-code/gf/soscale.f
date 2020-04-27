      subroutine soscale(mode,ldi,ldj,pfdim,offpi,offpj,P,g)
C- Scale 2x2 matrix g <- transpose(P) * g * P in SO case
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :not used now
Ci         :10s digit
Ci         :0 SR convention for m ordering
Ci         :1 FR convention for m ordering
Ci   ldi   :leading (row) dimension of g (may be a subblock of the entire g)
Ci   ldj   :second (column) dimension of g (may be a subblock of the entire g)
Ci   pfdim :leading dimension of P
Ci   offpi :offset in P to first row in g subblock
Ci         :P(offpi+1) points to the same row in the full hamiltonian  as g(1,:)
Ci   offpj :offset in P to first column in g subblock
Ci         :P(offpj+1) points to the same row in the full hamiltonian  as g(1,:)
Ci   P     :some potential function such as P^alp/Pgam or sqrt(Pdot) (mksopf.f)
Ci         :with the properties described in Remarks
Cio Inputs/Outputs
Cio  g     :Green's fumction (possibly diagonal subblock) to scale (overwritten)
Cio        :g(1,1) is the first element in the subblock
Cr Remarks
Cr   This routine does the matrix product
Cr     X^ss',ij = sum_kk',lm  (Ptrans)^sk_il g^kk'_lm P^k's'_mj
Cr     s,s' and k,k' are spin indices ranging from 1 to 2, ijlm are orbital indices
Cr   with P having selection rules (SO case)
Cr     P^ss_mj ~ d_mj    (d_il is Kroneker delta: P^ss is a diagonal matrix)
Cr     P^12_mj ~ d_m-1,j
Cr     P^21_mj ~ d_m+1,j
Cr   and P having selection rules (FR case)
Cr     P^ss_mj ~ d_mj
Cr     P^12_mj ~ d_m+1,j
Cr     P^21_mj ~ d_m-1,j
Cr   These elements are stored in p as:
Cr               SO                          FR
Cr     P^ss_jj   = p_j(s,s)        P^ss_jj   = p_j(s,s)
Cr     P^12_j+1j = p_j(1,2)        P^12_jj+1 = p_j(1,2)  => P^12_j-1j = p_j-1(1,2)
Cr     P^21_j-1j = p_j(2,1)        P^21_jj-1 = p_j(2,1)  => P^21_j+1j = p_j+1(2,1)
Cr   Note that (Ptrans)^sj_il = P^ks_li
Cr   The Kroneker deltas imply that the sum over m contains a single term so that
Cr      sum_m g^kk'_lm P^k's'_mj = g^kk'_lj' p_j(k's')
Cr   provided j' = j when k'=s', j' = j-1 when k'>s'  and j' = j+1 when k'<s'
Cr   Similarly
Cr      sum_lm P^ks_li g^kk'_lm P^k's'_mj = p_i(ks) g^kk'_l'j' p_j(k's')
Cr   provided also l' = l when k=s, l' = l-1 when k>s  and l' = l+1 when k<s
Cr   These rules apply to the SO case.
Cr   In compact form, SO case i' = i - k + s    and   j' = j - k' + s'
Cr   For the FP case,
Cr   In compact form, FR case i' = i + k - s    and   j' = j + k' - s'
Cr   There is also a registry shift in the storage of P by p (see above)
Cr   To accommodate this, copy to local array with registry shift.
Cr
Cr   g can be a subblock of the entire hamiltonian.
Cr   At present, the subblock must digaonal part of g (ldi=ldj)
Cr   g(1,1) points to the start of the subblock
Cu Updates
Cu  08 May 16 Added offpi, offpj
Cu  11 Mar 16 (MvS) adapted from frscale
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldi,ldj,offpi,offpj,pfdim
      complex(8), target :: P(pfdim,2,2),g(ldi,2,ldj,2)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: wk(:,:,:,:)
C     complex(8), pointer :: pt(:,:,:)
C ... Local parameters
      integer mod1,s,sp,k,kp,i,j,ip,jp
      complex(8), pointer :: Pr(:,:,:), Pl(:,:,:)

      call rxx(max(ldi,ldj) > pfdim,'soscale: wrong dimensions')
      call rxx(ldi /= ldj .or. offpi /= offpj,
     .  'soscale not ready for off-diagonal subblocks')

C      call zprm('g before trans',2,g,ldi*2,ldi*2,ldj*2)
C      print *, '--- START'

      mod1 = mod(mode/10,10)

      allocate(wk(ldi,2,ldj,2)); call dpzero(wk,2*size(wk))
      if (mod1 == 0) then
        Pl => P; Pr => P
      else
        allocate(Pl(pfdim,2,2)); call zcopy(size(Pl),P,1,Pl,1)
        Pl(1,1,2) = 0; Pl(ldi,2,1) = 0
        Pl(2:ldi,1,2)   = P(1:ldi-1,1,2)
        Pl(1:ldi-1,2,1) = P(2:ldi,2,1)

        if (offpi == offpj) then
          Pr => Pl
        else
          call rx('check this branch')
C       allocate(Pr(pfdim,2,2)); call zcopy(size(Pr),P,1,Pr,1)
C       Pr(2:ldj,1,2)   = P(1:ldj-1,1,2)
C       Pr(1:ldj-1,2,1) = P(2:ldj,2,1)
        endif
      endif

C     Comments from frscale ... not sure about them?
C ... Sum rule (conserved mu: i-s = j-k) (?)
C ... (spins are numbered -1/2, 1/2; m's are l,l-1,...,-l)
      do  s = 1, 2
      do  sp = 1, 2
        do  k = 1, 2
        do  kp = 1, 2
          do  i = 1, ldi
            ip = i - k + s;  if (mod1 == 1) ip = i + k - s
            if (ip < 1 .or. ip > ldi) cycle
            do  j = 1, ldj
              jp = j - kp + sp;  if (mod1 == 1) jp = j + kp - sp
              if (jp < 1 .or. jp > ldj) cycle
              wk(i,s,j,sp) = wk(i,s,j,sp) + Pl(offpi+i,k,s) * g(ip,k,jp,kp) * Pr(offpj+j,kp,sp)

C              if (abs(Pl(offpi+i,k,s) * g(ip,k,jp,kp) * Pr(offpj+j,kp,sp)) > 1d-9)
C     .        call info5(1,0,0,'%4:1i %2;11,6D %2;11,6D %2;11,6D',[i,j,s,sp],
C     .          Pl(offpi+i,k,s),g(ip,k,jp,kp),Pr(offpj+j,kp,sp),5)


C             if (i == 6 .and. s == 1 .and. j == 6 .and. sp == 1) then
C             if (i == 6 .and. s == 1 .and. j == 17 .and. sp == 2) then
C              if (i == 1 .and. s == 1 .and. j == 2 .and. sp == 2) then
C                print 321, i,s,j,sp,Pl(i,k,s),ip,k,jp,kp,g(ip,k,jp,kp),Pr(j,kp,sp),g(ip,k,jp,kp)*Pr(j,kp,sp)
C  321           format(4i3,2f15.10,2x,4i3,6f15.10)
C              endif

            enddo !j
          enddo !i
        enddo !kp
        enddo !k
      enddo !sp
      enddo !s

      call dcopy(2*size(g),wk,1,g,1)
C     call zprm('Ptrans g P ',2,g,ldi*2,ldi*2,ldj*2)
      deallocate(wk)

      end
