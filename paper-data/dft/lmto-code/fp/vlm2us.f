      subroutine vlm2us(lmaxu,rmt,idu,lmxa,iblu,vorb,ppnl,sab,vumm)
C- Rotate vorb from (phi,phidot) to (u,s) and store in vumm
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmaxu :dimensioning parameter for U matrix
Ci   rmt   :augmentation radius, in a.u.
Ci   idu   :1s digit
Ci         : 0 no U
Ci         : 1 around MF limit; see Petukhov, PRB 67, 153106 (2003)
Ci         : 2 around FLL limit; see PRB 52, R5467 (1995)
Ci         :10s digit (See Remarks)
Ci         : 0 put U on phi only
Ci         : 1 put U on phi,phidot,LO
Ci   lmxa  :augmentation l-cutoff
Ci   vorb  :orbital-dependent potential matrices
Ci   ppnl  :potential parameters
Ci   ppnl  :NMTO pot pars (but no backward extrapolation; <phi phi>=1)
Ci         :(1)  = Not used here
Ci         :(2)  = normalization of phi (= 1 here for now)
Ci         :(3)  = rmax * log derivative of phi at rmax
Ci         :(4)  = rmax * log derivative of phidot at rmax
Ci         :(5)  = phi at rmax
Ci         :(6)  = phidot at rmax
Ci         :(7)  = normalization of phidot
Ci         :(8)  = normalization <gz|gz> for local orbital
Ci         :(9)  = overlap <phi|gz>
Ci         :(10) = overlap <phidot|gz>
Ci   sab   :matrix elts of <u,s,gz | 1 | u,s,gz>:
Ci         : (1)=<ul|1|ul>   (5)=<ul|1|gz>   (8)=<gz|1|ul>
Ci         : (2)=<ul|1|sl>   (6)=<sl|1|gz>   (9)=<gz|1|sl>
Ci         : (3)=<sl|1|ul>   (7)=<gz|1|gz>
Ci         : (4)=<sl|1|sl>
Ci         :NB sab(1:9,:,:) have same ordering as vumm(:,:,1:9)
Cio Inputs/Outputs
Cio  iblu  :index to current LDA+U block
Cio        :on input, index to last LDA+U block that was accessed
Cio        :iblu will be incremented to from blocks at this site
Co Outputs
Co   vumm  :vorb for this site in (us) representation
Co         :vumm(m1,m2,1) = <u| vorb(m1,m2) |u>
Co         :vumm(m1,m2,2) = <u| vorb(m1,m2) |s>
Co         :vumm(m1,m2,3) = <s| vorb(m1,m2) |u>
Co         :vumm(m1,m2,4) = <s| vorb(m1,m2) |s>
Co         :vumm(m1,m2,5) = <u| vorb(m1,m2) |z>
Co         :vumm(m1,m2,6) = <s| vorb(m1,m2) |z>
Co         :vumm(m1,m2,7) = <z| vorb(m1,m2) |z>
Co         :vumm(m1,m2,8) = <z| vorb(m1,m2) |u>
Co         :vumm(m1,m2,9) = <z| vorb(m1,m2) |s>
Cr Remarks
Cr  Write functions (u,s) as (u_1,u_2), and (phi,phidot) as (p_1,p_2)
Cr
Cr  Rotation connecting (u,s) to (phi,phidot) is
Cr      ( u(r) )     ( u_1 )       ( phi(r)    )     ( p_1 )
Cr      (      )  =  (     )  =  M (           ) = M (    )
Cr      ( s(r) )     ( u_2 )       ( phidot(r) )     ( p_2 )
Cr  (See potpus.f for how to consruct M.)
Cr
Cr  1. 10s digit idu=0.  Interpret vorb as
Cr     < phi | U | phi > with < phi | U | dot > = < dot | U | dot > = 0
Cr  Then
Cr    <u_i | vorb | u_j> = Sum_kl delta_k1 Delta_l1 *
Cr                         <m_ik p_k | U | m_jl p_l>
Cr                       = <m_i1  p_1 | U | m_j1 p_1>
Cr                       = m_i1 m_j1 vorb                  (1)
Cr  2. 10s digit idu=1.  Interpret vorb as matrix element of a
Cr     constant potential for a particular (lm,l'm') pair,
Cr     where <lm,l'm'> is normalized.
Cr
Cr     Then <phi_1 | vorb | phi_2> = vorb S_12, S_12 being
Cr     the overlap between phi_1 and phi_2, for any of the three phi_1
Cr     in (lm) and phi_2 in (l'm').
Cu Updates
Cu   24 Dec 10 New 10's mode for idu
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   08 Jun 05 (MvS) extended to local orbitals
Cu   30 Apr 05 Lambrecht first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmaxu,lmxa,iblu,idu(4)
      double precision rmt
      integer n0,nppn,nab
      parameter (n0=10,nppn=12,nab=9)
      double precision ppnl(nppn,n0,2),sab(nab,n0,2)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,*),
     .               vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
C ... Local parameters
      integer m1,m2,l,i,k
      double precision phi,dlphi,phip,dlphip,dphi,dphip
C     double precision r12,m21,r11,r22,det
      double precision m11,m21,det
      double precision phz,dphz
      double complex vzz,vuz,vsz,vzu,vzs

C     call prmx('vorb',vorb(1,1,2,1),2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

C --- Rotate vorb from phi,phidot basis to u,s basis ---
      do  l = 0, min(lmxa,3)
        if (mod(idu(l+1),10) /= 0) then
          iblu = iblu+1
          do  i = 1, 2
C     ... U on phi,phidot, LO; see Remarks
          if (idu(l+1) > 10) then
            do  m1 = -l, l
              do  m2 = -l, l
                do  k = 1, 9
                  vumm(m1,m2,k,i,l) = vorb(m1,m2,i,iblu)*sab(k,l+1,i)
                enddo
              enddo
            enddo
C     ... U on phi only; see Eq (1) in Remarks
          elseif (idu(l+1) > 0) then
            dlphi  = ppnl(3,l+1,i)/rmt
            dlphip = ppnl(4,l+1,i)/rmt
            phi    = ppnl(5,l+1,i)
            phip   = ppnl(6,l+1,i)
            dphi   = phi*dlphi/rmt
            dphip  = dlphip/rmt*phip
            det = phi*dphip - dphi*phip
            m11 = dphip/det
C           r12 = -dphi/det
            m21 = -phip/det
C           r22 = phi/det
            do  m1 = -l, l
              do  m2 = -l, l
C               <u|vorb|u>, <u|vorb|s>, <s|vorb|u>, <s|vorb|s> from Eq. 1
                vumm(m1,m2,1,i,l) = vorb(m1,m2,i,iblu)*m11*m11
                vumm(m1,m2,2,i,l) = vorb(m1,m2,i,iblu)*m11*m21
                vumm(m1,m2,3,i,l) = vorb(m1,m2,i,iblu)*m21*m11
                vumm(m1,m2,4,i,l) = vorb(m1,m2,i,iblu)*m21*m21

                vumm(m1,m2,5,i,l) = 0
                vumm(m1,m2,6,i,l) = 0
                vumm(m1,m2,7,i,l) = 0
                vumm(m1,m2,8,i,l) = 0
                vumm(m1,m2,9,i,l) = 0
              enddo
            enddo

            phz  = ppnl(11,l+1,i)
            dphz = ppnl(12,l+1,i)

            if (phz /= 0) then
            do  m1 = -l, l
            do  m2 = -l, l
              vzz = phz**2*vumm(m1,m2,1,i,l) +
     .              phz*dphz*(vumm(m1,m2,2,i,l)+vumm(m1,m2,3,i,l)) +
     .              dphz**2*vumm(m1,m2,4,i,l)
              vuz = - phz*vumm(m1,m2,1,i,l) - dphz*vumm(m1,m2,2,i,l)
              vsz = - phz*vumm(m1,m2,3,i,l) - dphz*vumm(m1,m2,4,i,l)
              vzu = - phz*vumm(m1,m2,1,i,l) - dphz*vumm(m1,m2,3,i,l)
              vzs = - phz*vumm(m1,m2,2,i,l) - dphz*vumm(m1,m2,4,i,l)

              vumm(m1,m2,5,i,l) = vuz
              vumm(m1,m2,6,i,l) = vsz
              vumm(m1,m2,7,i,l) = vzz
              vumm(m1,m2,8,i,l) = vzu
              vumm(m1,m2,9,i,l) = vzs

            enddo
            enddo
            endif

C            if (l == 3) then
C              print *, 'vorb before rotation isp=',i,'l=',l
C              print ('(7f8.4)'),((vorb(m1,m2,i,iblu),m2=-l,l),m1=-l,l)
C              print *, 'rotated vumm'
C              print ('(7(2f8.4,1x))'),
C     .          (((vumm(m1,m2,k,i,l),m2=-l,l),m1=-l,l),k=1,4)
C              stop
C            endif

          endif
          enddo
        endif
      enddo
      end
