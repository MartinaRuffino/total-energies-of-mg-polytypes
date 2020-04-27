      subroutine potpm(lmaxp,nbas,nbasp,n0,ips,mpole,dpole,atid,
     .                 nxi,lxi,pos,nvi,pot0,qmom)
C- Potential from point multipoles; print out multipole moments
C ----------------------------------------------------------------------
Ci Inputs: lmax (maximum l of point multipoles)
Ci         nbas,nbasp (nbasp-nbas point multipoles from ctrl)
Ci         mpole,dpole (multipole moments unto l=1 from ctrl)
Ci
Co Outputs:
Co   pot0 :for sites nbas+1..nbasp is computed from the point multipoles
Co   qmom :are copied from mpole,dpole
Cr Remarks
Cr   'monopole' and 'dipole' taken from ctrl are in units of number of
Cr   electrons and a.u. of length. Here, they are converted to Jackson's
Cr   'multipole moments.'
Cr   So the monopole moment is (number of electrons)*Y_0
Cr   e.g., a positive charge of n|e| has Q_00=-n*Y0. If such a charge
Cr   is not at the origin, but at say a distance s/2 along z, then
Cr   its dipole moment is Q_10=(s/2)*sqrt(3)*Q_00. A dipole made up of
Cr   a -ve charge n|e| at 0,0,-s/2 and a +ve charge n|e| at 0,0,s/2
Cr   has dipole moment -sqrt(3/4pi)*s*n.
Cr   Remember that electrons are positive in this program!
Cu Updates
Cu   04 Jun 08 (ATP) exclude charge * pos vector from local moment printout
C ----------------------------------------------------------------------
      implicit none
C ... Passed Parameters
      integer lmaxp,nbas,nbasp,n0,nvi,ips(1),nxi(1),lxi(n0,1)
      double precision mpole(nbasp),dpole(3,nbasp),pot0(nvi),qmom(nvi),
     .                 pos(3,nbasp)
      character*8 atid(nbasp)
C ... Local parameters
      double precision qmv0(81),qmv1(81),sqmv0(81)
      integer ib,ilm,ipr,j,j0,l,lmax,ixyz(4),m,is,ie,lx
      double precision Y0,atepi,df,pi,r3,D,Cm
      data ixyz /0,2,3,1/ D /2.541748/ Cm /8.478358/

      call getpr(ipr)
      pi = 4d0*datan(1d0)
      Y0 = 1/sqrt(4d0*pi)
      atepi = 8d0*pi
      r3 = dsqrt(3d0)
      if(ipr >= 30)write(6,220)
  220 format(/' potpm: local multipole moments resolved by atom'/
     .  '   ib',2x,'spec',10x,'q',10x,'px',10x,'py',10x,'pz')

C --- Start loop over atoms ---
      j0 = 0
      call dpzero(sqmv0,81)
      do  ib = 1, nbasp
        if (ib > nbas) then
          lmax = lmaxp
        else
          lmax = 0
          is = ips(ib)
          do  ie = 1, nxi(is)
            lx = lxi(ie,is)
            lmax = max(lmax,lx)
          enddo
        endif
        ilm = 0
        df = 1d0
C       Moments for this atom, for its nucleus at the origin
        do  l = 0, min(lmax,1)
          df = df*(2*l+1)
          do   m = -l, l
            ilm = ilm+1
            j = j0+ilm
C           Poke asymtotic potential of point multipoles into pot0
            if (ib > nbas) then
              if (l == 0) then
                pot0(j) = (atepi/df)*mpole(ib)*Y0
                qmom(j) = mpole(ib)*Y0
              else
                pot0(j) = (atepi/df)*dpole(ixyz(ilm),ib)*r3*Y0
                qmom(j) = dpole(ixyz(ilm),ib)*r3*Y0
              endif
            endif
C           Asymtotic potential is pot0*H0(r).
C           pot0 * nabla^2 H0 = multipole moment * 8 pi / (2l+1)!!
            qmv0(ilm) = pot0(j) * df/atepi

C           l=1 : add R Q_R
            if (l == 1) then
              qmv1(ilm) = qmv0(ilm) + pos(ixyz(ilm),ib)*qmv0(1)*r3
            endif
            if (ib <= nbas) then
              if (l == 0) sqmv0(ilm) = sqmv0(ilm) + qmv0(ilm)
              if (l == 1) sqmv0(ilm) = sqmv0(ilm) + qmv1(ilm)
            endif
          enddo
        enddo
        j0 = j0 + (lmax+1)**2
        if(ipr >= 30)
     .    write(6,221) ib,atid(ib),qmv0(1),qmv0(4),qmv0(2),qmv0(3)
  221   format(i5,3x,a4,1x,4f12.6)
      enddo
        if (ipr >= 30) then
        if (nbas == nbasp) then
        write(6,100)
     .  sqmv0(1),sqmv0(4),sqmv0(2),sqmv0(3),
     .  -sqmv0(1)/Y0,
     .  -sqmv0(4)/(r3*Y0),-sqmv0(2)/(r3*Y0),-sqmv0(3)/(r3*Y0),
     .  -sqmv0(4)/(r3*Y0)*D,-sqmv0(2)/(r3*Y0)*D,-sqmv0(3)/(r3*Y0)*D,
     .  -sqmv0(4)/(r3*Y0)*Cm,-sqmv0(2)/(r3*Y0)*Cm,-sqmv0(3)/(r3*Y0)*Cm
        else
        write(6,110)
     .  sqmv0(1),sqmv0(4),sqmv0(2),sqmv0(3),
     .  -sqmv0(1)/Y0,
     .  -sqmv0(4)/(r3*Y0),-sqmv0(2)/(r3*Y0),-sqmv0(3)/(r3*Y0),
     .  -sqmv0(4)/(r3*Y0)*D,-sqmv0(2)/(r3*Y0)*D,-sqmv0(3)/(r3*Y0)*D,
     .  -sqmv0(4)/(r3*Y0)*Cm,-sqmv0(2)/(r3*Y0)*Cm,-sqmv0(3)/(r3*Y0)*Cm
        endif
        endif
  100   format('  --- '/' total',7x,4f12.6/'   (a.u.)',4x,4f12.6/
     .         '   Debye (10^-18 esu)',4x,3f12.6/
     .         '   10^-30 Cm',13x,3f12.6)
  110   format('  --- '/' total not PM',4f12.6/'   (a.u.)',4x,4f12.6/
     .         '   Debye (10^-18 esu)',4x,3f12.6/
     .         '   10^-30 Cm',13x,3f12.6)
      end


