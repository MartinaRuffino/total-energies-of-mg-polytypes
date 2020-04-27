      subroutine lfrog(nbas,amass,tstep,alat,f,ipc,temp,relax,
     .  bas,vel,wk,yk,zeta,ekin,tkin)
C- Velocity-Leapfrog algorithm for advancing atomic positions
C ----------------------------------------------------------------------
Ci Inputs:
Ci nbas; amass; tstep; alat; f: forces; iclas; temp: temperature (a.u.)
Ci relax: relaxation time of the thermostat (if bigger then 1e10 NO thermostat)
Ci zeta: "hoover viscosity"
Co Output:
Co bas: new positions; vel: new velocities; ekin: kinetic energy;
C
C --- The equation of motion are integrated using the Velocity-Leapfrog
C     algoritm. A Nose-Hoover thermostat is implemented so that different
C     ensembles can be simulated (NVE or NVT).
C     The velocities are calculated at the "half-step", but the kinetic
C     energy is calulated using the velocities "on the step".
C------------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nbas,ipc(nbas)
      double precision amass(nbas),f(3,nbas),bas(3,nbas),vel(3,nbas),
     .  alat,temp,ekin,tstep,zeta,wk(3,nbas),yk(3,nbas),relax
C Local Variables
      integer id,np,fopn,ibas,i,iprint
      double precision zetden,zetnum,ekin2,velsq,acc,tkin,
     .     zetnew

c --- Set up the variables for the Nose-Hoover thermostat. In case of the
c     Canonical Ensemble (NVE), "zetnum" and "zetden" are set to unity.
      if (relax > 1d10) zeta = 0d0
      zetden = 1d0 / ( 1d0 + 0.5d0*tstep*zeta )
      zetnum = 1d0 - 0.5d0*tstep*zeta

c --- Zero the accumulators
      ekin2 = 0d0
      ekin = 0d0

c --- Integrate the equations of motion.
      do  np = 1, nbas
        velsq = 0d0
        do  id = 1, 3
c ------- accelerations,velocities and positions
          acc = f(id,np) / amass(ipc(np))
          yk(id,np) = ( vel(id,np)*zetnum + tstep*acc ) * zetden
          wk(id,np) = bas(id,np) + tstep*yk(id,np)/alat

c ------- kinetc terms for thermostat and for energy
          ekin2 = ekin2 + amass(ipc(np))*yk(id,np)**2
          velsq = velsq + (0.5d0*( yk(id,np)+vel(id,np) ))**2
        enddo
c ----- kinetic energy
        ekin = ekin + 0.5d0*amass(ipc(np))*velsq
      enddo

c --- Compute the new value of the friction coefficient "zeta"
      tkin = ekin2 / (3.d0*(nbas-1))
      zetnew = zeta + tstep * (tkin/temp - 1d0) / relax**2

c --- Shuffle positions, velocities and friction
      zeta = zetnew
      do  np = 1, nbas
        do  id = 1, 3
          bas(id,np) = wk(id,np)
          vel(id,np) = yk(id,np)
        enddo
      enddo
      if (relax > 1d10) zeta = 0d0

C --- Printout
      if (iprint() >= 30) then
        write (*,1)
        write (fopn('LOG'),1)
        do  ibas = 1, nbas
          write (*,2) ibas,(bas(i,ibas),i=1,3)
          write (fopn('LOG'),2) ibas,(bas(i,ibas),i=1,3)
        enddo
      print *
      endif
    1 format(' LEAPFROG: New atom positions :')
    2 format(i4,3x,3F14.8)
      end
