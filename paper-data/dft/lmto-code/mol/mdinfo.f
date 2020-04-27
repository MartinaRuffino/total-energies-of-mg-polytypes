      subroutine mdinfo(nbas,amass,tstep,f,iclas,temp,relax,
     .                  bas,vel,zeta,zacc,zsqacc,itrlx,ekin,tkin,
     .                  etot)
C --- Print out some information about the MD.
C ---------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nbas,iclas(nbas),itrlx
      double precision amass(nbas),f(3,nbas),bas(3,nbas),
     .  vel(3,nbas),temp,ekin,tkin,etot,tstep,zeta,
     .  relax,zacc,zsqacc
C Local Variables
      integer id,np,fopn,i1mach
      double precision tzeta,tzetsq,totmom(3),totfrc(3),cons,
     .     autime,autmp
      data autime/0.048377d0/
      data autmp/0.6333328d-5/

c --- Check the health of the simulation
      zacc= zacc + zeta
      zsqacc= zsqacc + zeta*zeta
c --- 1) zeta should avarage to zero:
      tzeta= zacc / itrlx
c --- 2) mean fake kinetic energy:
      tzetsq= zsqacc / itrlx
c --- 3) conserved quantity:
      cons= etot + ekin + 3*(nbas-1)*temp*zacc*tstep
     .     + 0.5d0*3*(nbas-1)*temp*relax**2*zeta**2
c --- 4) center of mass (at the half step!) and total force:
      do 150 id=1,3
        totmom(id)=0.d0
        totfrc(id)=0.d0
        do 140 np=1,nbas
          totmom(id)= totmom(id) + amass(iclas(np))*vel(id,np)
          totfrc(id)= totfrc(id) + f(id,np)
 140    continue
 150  continue

      call awrit5(' MDINFO: iter=%i time=%;2dfs eh=%;5d ekin=%;4d'//
     .     ' sum=%;4d Ry',' ',120,i1mach(2),
     .     itrlx,itrlx*tstep*autime,etot,ekin,etot+ekin)

      call awrit3('         conserved quantity = %;4d Ry %N'//
     .            '         temperature        = %;0d K '//
     .            '  target temperature= %;0d K ','
     .       ',240,i1mach(2),cons,tkin/autmp
     .       , temp/autmp)
      call awrit6('        total momentum     = %;8d %;8d %;8d %N'//
     .            '        total force        = %;8d %;8d %;8d ','
     .       ',240,i1mach(2),totmom(1),totmom(2),totmom(3)
     .       ,totfrc(1),totfrc(2),totfrc(3))

      if (relax < 1d10)
     .       call awrit3('        zeta = %;8d '//
     .                   ' <zeta>= %;8d '//
     .                   ' <zeta>^2= %;8d ','
     .       ',240,i1mach(2),zeta,tzeta,tzetsq)

 100  format(f9.2,f10.2,6e12.3)
 200  format(f9.2,4f14.8,3f12.8)
 300  format(f9.2,9f11.6)
      end


