      module newtvC_common
      implicit none

      integer, parameter :: NP=4
      real(8) :: fvecC(NP), goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP)
      integer :: nn

      end module newtvC_common


!C      subroutine fdpp(enu1,enu2,ksop,shft,gmt,fmt,gmtde,fmtde,
!C     .  z,rmax,avw,l,lmx,imu,
!C     .  sr1,sr2,srav1,srav2,srdot1,srdot2,sravdot1,sravdot2,  ! not used now
!C     .  gsr1,gsr2,gsrav1,gsrav2,gsrdot1,gsrdot2,gsravdot1,gsravdot2, ! not used now
!C     .  pprel,gsmt)

!     module m_rdeq
!     implicit none
!     contains


      subroutine rdeqcore(enum1,enum2,ksop,z,v,Erseq,nmrseq,rofi,nr,nsp,a,b,l,lmx, &
                   nodes,imu,scalede,avw,grseq,gn,fn)!gmt,fmt,gmtde,fmtde,gsmt,pprel,Eout,gn,fn)
!C-Solves radial dirac equations, in the presence of an effective magnetic field
!C ----------------------------------------------------------------------
!Ci Inputs
!Ci   enum1 :linearization energy, spin 1
!Ci   enum2 :linearization energy, spin 2
!Ci   ksop  :spin orbit coupling parameters, used to adjust enu1,enu2
!Ci   z     :nuclear charge
!Ci   v     :spherical potential (atomsr.f)
!Ci   rofi  :radial mesh points
!Ci   nr    :number of radial mesh points
!Ci   nsp   :2 for spin-polarized case, otherwise 1
!Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!Ci   l     :l quantum number
!Ci   lmx   :dimensions ksop and pprel
!Ci   imu   :projection of j (positive integer, counts -l-1/2...l+1/2)
!Ci   scalede: numerical factor for energy differentiation
!Ci   Erseq : energy output of a scalar case rseq
!Ci   avw   : average WS radius
!Ci   nodes : number of nodes
!Co Outputs
!Co   gmt   :(2x2) g amplitudes x r at rmax
!Co   fmt   :(2x2) f amplitudes x r at rmax
!Co   gsmt  : ?
!Co   pprel :Relativistic potential parameters in kappa-mu representation
!Co         :pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),1:2,1:2)
!Co         :Not made here:
!Co         :pprel(1,:,:,1:2,1:2) = C
!Co         :pprel(2,:,:,1:2,1:2) = gamma
!Co         :pprel(3,:,:,1:2,1:2) = delta
!Co         :Made here:
!Co         :pprel(4,:,:,1:2,1:2) = small parameter p
!Co   gmtde :(2x2 matrix) gdot at rmax
!Co   fmtde :(2x2 matrix) fdot at rmax
!Cl Local variables
!Cr Remarks
!Cr   In this version enu depends only on kappa
!Cu Updates
!Cu   18 Oct 13 (Belashchenko) cleanup, bug fix
!Cu   18 Jun 04 (A Chantis) working version
!Cu   24 Jun 03 (A Chantis) adopted from V. Antropov
!C ----------------------------------------------------------------------
      implicit none
!C ... Passed parameters
!      integer today(3), now(3)
      integer nr,nsp,imu,l,lmx,nm, nctp, nctp0, nre, nF, nc, i, ifi
      integer j1, j2, ien, ntrial, nsh, nnod,nmrseq,nodes
      integer nnug, nnugzv, nclose
      double precision v(nr,nsp),rofi(nr),z,a,b,enum1,enum2,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision pprel(4,0:lmx,2*(lmx+1),2,2)
      double precision gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),gsmt(2,2)
      double precision grseq(nr,2)
!C ... Local parameters
      logical check, maxi
      integer n,np,nrmx,ie,alpha,a1,ncut
!C     integer k1,k2,i,alp,ifi,fopn,ie
!C     character*(10) outs*80, fil*12
      double precision J, &   !Jacobian
       c,csq,            &   !Speed of light 274.072d0
       de,               &   !Finite energy difference for num diff.
       eav,enu1,enu2, gamma,hoc,mu,norm,pr,q,r,rm,  &  !,prodr,prodr2
       r1,rmax,smalpa(2,2),sqru,u,up,upp,w0,w1,w2,w3,x1,x2,xk,   &
       xk1,ff, e1, e2, re, En, enu1my, enu2my, E00

      double precision  g1out1(nr),g2out1(nr),f1out1(nr),f2out1(nr),   &
       g1sout1(nr),g2sout1(nr),f1sout1(nr),f2sout1(nr),g1in1(nr),      &
       g2in1(nr),f1in1(nr),f2in1(nr),g1sin1(nr),g2sin1(nr),f1sin1(nr), &
       f2sin1(nr),nug1(nr),nug2(nr) ,nuf1(nr),nuf2(nr),bm(nr),         &
       vm(nr),psiout1(nr,4,2,5),psiin1(nr,4,2,5),psd(nr,4,2),Am, Aout, &
       Ain, psiout2(nr,4,2,5),psiin2(nr,4,2,5),psi(nr,4,2,5),          &
       F(4),Pv(4),Erseq,                                               &
       g1out2(nr),g2out2(nr),f1out2(nr),f2out2(nr),                    &
       g1sout2(nr),g2sout2(nr),f1sout2(nr),f2sout2(nr),g1in2(nr),      &
       g2in2(nr),f1in2(nr),f2in2(nr),g1sin2(nr),g2sin2(nr),f1sin2(nr), &
       f2sin2(nr), g1final(nr),gfinal(nr),gfinalb(nr),ffinal(nr),      &
       ffinalb(nr),fvecC(4),Pv0(4),fvec0(4),Pvpl(4),Pv02(4),fsave(4),   &
       g1finalE(nr),gfinalE(nr),gfinalbE(nr),ffinalE(nr),ffinalbE(nr)

      double precision gn(2,2,nr), fn(2,2,nr),gn1(2,2,nr),&
       fn1(2,2,nr), gn2(2,2,nr), fn2(2,2,nr), Pvsh(4)

      double precision DeltaE, Eni(5), gp(2,2,5,nr), fp(2,2,5,nr),     &
       Eip(5), fdot(2,2,nr), gdot(2,2,nr), E0, sr, enu, avw, coeff

      double precision gdotT(2,2), gnr(2,2), fnr(2,2), gnrT(2,2),      &
      Temp1(2,2), Temp2(2,2),  Anu(2,2), kappa(2,2), Temp3(2,2),       &
      Unit2(2,2),gammal(2,2), Cl(2,2), Wl(2,2), Dnu(2,2),Fsq

      double precision tolx, tolf

      double precision Epl(1:61), Delpl, Npl

      double precision psifinal(nr,4,2),psifinalE(nr,4), normf, normfE

      procedure(real(8)) :: prodr,prodr2

      double precision  g1out10(nr),g2out10(nr),f1out10(nr),f2out10(nr),    &
       g1sout10(nr),g2sout10(nr),f1sout10(nr),f2sout10(nr),g1in10(nr),      &
       g2in10(nr),f1in10(nr),f2in10(nr),g1sin10(nr),g2sin10(nr),f1sin10(nr),&
       f2sin10(nr),psiout10(nr,4,2,5),psiin10(nr,4,2,5),psiout20(nr,4,2,5), &
       psiin20(nr,4,2,5),psi0(nr,4,2,5),                                    &
       g1out20(nr),g2out20(nr),f1out20(nr),f2out20(nr),                     &
       g1sout20(nr),g2sout20(nr),f1sout20(nr),f2sout20(nr),g1in20(nr),      &
       g2in20(nr),f1in20(nr),f2in20(nr),g1sin20(nr),g2sin20(nr),f1sin20(nr),&
       f2sin20(nr)

      double precision Pvsh2(4),E1sh, E2sh, Ash, dEsh, Efact, Vtot,normsave1,normsave2
      double precision Vnug(nr), Vnugzv(nr),dVnug(nr)
      double precision Eout((lmx+1),(2*lmx+2)),product


!C ... External calls
      external rx !, fctp0, fctp !dinv22,mul22,
!C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      eav = (enum1 + enum2)/2.0_8
!      enu1 =  -(l+1) * ksop(l,1,1,1)/2.0_8 + eav
!      enu2 =  l    * ksop(l,1,1,1)/2.0_8 + eav

      mu  = imu - l - 1.5d0


      if (z == 0) return

      open(1000,FILE='fvecs.dat')
      csq = c*c
      nrmx = nr
      hoc = (2*z/c)
      de = 0.03d0*scalede ;
      ff = 1
      nc = nmrseq

      rmax = rofi(nr)
      vm(:) = (v(:,1) + v(:,2))/2
      bm(:) = ((v(:,2) - v(:,1))/2)

      if (imu == 1 .or. imu == 2*l+2) then
        xk = -l-1
      else
        xk = l
      endif
      xk1 = -l-1


      kappa(1,2) = 0.0_8; kappa(2,1)=0.0_8
      if (imu == 1 .or. imu == 2*l+2) then
        kappa(1,1) = -l-1.0_8
      else
        kappa(1,1) = l
      endif
      kappa(2,2) = -l-1.0_8

      up = mu/(l-0.5d0) ; u = mu/(l+0.5d0) ; upp = mu/(l+1.5d0)
      sqru = dsqrt(1.0_8-u*u)

!C ... Added to take care of the abs(mu)=l+0.5 case
      if (imu == 1 .or. imu == 2*l+2) up = upp

        ie = 1
!        e2 = enu2
!        e1 = enu1

!        En = Erseq

        En = eav
!        if ((l==2).and.(imu==5)) En = -5

!________________Adjusting the number of nodes_______________________
!C --Adjusting the number of nodes in the inward integrated g11
!C -- since it's supposed to be eq. to n-l-1 (or lmx-l). Number
!C -- decreases with energy (which is negative)


!        call integrationC(1, En, xk, xk1, mu, up, u,                      &
!            upp, sqru, vm,                                                &
!            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,              &
!            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,          &
!            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,        &
!            f2sin1,psiout1, psiin1, nc, nnod)

!        if (nnod /= (lmx-l)) then
!           Efact = 1.1d0
!           do while (nnod /= (lmx-l))
!               if (nnod > (lmx-l)) then
!                    En=En*Efact
!                  else
!                    En=En/Efact
!               endif
!               call integrationC(1, En, xk, xk1, mu, up, u,               &
!                upp, sqru, vm,                                            &
!                bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,          &
!                g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,      &
!                f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,    &
!                f2sin1,psiout1, psiin1, nc, nnod)
!           enddo
!        endif
!        print*,'Afder NodAdjust E = ',En,'; Num of nods =',nnod

!        if ((l==0).and.(imu==1)) then
!           ifi = 5555
!           write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!             write(ifi,'(9(x,g16.8))') rofi(n),g1out1(n),g1in1(n),      &
!                  g2out2(n),g2in2(n),f1out1(n),f1in1(n),f2out2(n),      &
!                  f2in2(n)
!           enddo
!           close(ifi)
!        endif

!_______________End adjusting the number of nodes_________________________

        print*,'imu =',imu,' l = ',l,'nc = ',nmrseq

!C---First intergation to get A=g11out(nc)/g11in(nc):

!        En = eav
        En = Erseq

        call integrationC(1, En, xk, xk1, mu, up, u,               &
           upp, sqru, vm,                                          &
           bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede, grseq, &
           g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,    &
           f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,  &
           f2sin1, nmrseq, nnod,normsave1)
!        print*,'nnod (alpha=1) = ',nnod,'rseq: ',nodes

        call integrationC(2, En, xk, xk1, mu, up, u,               &
           upp, sqru, vm,                                          &
           bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede, grseq, &
           g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
           f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
           f2sin2, nmrseq, nnod,normsave2)
!        print*,'nnod (alpha=2) = ',nnod,'rseq: ',nodes

!_________________Print_____________________

        ifi = 1000+100*l+imu
        write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
        write(ifi,*)'# l:',l,'imu:',imu
        do n=1,nr
              write(ifi,'(9(x,g16.8))') rofi(n),g1out1(n),g1in1(n),  &
                   g2out2(n),g2in2(n),f1out1(n),f1in1(n),f2out2(n),  &
                   f2in2(n)
        enddo
        close(ifi)

!__________________End Print_________________________________

!_________________Print_____________________

!        ifi = 1000*INT(abs(En))+100*l+imu
!        write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
!        write(ifi,*)'# l:',l,'imu:',imu
!        do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),g1out1(n),g1in1(n),  &
!                   g2out2(n),g2in2(n),f1out1(n),f1in1(n),f2out2(n),  &
!                   f2in2(n)
!        enddo
!        close(ifi)

!__________________End Print_________________________________

!C--Building initial vector P:-------------------------------
        nc = nmrseq
        Am = g1out1(nc)/g1in1(nc)
!        Aout = 0.001d0
!        Ain = 0.001d0
        Aout= 0.000001d0
        Ain = Aout*g2out2(nc)/Am/g2in2(nc)


!        Aout = Aout*(1+0.1)
!        Ain = Ain*(1+0.1)

        Pv(1)=En; Pv(2)=Am; Pv(3)=Aout; Pv(4)=Ain

        print*,'E initial in rdeqcore =',Pv(1)
        Pv0(:) = Pv(:)
        Pv02(:) = Pv(:)

!________________________First F=0 solver_______________________
!C-- Seems to be working worse than Broydn

 214    nF = 4;
        call newtC(Pv, nF, xk, xk1, mu, up, u,                   &
             upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,     &
             lmx, imu,scalede, grseq,check,fvecC,maxi)
        nc = nmrseq

        print*,'fvecC = ',fvecC
        print*,'E after newt = ',Pv(1),' nmrseq = ',nmrseq,'nc = ',nc


!____________________End First F=0 solver______________________

!__________Alternative F=0 solver__________________________________
!C--   mnewt from Numerical Recipes. Often finds solution too far   &
!C--   from the initial (scalar) energy value, though sometimes is  &
!C--   more effective.

       if (maxi==.true.) then
          ntrial = 100
          tolf = 1e-5_8
          tolx = 1e-5_8
          call mnewtC(ntrial,Pv0,nF,tolx,tolf, xk, xk1, mu, up, u,    &
            upp, sqru, vm,                                         &
            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq,fvec0)

          Pv(:)=Pv0(:)
          print*,'ALTERN METHOD F = ',fvec0
          print*,'E = ',Pv0(1),' nc = ',nc
!          call pause('MAXITS exceeded in newtC')
       endif
!      call pause('MAXITS exceeded in newtC')
!__________End Alternative F=0 solver________________________________


!C-------Solving F=0, Pv02 - input and output; pv(1)-energy---------
!C-------pv(2)=A, pv(3)=Aout, pv(4)=Aind----------------------------
!      Pv02(:) = Pv(:)

!      call broydnC(Pv02,nF,check, xk, xk1, mu, up, u,          &
!            upp, sqru, vm,                                    &
!            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,  grseq, &
!            fsave)
        nc = nmrseq
!        print*,' '
!        print*,'BROYDN F = ',fsave
!        print*,'E = ',Pv02(1),' nc = ',nc
!        Pv(:) = Pv02(:)
!        print*,' '
!        print*,'Is this a local minimum: ',check

!__________________Book approximation______________________________
!        nsh = 1
!         E1sh=eav
!!        E1sh = En
!        dEsh=1 ! Random above 1e-5
!        do while ((dEsh>=1e-5).and.(nsh<200))
!           call integrationC(1, E1sh, xk, xk1, mu, up, u,             &
!              upp, sqru, vm,                                        &
!              bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,      &
!              g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,  &
!              f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,&
!              f2sin1,psiout1, psiin1, nc, nnod)

!           call integrationC(2, E1sh, xk, xk1, mu, up, u,             &
!              upp, sqru, vm,                                        &
!              bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,      &
!              g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,  &
!              f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,&
!              f2sin2, psiout2, psiin2, nc, nnod)
!           Ash = (g1out1(nc)/g1in1(nc))
!           E2sh = E1sh + abs(g1out1(nc))*(Ash*f1in1(nc)-f1out1(nc))
!           nsh = nsh + 1
!           dEsh = dabs(E2sh - E1sh)
!           E1sh = E2sh
!           if ((nsh >= 197).or.(nsh <= 3)) print*,'dE = ', dEsh
!        enddo
!        Ash = g1out1(nc)/g1in1(nc)
!        print*,'BOOK APP: E = ',E2sh,' nsh = ',nsh

!        Aout = 0.001d0
!        Ain = 0.001d0

!        Pvsh(1)=E2sh; Pvsh(2)=Ash; Pvsh(3)=Aout; Pvsh(4)=Ain
!        Pvsh2 = Pvsh

!        call broydnC(Pvsh,nF,check, xk, xk1, mu, up, u,          &
!            upp, sqru, vm,                                      &
!            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,  grseq, &
!            fsave)

!        print*,' '
!        print*,'BROYDN BOOK F = ',fsave
!        print*,'E = ',Pvsh(1)
!        print*,' '
!        print*,'Is this a local minimum: ',check

!        call mnewt(ntrial,Pvsh2,nF,tolx,tolf, xk, xk1, mu, up, u,    &
!            upp, sqru, vm,                                         &
!            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,fvec0)

!        print*,' '
!        print*,'ALTERN METHOD F = ',fvec0
!        print*,'E = ',Pvsh2(1)
!        print*,' '

!__________________End Book approximation__________________________

!C-----Integrating with the energy found by Broydn: ---------------
!         Pv(:)=Pvsh(:)
!         Pv(:)=Pv02(:)
         Eout(l+1,imu) = Pv(1)
!         print*,'Eout =',Eout
         nc = nmrseq
!         e2 =  -(l+1) * ksop(l,1,1,1)/2.0_8 + Pv(1)
!         e1 =  l    * ksop(l,1,1,1)/2.0_8 + Pv(1)

         call integrationC(1, Pv(1), xk, xk1, mu, up, u,            &
            upp, sqru, vm,                                          &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede, grseq, &
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,    &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,  &
            f2sin1, nmrseq, nnod,normsave1)

         call integrationC(2, Pv(1), xk, xk1, mu, up, u,            &
            upp, sqru, vm,                                          &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede, grseq, &
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
            f2sin2, nmrseq, nnod,normsave2)

!         print*,'ln407, g1in1(nr) after F=0:', g1in1(nr)

!_________________Print_____________________


!           ifi = 2000+100*l+imu
!           write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),g1out1(n),g1in1(n),  &
!                   g2out2(n),g2in2(n),f1out1(n),f1in1(n),f2out2(n),  &
!                   f2in2(n)
!           enddo
!           close(ifi)

!__________________End Print_________________________________


!C---Building functions from Inwards and Outwars:--------------------


!       gfinal(1:nc)=(g1out1(1:nc)+Pv(3)*g1out2(1:nc))
!       gfinalb(1:nc)=(g2out1(1:nc)+Pv(3)*g2out2(1:nc))
!       ffinal(1:nc)=(f1out1(1:nc)+Pv(3)*f1out2(1:nc))
!       ffinalb(1:nc)=(f2out1(1:nc)+Pv(3)*f2out2(1:nc))

!       gfinal((nc+1):nr)=(g1in1((nc+1):nr)+Pv(4)*g1in2((nc+1):nr))*Pv(2)
!       gfinalb((nc+1):nr)=(g2in1((nc+1):nr)+Pv(4)*g2in2((nc+1):nr))*Pv(2)
!       ffinal((nc+1):nr)=(f1in1((nc+1):nr)+Pv(4)*f1in2((nc+1):nr))*Pv(2)
!       ffinalb((nc+1):nr)=(f2in1((nc+1):nr)+Pv(4)*f2in2((nc+1):nr))*Pv(2)

!_________________Print_____________________


!           ifi = 3000+100*l+imu
!           write(ifi,*)'r,gfin,ffin,gfin2,ffin2,same'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),gfinal(n),ffinal(n),     &
!                  gfinalb(n),ffinalb(n),gfinal(n),ffinal(n),gfinalb(n), &
!                   ffinalb(n)
!           enddo
!           close(ifi)

!__________________End Print_________________________________


!_______________________Normalization of the whole___________________

!       psifinal(:,1)=gfinal(:)
!       psifinal(:,2)=ffinal(:)
!       psifinal(:,3)=gfinalb(:)
!       psifinal(:,4)=ffinalb(:)
!       call productn4(psifinal(:,:),psifinal(:,:),a,b,rofi,nr,normf)
!       do j1=1,4
!          psifinal(:,j1)=psifinal(:,j1)/dsqrt(normf)
!       enddo

!___________________End Normalization of the whole___________________


!       call printdata2('gfinal', 'rofi gfinal 12 02 2014', nr-1, rofi(2:), gfinal(2:))
!       call printdata2('ffinal', 'rofi ffinal 29 01 2014', nr-1, rofi(2:), ffinal(2:))
!       call printdata2('gfinalb', 'rofi gfinalb 29 01 2014', nr-1, rofi(2:), gfinalb(2:))
!       call printdata2('ffinalb', 'rofi ffinalb 29 01 2014', nr-1, rofi(2:), ffinalb(2:))

!________________Sewing the inward/outward____________________________

       nc = nmrseq

!       gn(1,1,1:nc)=g1out1(1:nc)
!       gn(2,2,1:nc)=g2out2(1:nc)
!       gn(2,1,1:nc)=g2out1(1:nc)
!       gn(1,2,1:nc)=g1out2(1:nc)
!       fn(1,1,1:nc)=f1out1(1:nc)
!       fn(2,2,1:nc)=f2out2(1:nc)
!       fn(2,1,1:nc)=f2out1(1:nc)
!       fn(1,2,1:nc)=f1out2(1:nc)

!       gn(1,1,(nc+1):nr)=g1in1((nc+1):nr)*Pv(2)
!       gn(2,2,(nc+1):nr)=g2in2((nc+1):nr)*Pv(2)*Pv(4)/Pv(3)
!       gn(2,1,(nc+1):nr)=g2in1((nc+1):nr)*Pv(2)*Pv(4)/Pv(3)
!       gn(1,2,(nc+1):nr)=g1in2((nc+1):nr)*Pv(2)
!       fn(1,1,(nc+1):nr)=f1in1((nc+1):nr)*Pv(2)
!       fn(2,2,(nc+1):nr)=f2in2((nc+1):nr)*Pv(2)*Pv(4)/Pv(3)
!       fn(2,1,(nc+1):nr)=f2in1((nc+1):nr)*Pv(2)*Pv(4)/Pv(3)
!       fn(1,2,(nc+1):nr)=f1in2((nc+1):nr)*Pv(2)

       gn(1,1,1:nc)=g1out1(1:nc)
       gn(2,2,1:nc)=g2out2(1:nc)*Pv(3)
       gn(2,1,1:nc)=g2out1(1:nc)
       gn(1,2,1:nc)=g1out2(1:nc)*Pv(3)
       fn(1,1,1:nc)=f1out1(1:nc)
       fn(2,2,1:nc)=f2out2(1:nc)*Pv(3)
       fn(2,1,1:nc)=f2out1(1:nc)
       fn(1,2,1:nc)=f1out2(1:nc)*Pv(3)

       gn(1,1,(nc+1):nr)=g1in1((nc+1):nr)*Pv(2)
       gn(2,2,(nc+1):nr)=g2in2((nc+1):nr)*Pv(2)*Pv(4)
       gn(2,1,(nc+1):nr)=g2in1((nc+1):nr)*Pv(2)
       gn(1,2,(nc+1):nr)=g1in2((nc+1):nr)*Pv(2)
       fn(1,1,(nc+1):nr)=f1in1((nc+1):nr)*Pv(2)
       fn(2,2,(nc+1):nr)=f2in2((nc+1):nr)*Pv(2)*Pv(4)
       fn(2,1,(nc+1):nr)=f2in1((nc+1):nr)*Pv(2)
       fn(1,2,(nc+1):nr)=f1in2((nc+1):nr)*Pv(2)

!_______________________Normalization of the whole___________________

       psifinal(1:nr,1,1)=gn(1,1,1:nr)
       psifinal(1:nr,2,1)=fn(1,1,1:nr)
       psifinal(1:nr,3,1)=gn(2,1,1:nr)
       psifinal(1:nr,4,1)=fn(2,1,1:nr)
       psifinal(1:nr,1,2)=gn(1,2,1:nr)
       psifinal(1:nr,2,2)=fn(1,2,1:nr)
       psifinal(1:nr,3,2)=gn(2,2,1:nr)
       psifinal(1:nr,4,2)=fn(2,2,1:nr)

       do j1=1,2
         call productn4(psifinal(:,:,j1),psifinal(:,:,j1),a,b,rofi,nr,normf)
         do j2=1,4
            psifinal(:,j2,j1)=psifinal(:,j2,j1)/dsqrt(normf)
         enddo
       enddo

!--------------------Orthoganalization--------------------------------

       call productn4(psifinal(:,:,1),psifinal(:,:,2),a,b,rofi,nr,product)
       do j1=2,nr
          psifinal(j1,1,2) = psifinal(j1,1,2) - product*psifinal(j1,1,1)
          psifinal(j1,2,2) = psifinal(j1,2,2) - product*psifinal(j1,2,1)
          psifinal(j1,3,2) = psifinal(j1,3,2) - product*psifinal(j1,3,1)
          psifinal(j1,4,2) = psifinal(j1,4,2) - product*psifinal(j1,4,1)
       enddo

!----------------------Multiply the second solution (alpha=2) by
!---------------------- a phase factor (-1)

       if ((psifinal(nr,1,1)*psifinal(nr,3,2)) < 0) then
          psifinal(:,:,2)=-psifinal(:,:,2)
       endif

       gn(1,1,1:nr) = psifinal(1:nr,1,1)        ! gn(#,alpha,1:nr)
       fn(1,1,1:nr) = psifinal(1:nr,2,1)
       gn(2,1,1:nr) = psifinal(1:nr,3,1)
       fn(2,1,1:nr) = psifinal(1:nr,4,1)
       gn(1,2,1:nr) = psifinal(1:nr,1,2)
       fn(1,2,1:nr) = psifinal(1:nr,2,2)
       gn(2,2,1:nr) = psifinal(1:nr,3,2)
       fn(2,2,1:nr) = psifinal(1:nr,4,2)

!_________________Print_____________________


           ifi = 5000+100*l+imu
           write(ifi,*)'r,gfin,ffin,gfin2,ffin2,same'
           write(ifi,*)'# l:',l,'imu:',imu
           do n=1,nr
              write(ifi,'(9(x,g16.8))') rofi(n),gn(1,1,n),fn(1,1,n),     &
                  gn(2,2,n),fn(2,2,n),gn(1,2,n),fn(1,2,n),gn(2,1,n), &
                   fn(2,1,n)
           enddo
           close(ifi)

!__________________End Print_________________________________

!      call pause('One step done')

!___________________End Normalization, orthog.________________________


!       do j1=1,2
!          do j2=1,2
!             norm = dsqrt(prodr2(gn(j1,j2,:),gn(j1,j2,:),a,b,rofi,nr))
!             gn1(j1,j2,:) = gn(j1,j2,:)/norm
!             norm = dsqrt(prodr2(fn(j1,j2,:),fn(j1,j2,:),a,b,rofi,nr))
!             fn1(j1,j2,:) = fn(j1,j2,:)/norm
!          enddo
!       enddo

!_______________Normalization____________________________________________

!       do j1=1,2
!          do j2=1,2
!             norm = dsqrt(prodr(0,gn(1,j2,1),gn(1,j2,1),a,b,rofi,nr))
!             gn(j1,j2,:) = gn(j1,j2,:)/norm
!             fn(j1,j2,:) = fn(j1,j2,:)/norm
!          enddo
!       enddo

!C____________Orthonormalize alpha 1 and 2_____??????????????????????????

!       pr = prodr(0,gn(1,1,1),gn(1,2,1),a,b,rofi,nr)
!       gn(:,2,:) = gn(:,2,:) - pr*gn(:,1,:)
!       fn(:,2,:) = fn(:,2,:) - pr*gn(:,1,:)
!       norm = dsqrt(prodr(0,gn(1,2,1),gn(1,2,1),a,b,rofi,nr))
!       gn(:,2,:) = gn(:,2,:)/norm
!       fn(:,2,:) = fn(:,2,:)/norm

!       call printdata2('gnorm2', 'rofi gnorm 19 02 2014', nr-1, rofi(2:), gn(1,1,2:))

!_______________Normalization end_________________________________________


!          gnr(1,1) = gn(1,1,nr); gnr(1,2) = gn(1,2,nr)
!          gnr(2,1) = gn(2,1,nr); gnr(2,2) = gn(2,2,nr)
!          fnr(1,1) = fn(1,1,nr); fnr(1,2) = fn(1,2,nr)
!          fnr(2,1) = fn(2,1,nr); fnr(2,2) = fn(2,2,nr)

!          gmt(:,:) = gnr(:,:)*rofi(nr)
!          fmt(:,:) = fnr(:,:)*rofi(nr)

!_______________Energy derivative__________________________________

!C-------- Calculating gs, fs at 5 different energy points---------
!C-------- to find energy derivatives------------------------------


!       DeltaE = 5e-3_8
!       E0 = Erseq
!       E0 = Pv(1)
!       enu = E0

!       Eni(1) = E0 - 2*DeltaE
!       Eni(2) = E0 - DeltaE
!       Eni(3) = E0
!       Eni(4) = E0 + DeltaE
!       Eni(5) = E0 + 2*DeltaE

!        do ien = 1,5
!          call integrationC(1, Eni(ien), xk, xk1, mu, up, u,       &
!            upp, sqru, vm,                                         &
!            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,       &
!            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,   &
!            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1, &
!            f2sin1,psiout1, psiin1, nc, nnod)
!          call integrationC(2, Eni(ien), xk, xk1, mu, up, u,       &
!            upp, sqru, vm,                                         &
!            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,       &
!            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,   &
!            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2, &
!            f2sin2, psiout2, psiin2, nc, nnod)


!           Am = g1out1(nc)/g1in1(nc)
!           Aout = 0.01d0
!           Ain = 0.01d0

!           Pv(1)=Eni(ien); Pv(2)=Am; Pv(3)=Aout; Pv(4)=Ain

!           nF = 4;
!           call newt(Pv, nF, xk, xk1, mu, up, u,                   &
!            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,         &
!            lmx, imu,scalede,check,fvecC)

!           Fsq = sqrt(fvecC(1)**2+fvecC(2)**2+fvecC(3)**2+fvecC(4)**2)
!           print*,'MODULE F = ', fsq,' , E =',Pv(1)

!           call broydnC(Pv,nF,check, xk, xk1, mu, up, u,           &
!             upp, sqru, vm,                                       &
!             bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,  grseq, &
!             fsave)

!!          print*,' F = '
!!          print*,fvecC
!          write(1000,*)' F = ',fvecC(:)
!!          print*,'Is this a local minimum2: ',check
!!          print*,' '
!!          print*, 'Energy=',Pv(1),' A = ',Pv(2),
!!     . ' Aout = ', Pv(3),' Ain = ', Pv(4)

!           call integrationC(1, Pv(1), xk, xk1, mu, up, u,          &
!            upp, sqru, vm,                                          &
!            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
!            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,    &
!            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,  &
!            f2sin1,psiout1, psiin1, nc, nnod)

!           call integrationC(2, Pv(1), xk, xk1, mu, up, u,          &
!            upp, sqru, vm,                                          &
!            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
!            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
!            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
!            f2sin2, psiout2, psiin2, nc, nnod)

!_______________________Normalization of the whole___________________

!           gfinalE(1:nc)=(g1out1(1:nc)+Pv(3)*g1out2(1:nc))
!           gfinalbE(1:nc)=(g2out1(1:nc)+Pv(3)*g2out2(1:nc))
!           ffinalE(1:nc)=(f1out1(1:nc)+Pv(3)*f1out2(1:nc))
!           ffinalbE(1:nc)=(f2out1(1:nc)+Pv(3)*f2out2(1:nc))

!           gfinalE((nc+1):nr)=(g1in1((nc+1):nr)+Pv(4)*              &
!                g1in2((nc+1):nr))*Pv(2)
!           gfinalbE((nc+1):nr)=(g2in1((nc+1):nr)+Pv(4)*             &
!                g2in2((nc+1):nr))*Pv(2)
!           ffinalE((nc+1):nr)=(f1in1((nc+1):nr)+Pv(4)*              &
!                f1in2((nc+1):nr))*Pv(2)
!           ffinalbE((nc+1):nr)=(f2in1((nc+1):nr)+Pv(4)*             &
!                f2in2((nc+1):nr))*Pv(2)
!           psifinalE(:,1)=gfinalE(:)
!           psifinalE(:,2)=ffinalE(:)
!           psifinalE(:,3)=gfinalbE(:)
!           psifinalE(:,4)=ffinalbE(:)
!           call productn4(psifinalE(:,:),psifinalE(:,:),a,b,rofi,   &
!                nr,normfE)

!___________________End Normalization of the whole___________________

!           gn(1,1,1:nc)=g1out1(1:nc)      /dsqrt(normfE)
!           gn(2,2,1:nc)=g2out2(1:nc)*Pv(3)/dsqrt(normfE)
!           gn(2,1,1:nc)=g2out1(1:nc)      /dsqrt(normfE)
!           gn(1,2,1:nc)=g1out2(1:nc)*Pv(3)/dsqrt(normfE)
!           fn(1,1,1:nc)=f1out1(1:nc)      /dsqrt(normfE)
!           fn(2,2,1:nc)=f2out2(1:nc)*Pv(3)/dsqrt(normfE)
!           fn(2,1,1:nc)=f2out1(1:nc)      /dsqrt(normfE)
!           fn(1,2,1:nc)=f1out2(1:nc)*Pv(3)/dsqrt(normfE)

!           gn(1,1,(nc+1):nr)=g1in1((nc+1):nr)*Pv(2)      /dsqrt(normfE)
!           gn(2,2,(nc+1):nr)=g2in2((nc+1):nr)*Pv(4)*Pv(2)/dsqrt(normfE)
!           gn(2,1,(nc+1):nr)=g2in1((nc+1):nr)*Pv(2)      /dsqrt(normfE)
!           gn(1,2,(nc+1):nr)=g1in2((nc+1):nr)*Pv(4)*Pv(2)/dsqrt(normfE)
!           fn(1,1,(nc+1):nr)=f1in1((nc+1):nr)*Pv(2)      /dsqrt(normfE)
!           fn(2,2,(nc+1):nr)=f2in2((nc+1):nr)*Pv(4)*Pv(2)/dsqrt(normfE)
!           fn(2,1,(nc+1):nr)=f2in1((nc+1):nr)*Pv(2)      /dsqrt(normfE)
!           fn(1,2,(nc+1):nr)=f1in2((nc+1):nr)*Pv(4)*Pv(2)/dsqrt(normfE)



!C____________Orthonormalize alpha 1 and 2_____??????????????????????????


!           pr = prodr(0,gn(1,1,1),gn(1,2,1),a,b,rofi,nr)
!           gn(:,2,:) = gn(:,2,:) - pr*gn(:,1,:)
!           fn(:,2,:) = fn(:,2,:) - pr*gn(:,1,:)


!           do j1=1,2
!              psi(:,1,j1,ie) = gn(1,j1,:)
!              psi(:,2,j1,ie) = fn(1,j1,:)
!              psi(:,3,j1,ie) = gn(2,j1,:)
!              psi(:,4,j1,ie) = fn(2,j1,:)
!           enddo

!           gp(1,1,ien,:) = gn(1,1,:); gp(1,2,ien,:) = gn(1,2,:)
!           gp(2,1,ien,:) = gn(2,1,:); gp(2,2,ien,:) = gn(2,2,:)
!           fp(1,1,ien,:) = fn(1,1,:); fp(1,2,ien,:) = fn(1,2,:)
!           fp(2,1,ien,:) = fn(2,1,:); fp(2,2,ien,:) = fn(2,2,:)
!           Eip(ien) = Pv(1)

!        enddo

!C____Energy Derivatives of g and f at nr__________________________

!        do j1 = 1,2
!           do j2 = 1,2

!              gdot(j1,j2,:) = gp(j1,j2,1,:)*(Eip(3)-Eip(2))                   &
!     *(Eip(3)-Eip(4))*                                                        &
!     (Eip(3)-Eip(5))/(Eip(1)-Eip(2))/(Eip(1)-Eip(3))/(Eip(1)-Eip(4))/         &
!     (Eip(1)-Eip(5)) +                                                        &
!                           gp(j1,j2,2,:)*(Eip(3)-Eip(1))*                     &
!      (Eip(3)-Eip(4))*                                                        &
!     (Eip(3)-Eip(5))/(Eip(2)-Eip(1))/(Eip(2)-Eip(3))/(Eip(2)-Eip(4))/         &
!     (Eip(2)-Eip(5))+                                                         &
!                           gp(j1,j2,3,:)*((Eip(3)-Eip(2))*                    &
!     (Eip(3)-Eip(4))*                                                         &
!     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(4))*                       &
!     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
!     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
!     (Eip(3)-Eip(4)))                                                         &
!     /(Eip(3)-Eip(1))/(Eip(3)-Eip(2))/(Eip(3)-Eip(4))/                        &
!     (Eip(3)-Eip(5))+                                                         &
!                           gp(j1,j2,4,:)*(Eip(3)-Eip(1))                      &
!      *(Eip(3)-Eip(2))*                                                       &
!     (Eip(3)-Eip(5))/(Eip(4)-Eip(1))/(Eip(4)-Eip(2))/(Eip(4)-Eip(3))/         &
!     (Eip(4)-Eip(5))+                                                         &
!                           gp(j1,j2,5,:)*(Eip(3)-Eip(1))                      &
!      *(Eip(3)-Eip(2))*                                                       &
!     (Eip(3)-Eip(4))/(Eip(5)-Eip(1))/(Eip(5)-Eip(2))/(Eip(5)-Eip(3))/         &
!     (Eip(5)-Eip(4))


!             fdot(j1,j2,:) = fp(j1,j2,1,:)*(Eip(3)-Eip(2))                    &
!      *(Eip(3)-Eip(4))*                                                       &
!     (Eip(3)-Eip(5))/(Eip(1)-Eip(2))/(Eip(1)-Eip(3))/(Eip(1)-Eip(4))/         &
!     (Eip(1)-Eip(5)) +                                                        &
!                           fp(j1,j2,2,:)*(Eip(3)-Eip(1))                      &
!      *(Eip(3)-Eip(4))*                                                       &
!     (Eip(3)-Eip(5))/(Eip(2)-Eip(1))/(Eip(2)-Eip(3))/(Eip(2)-Eip(4))/         &
!     (Eip(2)-Eip(5))+                                                         &
!                           fp(j1,j2,3,:)*((Eip(3)-Eip(2))                     &
!      *(Eip(3)-Eip(4))*                                                       &
!     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(4))*                       &
!     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
!     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
!     (Eip(3)-Eip(4)))                                                         &
!     /(Eip(3)-Eip(1))/(Eip(3)-Eip(2))/(Eip(3)-Eip(4))/                        &
!     (Eip(3)-Eip(5))+                                                         &
!                           fp(j1,j2,4,:)*(Eip(3)-Eip(1))                      &
!      *(Eip(3)-Eip(2))*                                                       &
!     (Eip(3)-Eip(5))/(Eip(4)-Eip(1))/(Eip(4)-Eip(2))/(Eip(4)-Eip(3))/         &
!     (Eip(4)-Eip(5))+                                                         &
!                           fp(j1,j2,5,:)*(Eip(3)-Eip(1))                      &
!     *(Eip(3)-Eip(2))*                                                        &
!     (Eip(3)-Eip(4))/(Eip(5)-Eip(1))/(Eip(5)-Eip(2))/(Eip(5)-Eip(3))/         &
!     (Eip(5)-Eip(4))

!           enddo
!        enddo

!C --- Orthogonalize phidot to phi ---?????????????????????????????????
!         do  alpha = 1, 2
!C ...      Same alpha (already orthogonal; do to avoid roundoff error)
!           pr = prodr(0,gn(1,alpha,1),gdot(1,alpha,1),a,b,rofi,nr)
!           gdot(:,alpha,:) = gdot(:,alpha,:) - pr*gn(:,alpha,:)
!           fdot(:,alpha,:) = fdot(:,alpha,:) - pr*fn(:,alpha,:)
!C ...      Different alphas
!           pr = prodr(0,gn(1,3-alpha,1),gdot(1,alpha,1),a,b,rofi,nr)
!           gdot(:,alpha,:) = gdot(:,alpha,:) - pr*gn(:,3-alpha,:)
!           fdot(:,alpha,:) = fdot(:,alpha,:) - pr*fn(:,3-alpha,:)
!        enddo

!        gmtde(:,:) = gdot(:,:,nr)
!        fmtde(:,:) = fdot(:,:,nr)

!        print*,'gmtde = '
!        print*,gmtde

!!!        print*,' '
!!!        print*,'Energies at 5 points for energy
!!!     .  derivative: '
!!!        print*, Eip
!        print*,'Derivatives g(nr) : '
!        print*, gdot(:,:,nr)
!        print*,'Derivatives f(nr) : '
!        print*, fdot(:,:,nr)


!_______________Energy derivative end_______________________________

!__________________Potemtial parameters_____________________________



!        gdotT(1,1) = gdot(1,1,nr); gdotT(2,2) = gdot(2,2,nr) ! Deriv g transposed  nr
!        gdotT(2,1) = gdot(1,2,nr); gdotT(1,2) = gdot(2,1,nr)

!        gnrT(1,1) = gnr(1,1); gnrT(2,2) = gnr(2,2) ! g transposed at nr
!        gnrT(2,1) = gnr(1,2); gnrT(1,2) = gnr(2,1)


!        call mul22(gnr, gdotT, Temp1)            ! Temp1 = g*gdot^T
!        sr = rofi(nr)
!        call matrconst(Temp1,sr,Temp2)           ! Temp2 = g*gdot^T * s
!        call inver(Temp2, Temp1)                 ! Temp1 =(g*gdot^T*s)^-1

!        Anu(:,:) = -Temp1(:,:)

!        print*, Temp1
!        print*, Anu

!        call inver(gnr, Temp1)                   ! Temp1 = g^-1
!        call mul22(fnr, Temp1, Temp2)            ! Temp2 = f(nr)*g^-1
!        call matrconst(Temp2, sr, Temp1)         ! Temp1 = s*f(nr)*g^-1

!        Unit2(1,2) = 0.0_8; Unit2(2,1) = 0.0_8
!        Unit2(1,1) = 1.0_8; Unit2(2,2) = 1.0_8

!        Dnu(:,:) = Temp1(:,:) - kappa(:,:) - Unit2(:,:)

 !       print*, Temp1
 !       print*, (Temp1(:,:)-Unit2(:,:))
!        print*,'D = '
!        print*,Dnu

!        Temp1(:,:) = Dnu(:,:) + Anu(:,:) - l*Unit2(:,:) ! Temp1=D+A-l
!        Temp2(:,:) = Dnu(:,:) + Anu(:,:) + (l+1)*Unit2(:,:) ! T2=D+A+l+1
!        call inver(Temp2,Temp3)                             ! T3=(D+A+l+1)^-1
!        call mul22(Temp1,Temp3,Temp2)  ! T2=(D+A-l)(D+A+l+1)^-1
!        coeff = (sr/avw)**(2*l+1)/real(2*(2*l+1),8)
!        gammal(:,:) = coeff*Temp2

!        call mul22(Anu,Temp3,Temp1)                !T1=A*(D+A+l+1)^-1
!        call mul22(gnrT,Temp1,Temp2)               !T2=g^T*A*(D+A+l+1)^-1
!        coeff = sqrt(avw/2.0_8)*(sr/avw)**(l+1)
!        Wl(:,:) = coeff*Temp2(:,:)

!        call mul22(Anu,gnr,Temp1)                  !T1=A*g
!        call mul22(Temp2,Temp1,Temp3)              !T3=g^T*A*(D+A+l+1)^-1 *A*g
!        call mul22(gnrT,Temp1,Temp2)               !T2=g^T*A*g
!        Cl(:,:) = Enu*Unit2(:,:)+sr*Temp2(:,:)-sr*Temp3(:,:)


!        print*,'Parameters : '
!        print*,'Gamma = '
!        print*,gammal
!        print*,'W, C = ',Wl,Cl


!Co      pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),1:2,1:2)

!        pprel(1,l,imu,1:2,1:2) = Cl(:,:)
!        pprel(2,l,imu,1:2,1:2) = gammal(:,:)
!        pprel(3,l,imu,1:2,1:2) = Dnu(:,:)
!        pprel(4,l,imu,1:2,1:2) = Wl(:,:)

!________________________Calculaing small p________________________

!        gdot1(:,:) = gdot(1,j2,:)
!        gdot2(:,:) = gdot(2,j2,:)

!        do j1=1,2
!           psd(:,1,j1) = gdot(1,j1,:)
!           psd(:,2,j1) = fdot(1,j1,:)
!           psd(:,3,j1) = gdot(2,j1,:)
!           psd(:,4,j1) = fdot(2,j1,:)
!        enddo

!C --- Small parameter p=<gdot|gdot> by trapezium rule ---
!       smalpa = 0.0_8
!       nrmx = nr
!       do  n = 2, nrmx-1
!         np = n + 1 ; r = rofi(n) ; r1 = rofi(np)
!         w0 = (1 + xk *(xk+1) /((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
!         w1 = (1 + xk *(xk+1) /((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
!         w2 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
!         w3 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
!         x1 = a*(r+b) ; x2 = a*(r1+b)
!         do  a1 = 1, 2
!           smalpa(a1,:) = smalpa(a1,:)                                   &
!            + x1*(w0*psd(n ,1,a1)*psd(n ,1,:)+psd(n ,2,a1)*psd(n ,2,:))  &
!            + x2*(w1*psd(np,1,a1)*psd(np,1,:)+psd(np,2,a1)*psd(np,2,:))  &
!            + x1*(w2*psd(n ,3,a1)*psd(n ,3,:)+psd(n ,4,a1)*psd(n ,4,:))  &
!            + x2*(w3*psd(np,3,a1)*psd(np,3,:)+psd(np,4,a1)*psd(np,4,:))
!         enddo
!         smalpa = smalpa/2.0_8
!       enddo

!       pprel(4,l,imu,:,:) = smalpa(:,:)
!!!       print*,'Small p = ',smalpa

!________________________Small p end_______________________________

!      if (imu == 1 .or. imu == 2*l+2) then
!        nug1(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq
!      else
!        nug1(nr) = 1d0 + 2*z/rmax/csq + (enu1 - vm(nr) + up *bm(nr))/csq
!      endif

!      nug2(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq

!      gsmt(1,:) = (c*nug1(nr)*psi(nr,2,:,3) - xk /rmax*psi(nr,1,:,3))
!      gsmt(2,:) = (c*nug2(nr)*psi(nr,4,:,3) - xk1/rmax*psi(nr,3,:,3))

!      if (psi(nr,1,1,3)*psi(nr,3,2,3) < 0.0_8) then
!        psi(:,:,2,:) = - psi(:,:,2,:)
!      endif


!__________________Potential parameters end_________________________

!C --- Orthonormalize Psi_1 and Psi_2
!C      do  ie = 1, 5
!C        pr = prodr(0,psiout(1,1,1,ie),psiout(1,1,2,ie),a,b,rofi,nr)
!C        psiout(:,:,2,ie) = psiout(:,:,2,ie) - pr * psiout(:,:,1,ie)
!CC ...   Normalization of Psi_2 (missing in the previous versions)
!C        norm = dsqrt(prodr(0,psiout(1,1,2,ie),psiout(1,1,2,ie),a,b,rofi,nr))
!C        psiout(:,:,2,ie) = psiout(:,:,2,ie)/norm
!C      enddo
!C
!CC --- Orthonormalize Psi_1 and Psi_2   (Inwards part)
!C      do  ie = 1, 5
!C        pr = prodr(0,psiin(1,1,1,ie),psiin(1,1,2,ie),a,b,rofi,nr)
!C        psiin(:,:,2,ie) = psiin(:,:,2,ie) - pr * psiin(:,:,1,ie)
!CC ...   Normalization of Psi_2 (missing in the previous versions)
!C        norm = dsqrt(prodr(0,psiin(1,1,2,ie),psiin(1,1,2,ie),a,b,rofi,nr))
!C        psiin(:,:,2,ie) = psiin(:,:,2,ie)/norm
!C      enddo

!C --- Write wavefunctions to disk for testing purposes
!c     if (l == 2) then
!c       fil ='wfrel'
!c       call awrit1('%a%i.dat',fil,11,0,imu)
!c       print *,'filename:',fil, ':'
!c       ifi = fopn(fil)
!c       do  n = 1, nrmx
!c         write(ifi,923)rofi(n),((psiout(n,i,alp,3),i=1,4),alp=1,2)
!c       enddo
!c923    format(e12.4,8e15.6)
!c       call fclr(fil,ifi)
!c       print *,'Wrote (l,imu)=',l,imu
!c       call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
!c       read(*,'(a80)') outs
!c       if (outs == 'q') call rx0('quit in rdeq')
!c       if (imu == 2*(l+1)) stop
!c     endif


 9999 continue



!C --- Save radial derivative at rmax ---
!      if (imu == 1 .or. imu == 2*l+2) then
!        nug1(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq
!      else
!        nug1(nr) = 1d0 + 2*z/rmax/csq + (enu1 - vm(nr) + up *bm(nr))/csq
!      endif

!      nug2(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq

!      gsmt(1,:) = (c*nug1(nr)*psiout1(nr,2,:,3) - xk /rmax*psiout1(nr,1,:,3))
!      gsmt(2,:) = (c*nug2(nr)*psiout1(nr,4,:,3) - xk1/rmax*psiout1(nr,3,:,3))

!      if (psiout1(nr,1,1,3)*psiout1(nr,3,2,3) < 0) then
!        psiout1(:,:,2,:) = - psiout1(:,:,2,:)
!      endif

!C --- Save g and f at rmax ---
!      gmt(1,:) = psiout1(nr,1,:,3)/rmax ; gmt(2,:) = psiout1(nr,3,:,3)/rmax
!      fmt(1,:) = psiout1(nr,2,:,3)/rmax ; fmt(2,:) = psiout1(nr,4,:,3)/rmax

!C --- Energy derivatives of rg,rf ---
!      psd(:,:,:) = (-8*(psiout1(:,:,:,2) - psiout1(:,:,:,4)) &
!                      + psiout1(:,:,:,1) - psiout1(:,:,:,5))/12.0_8/de

!C --- Orthogonalize phidot to phi ---
!      do  alpha = 1, 2
!C ...   Same alpha (already orthogonal; do to avoid roundoff error)
!        pr = prodr(0,psiout1(1,1,alpha,3),psd(1,1,alpha),a,b,rofi,nr)
!        psd(:,:,alpha) = psd(:,:,alpha) - pr * psiout1(:,:,alpha,3)
!C ...   Different alphas
!        pr = prodr(0,psiout1(1,1,3-alpha,3),psd(1,1,alpha),a,b,rofi,nr)
!        psd(:,:,alpha) = psd(:,:,alpha) - pr * psiout1(:,:,3-alpha,3)
!      enddo

!C --- Save gdot and fdot at rmax ---
!      gmtde(1,:) = psd(nr,1,:)/rmax; gmtde(2,:) = psd(nr,3,:)/rmax
!      fmtde(1,:) = psd(nr,2,:)/rmax; fmtde(2,:) = psd(nr,4,:)/rmax

!      print*,'gmtde old'
!      print*,gmtde

!C --- Small parameter p=<gdot|gdot> by trapezium rule ---
!      smalpa = 0
!      do  n = 2, nrmx-1
!        np = n + 1 ; r = rofi(n) ; r1 = rofi(np)
!        w0 = (1 + xk *(xk+1) /((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
!        w1 = (1 + xk *(xk+1) /((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
!        w2 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
!        w3 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
!        x1 = a*(r+b) ; x2 = a*(r1+b)
!        do  a1 = 1, 2
!          smalpa(a1,:) = smalpa(a1,:)
!     .      + x1*(w0*psd(n ,1,a1)*psd(n ,1,:)+psd(n ,2,a1)*psd(n ,2,:))
!     .      + x2*(w1*psd(np,1,a1)*psd(np,1,:)+psd(np,2,a1)*psd(np,2,:))
!     .      + x1*(w2*psd(n ,3,a1)*psd(n ,3,:)+psd(n ,4,a1)*psd(n ,4,:))
!     .      + x2*(w3*psd(np,3,a1)*psd(np,3,:)+psd(np,4,a1)*psd(np,4,:))
!        enddo
!        smalpa = smalpa/2
!      enddo

!      pprel(4,l,imu,:,:) = smalpa(:,:)

!C --- Overlap integrals for the calculation of averages [Turek (6.141)]
!C     Should be recorded for future use along with those for different mu's
!C     (for now just print out)
!C ... psi(n,kappa,alpha,ie)

!c     write(6,900)l,mu
!c900  format('Overlap parameters R for l=',i1,' , mu= ',F4.1)
!c     do a1 = 1, 2
!c       do a2 = 1, 2
!c         do k1 = 1, 4
!c           do k2 = 1, 4
!c             if (mod(k1,2) /= mod(k2,2)) cycle ! don't need <g|f>
!c             do ie = 1, 5
!c              pr = prodr(1,psiout(1,k1,a1,ie),psiout(1,k2,a2,ie),a,b,rofi,nr)
!c              if (mod(k1,2) == 1) then
!c                write(6,901)ie,k1/2+1,a1,k2/2+1,a2,pr
!c              else
!c                write(6,902)ie,k1/2,a1,k2/2,a2,pr
!c              endif
!c             enddo
!c           enddo
!c         enddo
!c       enddo
!c     enddo
!c901  format(1x,'ie=',i1,' Rgg(',i1,',',i1,';',i1,',',i1,')=',F11.8)
!c902  format(1x,'ie=',i1,' Rff(',i1,',',i1,';',i1,',',i1,')=',F11.8)




!      print*,' '
!      print*,'exiting rdeq'
!      print*,'Pv = ', Pv
!      print*,'enum1 = ',enum1,'; enum2 = ',enum2,
!     .'; ksop (1,1,1,4) = ', ksop(1,1,1,4),
!     . '; z = ',z,'; v (100,1) = ',v(100,1),' ; rofi(100) = ',
!     .  rofi(100),
!     . '; nr = ',nr,' ; nsp = ',nsp,' ; a = ',a,' ; b = ',b,
!     .' ; l = ',l,' ; lmx = ',lmx,' ; imu = ',imu,
!     . ' ; scaldele = ',scalede,' ; gmt = ',gmt,
!     . ' ; fmt = ',fmt,' ; gmtde = ',gmtde,
!     . ' ; fmtde = ',fmtde,' ; gsmt = ',gsmt,
!     . ' ; pprel(1,1,1,2,2) = ',pprel(1,1,1,2,2)
!      close(1)

      end subroutine rdeqcore


!AV----------------------- Integration inward-outward_______________

      subroutine integrationC(alpha, En, xk, xk1, mu, up, u,           &
          upp, sqru, vm,                                               &
          bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,     &
          g1out,f1out,g2out,f2out,g1sout,f1sout,g2sout,f2sout,         &
          g1in,f1in,g2in,f2in,g1sin,f1sin,g2sin,f2sin,                 &
          nmrseq,nnod,normsave)
      implicit none
!C ... Passed parameters
      integer nr,nsp,imu,l,lmx,ie, jalp,j1,j2,nfr,ncsave,nnod, ncut,nmrseq
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,9)
      double precision  grseq(nr,2)
!C ... Local parameters
      integer n,np,alpha,nrk,it,nrmx, n0, nprint, nit, nm, nc,deltacut,&
              ncut1, ncut2
!C     integer k1,k2,i,alp,ifi,fopn, ie
!C     character*(10) outs*80, fil*12
      double precision J, &   !Jacobian
        c,csq,                                                        & !Speed of light 274.072d0 !     .  de,                  !Finite energy difference for num diff.
       bmint,df11,df12,df13,df14,df21,df22,df23,df24,dg11,dg12,dg13,  &
       dg14,dg21,dg22,dg23,dg24,dx,e1,e2,En,f1c,f1p,f1sn,             &
       f1snf1,f1snm1,f1snm2,f1snm3,f1snm4,f2c,f2p,f2sn,f2snf2,f2snm1, &
       f2snm2,f2snm3,f2snm4,ff,g1c,g1p,g1sn,g1sng1,g1snm1,g1snm2,     &
       g1snm3,g1snm4,g2c,g2p,g2sn,g2sng2,g2snm1,g2snm2,g2snm3,g2snm4, &
       gamma,hoc,mu,norm,nuf1int,nuf2int,nug1int,nug2int,pr,prodr,q,r,&
       r1,rmax,smalpa(2,2),sqru,u,up,upp,vmint,w0,w1,w2,w3,x1,x2,xk,  &
       xk1, eps1, eps2, eps3, eps4, Delt, r0, mu1, q01, mu2, q02,rm

      double precision g1out(nr),g2out(nr),f1out(nr),f2out(nr),   &
       g1sout(nr),g2sout(nr),f1sout(nr),f2sout(nr),nug1(nr),      &
      nug2(nr),nuf1(nr),nuf2(nr),bm(nr),vm(nr),psiout(nr,4,2,5),  &
      psiin(nr,4,2,5),psd(nr,4,2), g1in(nr),g2in(nr),             &
       f1in(nr),f2in(nr),g1sin(nr),g2sin(nr),                     &
       f1sin(nr),f2sin(nr), psiinn(nr,4), psioutn(nr,4)
      double precision fr, dg1new, dg1old, g01,g02,normn,norm1,   &
                       norm2,normni

      double precision lamb1, lamb2, sq11, sq12, sq21, sq22, khi, &
                       rin, sigma1, sigma2, b11, b12, a11, a12,   &
                       g01i, g02i, f01i, f02i,la1,la2,si1,si2
      double precision g1inmax, g1outmax,normsave,frac,nodsig

!C ... External calls
!      external dinv22,mul22,rx
!C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

!______________________Zeros_________

         g1out(:)=0d0; f1out(:)=0d0; g2out(:)=0d0; f2out(:)=0d0
         g1sout(:)=0d0; f1sout(:)=0d0; g2sout(:)=0d0; f2sout(:)=0d0
         g1in(:)=0d0; f1in(:)=0d0; g2in(:)=0d0; f2in(:)=0d0
         g1sin(:)=0d0; f1sin(:)=0d0; g2sin(:)=0d0; f2sin(:)=0d0

!______________________EndZeroz______


      csq = c*c
      hoc = (2*z/c)
      dx = 1d0
      ff = 1.0_8
      rmax = rofi(nr)

!C --- Integration inwards part ---
      nprint = 0
      n0 = 30
      r0 = rofi(n0)
      nnod = 0

      e2 =  En ! -(l+1) * ksop(l,1,1,1)/2.0_8
      e1 =  En !+ l    * ksop(l,1,1,1)/2.0_8

!      e1 = En
!      e2 = En


      if (imu == 1 .or. imu == 2*l+2) then
        xk = -l-1
      else
        xk = l
      endif
      xk1 = -l-1

      if (imu == 1 .or. imu == 2*l+2) then
        nug1(2:nr) = 1d0 + 2*z/rofi(2:nr)/csq + ff*(e2 - vm(2:nr) - upp*bm(2:nr))/csq
        nuf1(2:nr) = - (e2 + 2*z/rofi(2:nr) - vm(2:nr) - u  *bm(2:nr))
      else
        nug1(2:nr) = 1d0 + 2*z/rofi(2:nr)/csq + ff*(e1 - vm(2:nr) + up*bm(2:nr))/csq
        nuf1(2:nr) = - (e1 + 2*z/rofi(2:nr) - vm(2:nr) + u *bm(2:nr))
      endif

      nug2(2:nr) = 1d0 + 2*z/rofi(2:nr)/csq + ff*(e2 - vm(2:nr) - upp*bm(2:nr))/csq
      nuf2(2:nr) = - (e2 + 2*z/rofi(2:nr) - vm(2:nr) - u  *bm(2:nr))

!______________Other initial conditions__________________

!        lamb1 =  dsqrt(dabs(-e1-e1**2/csq))
!        lamb2 =  dsqrt(dabs(-e2-e2**2/csq))
!        sq11 = dsqrt(dabs((csq+e1**2)/2/csq))
!        sq12 = dsqrt(dabs((csq+e2**2)/2/csq))
!        sq21 = dsqrt(dabs((csq-e1**2)/2/csq))
!        sq22 = dsqrt(dabs((csq-e2**2)/2/csq))
!        khi = 1
!        rin = rofi(nr)
!        sigma1 = e1*khi/csq/lamb1
!        sigma2 = e2*khi/csq/lamb2
!        if (alpha==1) then
!           b11 = (xk+khi/lamb1)/2/c
!           a11 = (xk+(1-sigma1)*e1/csq-khi*lamb1/csq)*c/lamb1
!           b12 = (xk+khi/lamb2)/2/c
!           a12 = (xk+(1-sigma2)*e2/csq-khi*lamb2/csq)*c/lamb2
!        else
!           b11 = (xk1+khi/lamb1)/2/c
!           a11 = (xk1+(1-sigma1)*e1/csq-khi*lamb1/csq)*c/lamb1
!           b12 = (xk1+khi/lamb2)/2/c
!           a12 = (xk1+(1-sigma2)*e2/csq-khi*lamb2/csq)*c/lamb2
!        endif
!        g01i = rin**sigma1*exp(-lamb1*rin)*(sq11*(1+a11/rin)+  &
!             sq21*(b11/rin))
!        g02i = rin**sigma2*exp(-lamb2*rin)*(sq12*(1+a12/rin)+  &
!             sq22*(b12/rin))
!        f01i = rin**sigma1*exp(-lamb1*rin)*(sq11*(1+a11/rin)+  &
!             sq21*(b11/rin))
!        f02i = rin**sigma2*exp(-lamb2*rin)*(sq12*(1+a12/rin)+  &
!             sq22*(b12/rin))
!        print*,'g01, f01, g02, f02: ',g01i, f01i, g02i, f02i

!        g1outmax = 0
!        g1inmax = 0


!______________End other initial conditions_______________



      nrk = 6 ; nrmx = nr
      nrmx = nr
      g1out(1) = 0.0_8 ; f1out(1) = 0.0_8
      g2out(1) = 0.0_8 ; f2out(1) = 0.0_8

!C --- Initial conditions at r-->0: V = V0 - 2Z/r, B = const
!C     The wave function is a polynomial g0i = r^(gi-1)*Sum_over_v(g1v*r^v)

      if (alpha == 1) then
        gamma = dsqrt(xk*xk - hoc*hoc)
        g1out(2) = grseq(2,1)!1d0
        f1out(2) = grseq(2,2)!(xk + gamma)/hoc
        g2out(2) = 0.0_8
        f2out(2) = 0.0_8
      else
        gamma = dsqrt(xk1*xk1 - hoc*hoc)
        g1out(2) = 0.0_8
        f1out(2) = 0.0_8
        g2out(2) = grseq(2,1)!1d0
        f2out(2) = grseq(2,2)!(xk1 + gamma)/hoc
      endif

!      e1 = En;
!      e2 = En;
 !     e2 = enu2 + (ie - 3)*de ; e1 = enu1 + (ie - 3)*de



!C     nug1(1) = nug1(2); nug2(1) = nug2(2)
!C     nuf1(1) = nuf1(2); nuf2(1) = nuf2(2)
!C     Extrapolate to the origin
!      call ratint(rofi(2),nug1(2),4,rofi(1),nug1(1),w0)
!      call ratint(rofi(2),nuf1(2),4,rofi(1),nuf1(1),w0)
!      call ratint(rofi(2),nug2(2),4,rofi(1),nug2(1),w0)
!      call ratint(rofi(2),nuf2(2),4,rofi(1),nuf2(1),w0)

!C --- Runge-Kutta for the first NKR points ------------------------------
      J = a*(rofi(2) + b) ; q = dexp(a/2) ; r = rofi(2)

      g1sout(2) = J*(c*nug1(2)*f1out(2) - xk /r*g1out(2))
      f1sout(2) = J*( (nuf1(2)*g1out(2) - sqru*bm(2)*g2out(2))/c + xk /r*f1out(2))
      g2sout(2) = J*(c*nug2(2)*f2out(2) - xk1/r*g2out(2))
      f2sout(2) = J*( (nuf2(2)*g2out(2) - sqru*bm(2)*g1out(2))/c + xk1/r*f2out(2))

!C     First point
      n = 2

!C     Re-entry for subsequent points
    2 continue
      g1c = g1out(n) ; f1c = f1out(n) ; g2c = g2out(n) ; f2c = f2out(n)
      dg11 = dx*J* (c*nug1(n)*f1c                     - xk /r*g1c)
      df11 = dx*J* ( (nuf1(n)*g1c - sqru*bm(n)*g2c)/c + xk /r*f1c)
      dg21 = dx*J* (c*nug2(n)*f2c                     - xk1/r*g2c)
      df21 = dx*J* ( (nuf2(n)*g2c - sqru*bm(n)*g1c)/c + xk1/r*f2c)

      g1c = g1out(n) + dg11/2.0_8 ; f1c = f1out(n) + df11/2.0_8
      g2c = g2out(n) + dg21/2.0_8 ; f2c = f2out(n) + df21/2.0_8

      J = J*q ; r = J/a - b

      vmint = (5*vm(n) + 2*vm(n+1) + vm(n+2))/8.0_8
      bmint = (5*bm(n) + 2*bm(n+1) + bm(n+2))/8.0_8

!c     vmint = (6*vm(n) + 5*vm(n+1) - 3*vm(n-1))/8
!c     bmint = (6*bm(n) + 5*bm(n+1) - 3*vm(n-1))/8

      if (imu == 1 .or. imu == 2*l+2) then
        nug1int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp * bmint)/csq
        nuf1int =     - 2*z/r     - (e2 - vmint - u   * bmint)
      else
        nug1int = 1d0 + 2*z/r/csq + ff*(e1 - vmint + up * bmint)/csq
        nuf1int =     - 2*z/r     - (e1 - vmint + u  * bmint)
      endif

      nug2int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp * bmint)/csq
      nuf2int =     - 2*z/r     - (e2 - vmint - u   * bmint)


      dg12 = dx*J*(c*nug1int*f1c                     - xk /r*g1c)
      df12 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + xk /r*f1c)
      dg22 = dx*J*(c*nug2int*f2c                     - xk1/r*g2c)
      df22 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + xk1/r*f2c)

      g1c = g1out(n) + dg12/2.0_8 ; f1c = f1out(n) + df12/2.0_8
      g2c = g2out(n) + dg22/2.0_8 ; f2c = f2out(n) + df22/2.0_8

      dg13 = dx*J*(c*nug1int*f1c                     - xk /r*g1c)
      df13 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + xk /r*f1c)
      dg23 = dx*J*(c*nug2int*f2c                     - xk1/r*g2c)
      df23 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + xk1/r*f2c)

      g1c = g1out(n) + dg13 ; f1c = f1out(n) + df13
      g2c = g2out(n) + dg23 ; f2c = f2out(n) + df23

      n = n + 1

      J = J*q ; r = J/a - b

      dg14 = dx*J*(c*nug1(n)*f1c                     - xk /r*g1c)
      df14 = dx*J*( (nuf1(n)*g1c - sqru*bm(n)*g2c)/c + xk /r*f1c)
      dg24 = dx*J*(c*nug2(n)*f2c                     - xk1/r*g2c)
      df24 = dx*J*( (nuf2(n)*g2c - sqru*bm(n)*g1c)/c + xk1/r*f2c)

      g1out(n) = g1out(n-1) + (dg11 + 2*(dg12+dg13) + dg14)/6.0_8
      f1out(n) = f1out(n-1) + (df11 + 2*(df12+df13) + df14)/6.0_8
      g2out(n) = g2out(n-1) + (dg21 + 2*(dg22+dg23) + dg24)/6.0_8
      f2out(n) = f2out(n-1) + (df21 + 2*(df22+df23) + df24)/6.0_8

!C ... Derivatives dP/dx and dQ/dx----------------------------

      g1sout(n) = J*(c*nug1(n)*f1out(n) - xk /r*g1out(n))
      f1sout(n) = J*( (nuf1(n)*g1out(n) - sqru*bm(n)*g2out(n))/c + xk /r*f1out(n))
      g2sout(n) = J*(c*nug2(n)*f2out(n) - xk1/r*g2out(n))
      f2sout(n) = J*( (nuf2(n)*g2out(n) - sqru*bm(n)*g1out(n))/c + xk1/r*f2out(n))

!C --- The rest of the integration by Milne's method

      if (n < nrk) goto 2
      if (n == nrmx) goto 12

      g1sn = g1sout(n) ; f1sn = f1sout(n)
      g2sn = g2sout(n) ; f2sn = f2sout(n)

      g1snm1 = g1sout(n-1) ; f1snm1 = f1sout(n-1)
      g2snm1 = g2sout(n-1) ; f2snm1 = f2sout(n-1)

      g1snm2 = g1sout(n-2) ; f1snm2 = f1sout(n-2)
      g2snm2 = g2sout(n-2) ; f2snm2 = f2sout(n-2)

      g1snm3 = g1sout(n-3) ; f1snm3 = f1sout(n-3)
      g2snm3 = g2sout(n-3) ; f2snm3 = f2sout(n-3)

      g1snm4 = g1sout(n-4) ; f1snm4 = f1sout(n-4)
      g2snm4 = g2sout(n-4) ; f2snm4 = f2sout(n-4)

!C   4 J = J*q*q ; r = J/a - b
    4 continue
      J = a*(rofi(n+1) + b); r = rofi(n+1)

      g1p = g1out(n-5) + (3*dx/10.0_8)*(11*g1sn-14*g1snm1+26*g1snm2-14*g1snm3+11*g1snm4)
      f1p = f1out(n-5) + (3*dx/10.0_8)*(11*f1sn-14*f1snm1+26*f1snm2-14*f1snm3+11*f1snm4)
      g2p = g2out(n-5) + (3*dx/10.0_8)*(11*g2sn-14*g2snm1+26*g2snm2-14*g2snm3+11*g2snm4)
      f2p = f2out(n-5) + (3*dx/10.0_8)*(11*f2sn-14*f2snm1+26*f2snm2-14*f2snm3+11*f2snm4)

      do  it = 1, 100
        g1sng1 = J*(c*nug1(n+1)*f1p                       - xk/r*g1p)
        f1snf1 = J*( (nuf1(n+1)*g1p - sqru*bm(n+1)*g2p)/c + xk/r*f1p)
        g2sng2 = J*(c*nug2(n+1)*f2p                       - xk1/r*g2p)
        f2snf2 = J*( (nuf2(n+1)*g2p - sqru*bm(n+1)*g1p)/c + xk1/r*f2p)

        g1c = g1out(n-3) + (2*dx/45.0_8)*(7*g1sng1+32*g1sn+12*g1snm1+32*g1snm2+7*g1snm3)
        f1c = f1out(n-3) + (2*dx/45.0_8)*(7*f1snf1+32*f1sn+12*f1snm1+32*f1snm2+7*f1snm3)
        g2c = g2out(n-3) + (2*dx/45.0_8)*(7*g2sng2+32*g2sn+12*g2snm1+32*g2snm2+7*g2snm3)
        f2c = f2out(n-3) + (2*dx/45.0_8)*(7*f2snf2+32*f2sn+12*f2snm1+32*f2snm2+7*f2snm3)

        if (dabs(g1c-g1p) <= dabs(g1c)*1d-12 .and.  &
            dabs(f1c-f1p) <= dabs(f1c)*1d-12 .and.  &
            dabs(g2c-g2p) <= dabs(g2c)*1d-12 .and.  &
            dabs(f2c-f2p) <= dabs(f2c)*1d-12) exit

!C ...   Check for convergence
        if (it == 100) then
          print *,abs(g1c-g1p),abs(f1c-f1p),abs(g2c-g2p),abs(f2c-f2p)
          call rx('Convergence not reached after 100 iterations')
        endif
        g1p = g1c ; f1p = f1c ; g2p = g2c ; f2p = f2c
      enddo

      n = n + 1
      if ((g1outmax <= dabs(g1c)).and.(n<1460)) g1outmax = dabs(g1c)

      g1out(n) = g1c ; f1out(n) = f1c
      g2out(n) = g2c ; f2out(n) = f2c
      if (alpha==1) then
           nodsig = g1out(n)*g1out(n-1)
      else
           nodsig = g2out(n)*g2out(n-1)
      endif
      if ((n <= nmrseq).and.(nodsig <= 0)) nnod = nnod + 1
      g1sout(n) = g1sng1; f1sout(n) = f1snf1
      g2sout(n) = g2sng2; f2sout(n) = f2snf2
      g1snm4 = g1snm3; f1snm4 = f1snm3
      g2snm4 = g2snm3; f2snm4 = f2snm3
      g1snm3 = g1snm2; f1snm3 = f1snm2
      g2snm3 = g2snm2; f2snm3 = f2snm2
      g1snm2 = g1snm1; f1snm2 = f1snm1
      g2snm2 = g2snm1; f2snm2 = f2snm1
      g1snm1 = g1sn ; f1snm1 = f1sn
      g2snm1 = g2sn ; f2snm1 = f2sn
      g1sn = g1sng1 ; f1sn = f1snf1
      g2sn = g2sng2 ; f2sn = f2snf2

      if (n < nrmx) goto 4

   12 continue

!___________________Remnant of the initial rdeq_________________
!C ... Save and normalize the Psi functions
!      do jalp=1,2
!      do ie=1,5
!        psiout(:,1,jalp,ie) = g1out(:)
!        psiout(:,2,jalp,ie) = f1out(:)
!        psiout(:,3,jalp,ie) = g2out(:)
!        psiout(:,4,jalp,ie) = f2out(:)
!      enddo
!      enddo

!___________________End of Remnant of the initial rdeq__________

!___________________For normalization___________________________

      psioutn(:,1)=g1out(:)
      psioutn(:,2)=f1out(:)
      psioutn(:,3)=g2out(:)
      psioutn(:,4)=f2out(:)

!__________________End _For normalization________________________
!C-------  Initial conditions (Ebert)---------------------

!        mu1 = dsqrt(dabs(-En-En**2/csq))
!        q01 = -mu1/(1d0+En/csq)/c
!        mu2 = dsqrt(dabs(-En-En**2/csq))
!        q02 = -mu2/(1d0+En/csq)/c

!        g01 = 1d0 !1e-3 !-1e-15    !e-3
!        g02 = g01
!        g02 = 1e-2
!        if ((grseq(nr,1)==0).or.(grseq(nr,2)==0)) then
!          g01 = g1out(nmrseq)!*10
!          g02 = g2out(nmrseq)!*10
!          mu1 = dsqrt(dabs(-e2-e2**2/csq))
!          q01 = -g01*mu1/(1d0+e2/csq)/c
!          mu2 = dsqrt(dabs(-e1-e1**2/csq))
!         q02 = -g02*mu2/(1d0+e1/csq)/c
!          la1=dsqrt(2*dabs(e1))
!          la2=dsqrt(2*dabs(e2))
!          si1=z/la1
!          si2=z/la2
!          r = rofi(nr)
!          if (alpha == 1) then
!            g1in(nr) =  g01*dexp(-mu1*r)    !r**si1*dexp(-la1*r)
!            f1in(nr) =  q01*dexp(-mu1*r)    !-la1*r**si1*dexp(-la1*r)
!            g2in(nr) = 0d0
!            f2in(nr) = 0d0
!          else
!            g1in(nr) = 0d0
!            f1in(nr) = 0d0
!            g2in(nr) =  g02*dexp(-mu2*r)    !r**si2*dexp(-la2*r)
!            f2in(nr) =  q02*dexp(-mu2*r)    !-la2*r**si2*dexp(-la2*r)
!          endif
!        else
!          if (imu == 1 .or. imu == 2*l+2) then
!            if (alpha == 1) then
!              g1in(nr) =   grseq(nr,1)*g1out(nmrseq)/grseq(nmrseq,1)
!              f1in(nr) =   grseq(nr,2)*f1out(nmrseq)/grseq(nmrseq,2)
!              g2in(nr) = 0d0
!              f2in(nr) = 0d0
!            else
!              g1in(nr) = 0d0
!              f1in(nr) = 0d0
!              g2in(nr) =   grseq(nr,1)*g2out(nmrseq)/grseq(nmrseq,1)
!              f2in(nr) =   grseq(nr,2)*f2out(nmrseq)/grseq(nmrseq,2)
!               endif
!          else
!            g01 = g1out(nmrseq)*10
!            g02 = g2out(nmrseq)*10
!            mu1 = dsqrt(dabs(-e1-e1**2/csq))
!            q01 = -g01*mu1/(1d0+e2/csq)/c
!            mu2 = dsqrt(dabs(-e2-e2**2/csq))
!            q02 = -g02*mu2/(1d0+e1/csq)/c
!            r = rofi(nr)
!            if (alpha == 1) then
!              g1in(nr) = g01*dexp(-mu1*r)
!              f1in(nr) = q01*dexp(-mu1*r)
!              g2in(nr) = 0d0
!              f2in(nr) = 0d0
!            else
!              g1in(nr) = 0d0
!              f1in(nr) = 0d0
!              g2in(nr) = g02*dexp(-mu2*r)
!              f2in(nr) = q02*dexp(-mu2*r)
!               endif
!          endif
!        endif

!        print*,' Init:',g1in(nr),f1in(nr)
!        print*,'mu*r',(mu1*rofi(nr))
            if (alpha == 1) then
              g1in(nr) = grseq(nr,1)
              f1in(nr) = grseq(nr,2)
              g2in(nr) = 0d0
              f2in(nr) = 0d0
            else
              g1in(nr) = 0d0
              f1in(nr) = 0d0
              g2in(nr) = grseq(nr,1)
              f2in(nr) = grseq(nr,2)
            endif

!!          if (imu == 1 .or. imu == 2*l+2) then
!!            if (alpha == 1) then
!!              g1in(nr) = grseq(nr,1)
!!              f1in(nr) = grseq(nr,2)
!!              g2in(nr) = 0d0
!!              f2in(nr) = 0d0
!!            else
!!              g1in(nr) = 0d0
!!              f1in(nr) = 0d0
!!              g2in(nr) = grseq(nr,1)
!!              f2in(nr) = grseq(nr,2)
!!            endif
!!          else
!!             g01 = g1out(nmrseq)!*10
!!             g02 = g2out(nmrseq)!*10
!!             mu1 = dsqrt(dabs(-e2-e2**2/csq))
!!             q01 = -g01*mu1/(1d0+e2/csq)/c
!!             mu2 = dsqrt(dabs(-e1-e1**2/csq))
!!             q02 = -g02*mu2/(1d0+e1/csq)/c
!!             r = rofi(nr)
!!             if (alpha == 1) then
!!                g1in(nr) =  g01*dexp(-mu1*r)    !r**si1*dexp(-la1*r)
!!                f1in(nr) =  q01*dexp(-mu1*r)    !-la1*r**si1*dexp(-la1*r)
!!                g2in(nr) = 0d0
!!                f2in(nr) = 0d0
!!             else
!!                g1in(nr) = 0d0
!!                f1in(nr) = 0d0
!!                g2in(nr) =  g02*dexp(-mu2*r)    !r**si2*dexp(-la2*r)
!!                f2in(nr) =  q02*dexp(-mu2*r)    !-la2*r**si2*dexp(-la2*r)
!!             endif
!!          endif


        J = a*(rofi(nr) + b) ; q = dexp(-a/2.0_8) ; r = rofi(nr)


        g1sin(nr) = J*(c*nug1(nr)*f1in(nr) - xk/r*g1in(nr))
        f1sin(nr) = J*( (nuf1(nr)*g1in(nr) - sqru*bm(nr)*g2in(nr))/c + xk/r*f1in(nr))
        g2sin(nr) = J*(c*nug2(nr)*f2in(nr) - xk1/r*g2in(nr))
        f2sin(nr) = J*( (nuf2(nr)*g2in(nr) - sqru*bm(nr)*g1in(nr))/c + xk1/r*f2in(nr))

!C   --- Runge-Kutta for the first NKR points ---
        do  n = nr, nr-5, -1
          g1c = g1in(n)
          f1c = f1in(n)
          g2c = g2in(n)
          f2c = f2in(n)

          dg11 = dx*J*(c*nug1(n)*f1c - xk/r*g1c)
          df11 = dx*J*((nuf1(n)*g1c - sqru*bm(n)*g2c)/c + xk/r*f1c)
          dg21 = dx*J*(c*nug2(n)*f2c - xk1/r*g2c)
          df21 = dx*J*((nuf2(n)*g2c - sqru*bm(n)*g1c)/c + xk1/r*f2c)

          g1c = g1in(n) - dg11/2.0_8
          f1c = f1in(n) - df11/2.0_8
          g2c = g2in(n) - dg21/2.0_8
          f2c = f2in(n) - df21/2.0_8

          J = J*q
          r = J/a - b

          vmint = (5*vm(n) + 2*vm(n-1) + vm(n-2))/8.0_8
          bmint = (5*bm(n) + 2*bm(n-1) + bm(n-2))/8.0_8

!c     vmint = (6*vm(n) + 5*vm(n-1) - 3*vm(n+1))/8
!c     bmint = (6*bm(n) + 5*bm(n-1) - 3*vm(n+1))/8

          if (imu == 1 .or. imu == 2*l+2) then
            nug1int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp*bmint)/csq
            nuf1int = -2*z/r - (e2 - vmint - u*bmint)
          else
            nug1int = 1d0 + 2*z/r/csq + ff*(e1 - vmint + up*bmint)/csq
            nuf1int = -2*z/r - (e1 - vmint + u*bmint)
          endif

          nug2int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp*bmint)/csq
          nuf2int = -2*z/r - (e2 - vmint - u*bmint)

          dg12 = dx*J*(c*nug1int*f1c - xk/r*g1c)
          df12 = dx*J*((nuf1int*g1c - sqru*bmint*g2c)/c + xk/r*f1c)
          dg22 = dx*J*(c*nug2int*f2c - xk1/r*g2c)
          df22 = dx*J*((nuf2int*g2c - sqru*bmint*g1c)/c + xk1/r*f2c)

          g1c = g1in(n) - dg12/2.0_8
          f1c = f1in(n) - df12/2.0_8
          g2c = g2in(n) - dg22/2.0_8
          f2c = f2in(n) - df22/2.0_8

          dg13 = dx*J*(c*nug1int*f1c - xk/r*g1c)
          df13 = dx*J*((nuf1int*g1c - sqru*bmint*g2c)/c + xk/r*f1c)
          dg23 = dx*J*(c*nug2int*f2c - xk1/r*g2c)
          df23 = dx*J*((nuf2int*g2c - sqru*bmint*g1c)/c + xk1/r*f2c)

          g1c = g1in(n) - dg13
          f1c = f1in(n) - df13
          g2c = g2in(n) - dg23
          f2c = f2in(n) - df23

          J = J*q
          r = J/a - b

          dg14 = dx*J*(c*nug1(n)*f1c - xk/r*g1c)
          df14 = dx*J*((nuf1(n)*g1c - sqru*bm(n)*g2c)/c + xk/r*f1c)
          dg24 = dx*J*(c*nug2(n)*f2c - xk1/r*g2c)
          df24 = dx*J*((nuf2(n)*g2c - sqru*bm(n)*g1c)/c + xk1/r*f2c)

          g1in(n-1) = g1in(n) - (dg11 + 2*(dg12+dg13) + dg14)/6.0_8
          f1in(n-1) = f1in(n) - (df11 + 2*(df12+df13) + df14)/6.0_8
          g2in(n-1) = g2in(n) - (dg21 + 2*(dg22+dg23) + dg24)/6.0_8
          f2in(n-1) = f2in(n) - (df21 + 2*(df22+df23) + df24)/6.0_8

!C --- Determine the derivatives dP/dx and dQ/dx ---

          g1sin(n-1) = J*(c*nug1(n-1)*f1in(n-1) - xk/r*g1in(n-1))
          f1sin(n-1) = J*((nuf1(n-1)*g1in(n-1) - sqru*bm(n-1)*g2in(n-1))/c +  xk/r*f1in(n-1))
          g2sin(n-1) = J*(c*nug2(n-1)*f2in(n-1) - xk1/r*g2in(n-1))
          f2sin(n-1) = J*((nuf2(n-1)*g2in(n-1) - sqru*bm(n-1)*g1in(n-1))/c +  xk1/r*f2in(n-1))
        enddo

!C---------------------------------------------------------------------
!C     The rest of the integration is performed using Milne's method.
!C---------------------------------------------------------------------
        n = nr-5

        g1sn = g1sin(n)
        f1sn = f1sin(n)
        g1snm1 = g1sin(n+1)
        f1snm1 = f1sin(n+1)
        g1snm2 = g1sin(n+2)
        f1snm2 = f1sin(n+2)
        g1snm3 = g1sin(n+3)
        f1snm3 = f1sin(n+3)
        g1snm4 = g1sin(n+4)
        f1snm4 = f1sin(n+4)

        g2sn = g2sin(n)
        f2sn = f2sin(n)
        g2snm1 = g2sin(n+1)
        f2snm1 = f2sin(n+1)
        g2snm2 = g2sin(n+2)
        f2snm2 = f2sin(n+2)
        g2snm3 = g2sin(n+3)
        f2snm3 = f2sin(n+3)
        g2snm4 = g2sin(n+4)
        f2snm4 = f2sin(n+4)

!        fr=1  ! Random but must be postive. For the first check in the following loop
!        n = nr-5
!        nfr = 0
!        if (alpha==1) nnod = 0

!________________alpha=1____________________________________-
        if (alpha == 1) then

        do n = nr-5, 2, -1
!        do while (n>1400)

!          if (l==2) print*,'n, nr, fr: ',n, nr,fr

          J = a*(rofi(n) + b); r = rofi(n)
!C          J = J*q*q
!C          r = J/a - b

          g1p = g1in(n+5) - (3*dx/10.0_8)*(11*g1sn - 14*g1snm1 + 26*g1snm2 - 14*g1snm3 + 11*g1snm4)
          f1p = f1in(n+5) - (3*dx/10.0_8)*(11*f1sn - 14*f1snm1 + 26*f1snm2 - 14*f1snm3 + 11*f1snm4)
          g2p = g2in(n+5) - (3*dx/10.0_8)*(11*g2sn - 14*g2snm1 + 26*g2snm2 - 14*g2snm3 + 11*g2snm4)
          f2p = f2in(n+5) - (3*dx/10.0_8)*(11*f2sn - 14*f2snm1 + 26*f2snm2 - 14*f2snm3 + 11*f2snm4)

          eps1 = dabs(g1p-g1c)
          eps2 = dabs(f1p-f1c)
          eps3 = dabs(g2p-g2c)
          eps4 = dabs(f2p-f2c)

          nit = 0
          Delt = 1d-12
          do while (eps1>dabs(g1c)*Delt .or. eps2>dabs(f1c)*Delt &
               .or. eps3>dabs(g2c)*Delt .or. eps4>dabs(f2c)*Delt)

             g1sng1 = J*(c*nug1(n-1)*f1p - xk/r*g1p)
             f1snf1 = J*((nuf1(n-1)*g1p - sqru*bm(n-1)*g2p)/c + xk/r*f1p)
             g2sng2 = J*(c*nug2(n-1)*f2p - xk1/r*g2p)
             f2snf2 = J*((nuf2(n-1)*g2p - sqru*bm(n-1)*g1p)/c + xk1/r*f2p)

             g1c = g1in(n+3) - (2*dx/45.0_8)*(7*g1sng1 + 32*g1sn + 12*g1snm1 + 32*g1snm2 + 7*g1snm3)
             f1c = f1in(n+3) - (2*dx/45.0_8)*(7*f1snf1 + 32*f1sn + 12*f1snm1 + 32*f1snm2 + 7*f1snm3)
             g2c = g2in(n+3) - (2*dx/45.0_8)*(7*g2sng2 + 32*g2sn + 12*g2snm1 + 32*g2snm2 + 7*g2snm3)
             f2c = f2in(n+3) - (2*dx/45.0_8)*(7*f2snf2 + 32*f2sn + 12*f2snm1 + 32*f2snm2 + 7*f2snm3)

             eps1 = dabs(g1c-g1p)
             eps2 = dabs(f1c-f1p)
             eps3 = dabs(g2c-g2p)
             eps4 = dabs(f2c-f2p)

!C...      Check for convergence
!          if (nit == 100) then
!            call info5(1,0,0,'rdq2 (warning): not converged for '//
!     .        'alpha=%i  ie=%i n=%i  kappa=%d  mu=%d',alpha,ie,n,xk,mu)
!C           nprint = nprint + 1
!            exit
!          endif

            g1p = g1c
            f1p = f1c
            g2p = g2c
            f2p = f2c

            nit = nit + 1
        enddo

!        if ((g1inmax <= dabs(g1c)).and.(n>1300)) g1inmax = dabs(g1c)
        g1in(n-1) = g1c
        f1in(n-1) = f1c
        nodsig = g1in(n-1)*g1in(n)
        if ((n >= nmrseq).and.(nodsig <= 0)) nnod = nnod + 1
        g1sin(n-1) = g1sng1
        f1sin(n-1) = f1snf1
        g1snm4 = g1snm3
        f1snm4 = f1snm3
        g1snm3 = g1snm2
        f1snm3 = f1snm2
        g1snm2 = g1snm1
        f1snm2 = f1snm1
        g1snm1 = g1sn
        f1snm1 = f1sn
        g1sn = g1sng1
        f1sn = f1snf1

!        dg1new = g1in(n-1) - g1in(n)
!        dg1old = g1in(n) - g1in(n+1)
!        fr = dg1new/dg1old

!        if ((fr<0).and.(nfr==0)) then
!           nfr=nfr+1
!           ncsave = n
!        endif

!        if ((alpha==1).and.(g1in(n) /= 0)) then
!           if ((g1in(n+1)/g1in(n))<=0) nnod = nnod + 1
!        endif

        g2in(n-1) = g2c
        f2in(n-1) = f2c
        g2sin(n-1) = g2sng2
        f2sin(n-1) = f2snf2
        g2snm4 = g2snm3
        f2snm4 = f2snm3
        g2snm3 = g2snm2
        f2snm3 = f2snm2
        g2snm2 = g2snm1
        f2snm2 = f2snm1
        g2snm1 = g2sn
        f2snm1 = f2sn
        g2sn = g2sng2
        f2sn = f2snf2
!        n = n-1
      enddo

!      if (nfr==0) ncsave=1450
!      nc = n + 1

       nc = 1360 !1030 !ncsave
!      nrk = 6 ; nrmx = nc
       nrk = 6 ; nrmx = nr

      else
!________________alpha=2____________________________________-

        do n = nr-5, 2, -1
          J = a*(rofi(n) + b); r = rofi(n)
!C          J = J*q*q
!C          r = J/a - b

          g1p = g1in(n+5) - (3*dx/10.0_8)*(11*g1sn - 14*g1snm1 + 26*g1snm2 - 14*g1snm3 + 11*g1snm4)
          f1p = f1in(n+5) - (3*dx/10.0_8)*(11*f1sn - 14*f1snm1 + 26*f1snm2 - 14*f1snm3 + 11*f1snm4)
          g2p = g2in(n+5) - (3*dx/10.0_8)*(11*g2sn - 14*g2snm1 + 26*g2snm2 - 14*g2snm3 + 11*g2snm4)
          f2p = f2in(n+5) - (3*dx/10.0_8)*(11*f2sn - 14*f2snm1 + 26*f2snm2 - 14*f2snm3 + 11*f2snm4)

          eps1 = dabs(g1p-g1c)
          eps2 = dabs(f1p-f1c)
          eps3 = dabs(g2p-g2c)
          eps4 = dabs(f2p-f2c)

          nit = 0
          Delt = 1d-12
          do while (eps1>dabs(g1c)*Delt .or. eps2>dabs(f1c)*Delt &
               .or. eps3>dabs(g2c)*Delt .or. eps4>dabs(f2c)*Delt)

             g1sng1 = J*(c*nug1(n-1)*f1p - xk/r*g1p)
             f1snf1 = J*((nuf1(n-1)*g1p - sqru*bm(n-1)*g2p)/c + xk/r*f1p)
             g2sng2 = J*(c*nug2(n-1)*f2p - xk1/r*g2p)
             f2snf2 = J*((nuf2(n-1)*g2p - sqru*bm(n-1)*g1p)/c + xk1/r*f2p)

             g1c = g1in(n+3) - (2*dx/45.0_8)*(7*g1sng1 + 32*g1sn + 12*g1snm1 + 32*g1snm2 + 7*g1snm3)
             f1c = f1in(n+3) - (2*dx/45.0_8)*(7*f1snf1 + 32*f1sn + 12*f1snm1 + 32*f1snm2 + 7*f1snm3)
             g2c = g2in(n+3) - (2*dx/45.0_8)*(7*g2sng2 + 32*g2sn + 12*g2snm1 + 32*g2snm2 + 7*g2snm3)
             f2c = f2in(n+3) - (2*dx/45.0_8)*(7*f2snf2 + 32*f2sn + 12*f2snm1 + 32*f2snm2 + 7*f2snm3)

             eps1 = dabs(g1c-g1p)
             eps2 = dabs(f1c-f1p)
             eps3 = dabs(g2c-g2p)
             eps4 = dabs(f2c-f2p)

! ...      Check for convergence
!          if (nit == 100) then
!            call info5(1,0,0,'rdq2 (warning): not converged for '//
!     .        'alpha=%i  ie=%i n=%i  kappa=%d  mu=%d',alpha,ie,n,xk,mu)
!           nprint = nprint + 1
!            exit
!          endif

            g1p = g1c
            f1p = f1c
            g2p = g2c
            f2p = f2c

!            nit = nit + 1
        enddo

!        if ((g1inmax <= dabs(g1c)).and.(n>1300)) g1inmax = dabs(g1c)
        g1in(n-1) = g1c
        f1in(n-1) = f1c

        g1sin(n-1) = g1sng1
        f1sin(n-1) = f1snf1
        g1snm4 = g1snm3
        f1snm4 = f1snm3
        g1snm3 = g1snm2
        f1snm3 = f1snm2
        g1snm2 = g1snm1
        f1snm2 = f1snm1
        g1snm1 = g1sn
        f1snm1 = f1sn
        g1sn = g1sng1
        f1sn = f1snf1

!        dg1new = g1in(n-1) - g1in(n)
!        dg1old = g1in(n) - g1in(n+1)
!        fr = dg1new/dg1old

        g2in(n-1) = g2c
        f2in(n-1) = f2c
        nodsig = g2in(n-1)*g2in(n)
        if ((n >= nmrseq).and.(nodsig <= 0)) nnod = nnod + 1
        g2sin(n-1) = g2sng2
        f2sin(n-1) = f2snf2
        g2snm4 = g2snm3
        f2snm4 = f2snm3
        g2snm3 = g2snm2
        f2snm3 = f2snm2
        g2snm2 = g2snm1
        f2snm2 = f2snm1
        g2snm1 = g2sn
        f2snm1 = f2sn
        g2sn = g2sng2
        f2sn = f2snf2
!        n = n-1
      enddo


!_______________________________End Integration inwards

!_______________________________Start integration outwards
!      nc = n + 1
!      print*,'nc = ',nc


!       nrk = 6 ; nrmx = nc



      endif
!________________end alpha=2____________________________________

!___________________Remnant of the initial rdeq_________________
!C ... Save and normalize the Psi functions
!      do jalp=1,2
!      do ie=1,5
!        psiin(:,1,jalp,ie) = g1in(:)
!        psiin(:,2,jalp,ie) = f1in(:)
!        psiin(:,3,jalp,ie) = g2in(:)
!        psiin(:,4,jalp,ie) = f2in(:)
!      enddo
!      enddo
!      do ie=1,5
!        psiinn(:,1,alpha,ie) = g1in(:)
!        psiinn(:,2,alpha,ie) = f1in(:)
!        psiinn(:,3,alpha,ie) = g2in(:)
!        psiin(:,4,alpha,ie) = f2in(:)
!      enddo



!___________________End of Remnant of the initial rdeq__________

!___________________For normalization___________________________

      psiinn(:,1)=g1in(:)
      psiinn(:,2)=f1in(:)
      psiinn(:,3)=g2in(:)
      psiinn(:,4)=f2in(:)

!__________________End _For normalization________________________

!__________________Normalization_________________________________

!!!      ncut=1400*nr/1501
!!!      ncut = INT(ncut)

       deltacut = 100*nr/1501
       deltacut = INT(deltacut)
       ncut1 = nmrseq - deltacut
       ncut2 = nmrseq + deltacut

!      call productn5(psioutn(:,:),psioutn(:,:),a,b,rofi,2,nmrseq,normn)
!      do j1=1,4
!         psioutn(:,j1)=psioutn(:,j1)/dsqrt(normn)
!      enddo

!!      do j1=1,4
!!        call productn11(psioutn(:,j1),psioutn(:,j1),a,b,rofi,nr,ncut1,ncut2,normn)
!!        if (normn /= 0) psioutn(:,j1)=psioutn(:,j1)/dsqrt(normn)
!!      enddo

      normsave = dsqrt(normn)

!!!      ncut = 1250*nr/1501
!!!      ncut = INT(ncut)

      frac = 1
      if (psiinn(nmrseq,1) /= 0) then
         if ((psioutn(nmrseq,1)/psiinn(nmrseq,1)) <= 0) frac=-1
      endif

!      call productn5(psiinn(:,:),psiinn(:,:),a,b,rofi,nmrseq,nr,normni)
!      do j1=1,4
!         psiinn(:,j1)=psiinn(:,j1)*frac/dabs(dsqrt(normni))
!      enddo

!      do j1=1,4
!        frac = 1
!        if (psiinn(nmrseq,j1) /= 0) then
!           if ((psioutn(nmrseq,j1)/psiinn(nmrseq,j1)) <= 0) frac=-1
!        endif
!        psiinn(:,j1)=psiinn(:,j1)*frac
!      enddo

      g1out(:)=psioutn(:,1)
      f1out(:)=psioutn(:,2)
      g2out(:)=psioutn(:,3)
      f2out(:)=psioutn(:,4)
      g1in(:)=psiinn(:,1)
      f1in(:)=psiinn(:,2)
      g2in(:)=psiinn(:,3)
      f2in(:)=psiinn(:,4)

!      call productn5(psiinn(:,:),psiinn(:,:),a,b,rofi,nc,nr,norm1)
!      call productn5(psioutn(:,:),psioutn(:,:),a,b,rofi,2,nc,norm2)
!      normn = norm1+norm2
!      do j1=1,4
!         psiinn(:,j1)=psiinn(:,j1)/dsqrt(normn)
!         psioutn(:,j1)=psioutn(:,j1)/dsqrt(normn)
!      enddo



!__________________End Normalization_______________________________


!__________________Scaling by g1(nc)________________________________
!      print*,'g1inmax = ',g1inmax
!      print*,'g1outmax = ',g1outmax

!      do j2=1,nr
!         do j1=1,4
!           if (g1inmax /= 0) psiinn(j2,j1)=psiinn(j2,j1)/g1inmax
!           if (g1outmax /= 0) psioutn(j2,j1)=psioutn(j2,j1)/g1outmax
!         enddo
!      enddo

!      g1out(:)=psioutn(:,1)
!      f1out(:)=psioutn(:,2)
!      g2out(:)=psioutn(:,3)
!      f2out(:)=psioutn(:,4)
!      g1in(:)=psiinn(:,1)
!      f1in(:)=psiinn(:,2)
!      g2in(:)=psiinn(:,3)
!      f2in(:)=psiinn(:,4)

!__________________End of Scaling by g1(nc)____________________________


!      if (alpha == 1) call printdata2('out_OutInt', 'rofi g1 12 12 2013', nr-1, rofi(2:), g1out(2:))

!      do ie=1,5
!        psiout(:,1,alpha,ie) = g1out(:)
!        psiout(:,2,alpha,ie) = f1out(:)
!        psiout(:,3,alpha,ie) = g2out(:)
!        psiout(:,4,alpha,ie) = f2out(:)
!        psiin(:,1,alpha,ie) = g1in(:)
!        psiin(:,2,alpha,ie) = f1in(:)
!        psiin(:,3,alpha,ie) = g2in(:)
!        psiin(:,4,alpha,ie) = f2in(:)
!      enddo

      end subroutine integrationC

!___________F(P)=0______________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      subroutine newtC(Pv,nF,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq,check,fvector,maxi)

      use newtvC_common
      implicit none
!C____________________________________________________________________________
!C Given an initial guess for Pv(E,A,Aout,Ain) (see Ebert,J.Phys.:Cond.M.1(1989)
!C  9111-9116) finds the root of F = 0 (Ebert) by a globally convergent
!C Newton's method. Taken from 'Numerical Recipes in fortran 77, p 379.
!C----------------------------------------------------------------------------

!      integer nF,MAXITS
      integer nF,MAXITS
      logical check, maxi
      double precision Pv(nF),TOLF,TOLMIN,TOLX,STPMX

      parameter (MAXITS=300,TOLF=1.e-6_8,TOLMIN=1.e-6_8, TOLX=1.e-7_8, STPMX=100._8)
!       common /newtv/ fvecC(NP),
!      . goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP),nn
!       save /newtv/
!CU    USES fdjacC,fminC,lnsrchC,lubksb,ludcmp
      integer i,its,j,indx(NP),nmrseq
      double precision d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),xold(NP), fvecold(NP) !,fminC
      integer nr,nsp,imu,l,lmx,nm
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision  grseq(nr,2)
      double precision mu,rm,sqru,u,up,upp,xk,xk1
      double precision  bm(nr),vm(nr)

      double precision deltaf(NP), fvector(NP), deltaE, Eold
      procedure(real(8)) :: fminC
!       EXTERNAL fminC

!       if (l == 2) print*,'newt - 1'
      maxi=.false.
      nn=nF
      f=fminC(Pv,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx,imu,scalede,grseq)
      test=0.0_8
      do i=1,nF
        if(abs(fvecC(i)) > test)test=abs(fvecC(i))
      enddo
      if(test < .01_8*TOLF)then
        check=.false.
        return
      endif
      sum=0.0_8
      do i=1,nF
        sum=sum+Pv(i)**2
      enddo
      stpmax=STPMX*max(sqrt(sum),float(nF))

      deltaE = 1.0_8
      deltaf(1) = 1000.0_8; deltaf(2) = 10000.0_8
      deltaf(3) = 1000.0_8; deltaf(4) = 20000.0_8

      do its=1,MAXITS
        call fdjacC(nF,Pv,fjac, xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,grseq)
!        print*,'Jacobian = ',fjac
        do i=1,nF
          sum=0.0_8
          do j=1,nF
            sum=sum+fjac(j,i)*fvecC(j)
          enddo
          g(i)=sum
        enddo
        do i=1,nF
          xold(i)=Pv(i)
          fvecold(i)=fvecC(i)
        enddo
        Eold = Pv(1)
        fold=f
        do i=1,nF
          p(i)=-fvecC(i)
        enddo

        call ludcmp(fjac,nF,NP,indx,d)
        call lubksb(fjac,nF,NP,indx,p)
        call lnsrchC(nF,xold,fold,g,p,Pv,f,stpmax,check,fminC,  &
            xk, xk1, mu, up, u,                               &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,    &
            lmx, imu,scalede,grseq)


        fvector(:) = fvecC(:)
!      if ((l==1).and.(imu==2)) then
!      print*,'? '
!      print*,'Final values'
!      print*, 'newt - end','Energy=',Pv(1),' A = ',Pv(2), ' Aout = ', Pv(3),' Ain = ', Pv(4)
!      print*,'ln 2273 newt F = ',fvecC
!      endif

        deltaE = Pv(1)- Eold
        do i=1,nF
          deltaf(i) = fvecC(i)-fvecold(i)
        enddo

        test=0.0_8
        do i=1,nF
          if(abs(fvecC(i)) > test)test=abs(fvecC(i))
        enddo
        if(test < TOLF)then
          write(1000,*)fvecC(:),', test = ',test,'TOL=',TOLF
          check=.false.
          return
        endif
        if(check)then
          test=0.0_8
          den=max(f,.5_8*nF)
          do i=1,nF
            temp=abs(g(i))*max(abs(Pv(i)),1.)/den
            if(temp > test)test=temp
          enddo
          if(test < TOLMIN)then
            check=.true.
          else
            check=.false.
            write(1000,*)fvecC(:),', test = ',test,'TOLM=',TOLMIN
          endif
          return
        endif
        test=0.0_8
        do i=1,nF
          temp=(abs(Pv(i)-xold(i)))/max(abs(Pv(i)),1.)
          if(temp > test)test=temp
        enddo
        if(test < TOLX)return
      enddo


      maxi =.true.
!      call pause('MAXITS exceeded in newtC')
      print*,'MAXITS exceeded in newtC'
      end subroutine newtC

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------
!      call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin0)
      subroutine lnsrchC(n,xold,fold,g,p,Pv,f,stpmax, &
        check,func, xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,ar,br,l, lmx, imu,scalede,&
        grseq)
        implicit none
!     n-dimensions
      integer n,nmrseq
      logical check
      double precision f,fold,stpmax,g(n),p(n),Pv(n),xold(n),func, ALF,TOLX
      parameter (ALF=1.e-4_8,TOLX=1.e-7_8)
      external func
!CU    USES func
      integer i
      double precision a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,   &
        slope,sum,temp,test,                                      &
        tmplam,                                                   &
        goutnm(n), foutnm(n), ginnm(n), finnm(n),                 &
        goutnms(n), foutnms(n), ginnms(n), finnms(n)
      double precision xk, xk1, mu, up, u, upp, sqru
      double precision vm(nr), bm(nr)
      integer nr,nsp,imu,l,lmx
      double precision v(nr,nsp),rofi(nr),z,ar,br,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision grseq(nr,2)

      check=.false.
      sum=0.0_8
      do i=1,n
        sum=sum+p(i)*p(i)
      enddo
      sum=sqrt(sum)
      if(sum > stpmax)then
        do i=1,n
          p(i)=p(i)*stpmax/sum
        enddo
      endif
      slope=0.0_8
      do i=1,n
        slope=slope+g(i)*p(i)
      enddo
      if(slope >= 0.0_8) call pause('roundoff problem in lnsrchC')
      test=0.0_8
      do i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.0_8)
        if(temp > test)test=temp
      enddo
      alamin=TOLX/test
      alam=1.0_8
1     continue

        do i=1,n
          Pv(i)=xold(i)+alam*p(i)
        enddo

        f=func(Pv,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,ar,br,l,lmx, imu,scalede,&
               grseq)

        if(alam < alamin)then
          do i=1,n
            Pv(i)=xold(i)
          enddo
          check=.true.

          return
        else if(f <= fold+ALF*alam*slope)then
          return
        else
          if(alam == 1.0_8)then
            tmplam=-slope/(2*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a == 0.0_8)then
              tmplam=-slope/(2.0_8*b)
            else
              disc=b*b-3*a*slope
              if(disc < 0.0_8)then
                tmplam=.5_8*alam
              else if(b <= 0.0_8)then
                tmplam=(-b+sqrt(disc))/(3.0_8*a)
              else
                tmplam=-slope/(b+sqrt(disc))
              endif
            endif
            if(tmplam > .5_8*alam)tmplam=.5_8*alam
          endif
        endif
        alam2=alam
        f2=f
        alam=max(tmplam,.1_8*alam)

      goto 1

      end subroutine lnsrchC

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      function fminC(Pv,xk, xk1, mu, up, u,upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,&
                grseq)

      use newtvC_common
      implicit none
!       integer n
      double precision fminC,Pv(*)
!       ,fvecC,
!      . goutnm, foutnm, ginnm, finnm
!       common /newtv/fvecC(NP),
!      . goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP),n
!       save /newtv/
      integer i,nmrseq
      double precision sum
      integer nr,nsp,imu,l,lmx,nm
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision grseq(nr,2)
      double precision mu,rm,sqru,u,up,upp,xk,xk1
      double precision  bm(nr),vm(nr)

!       if (l == 2) print*,'fminC - 1'

      call funcvC(nn,Pv,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l, lmx, imu,scalede,grseq,fvecC)!, goutnm, foutnm, ginnm, finnm)

!       if (l == 2) print*,'fminC - 2'

       sum=0.0_8
       do i=1,nn
         sum = sum + fvecC(i)**2
       enddo
       fminC=0.5_8*sum

       return
      end function fminC

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

!      subroutine ludcmp(a,n,np,indx,d)
!      implicit none
!      integer n,np,indx(n),NMAX
!      double precision d,a(np,np),TINY
!      parameter (NMAX=500,TINY=1.0e-20_8)
!      integer i,imax,j,k
!      double precision aamax,dum,sum,vv(NMAX)

!!      print*,'ludcmp - n,np', n,'; ',np
!      d=1.0_8
!      do i=1,n
!        aamax=0.0_8
!        do j=1,n
!          if (abs(a(i,j)) > aamax) aamax=abs(a(i,j))
!        enddo
!        if (aamax == 0.0_8) call pause('singular matrix in ludcmp')
!        vv(i)=1./aamax
!      enddo
!      do  j=1,n
!        do i=1,j-1
!          sum=a(i,j)
!          do k=1,i-1
!            sum=sum-a(i,k)*a(k,j)
!          enddo
!          a(i,j)=sum
!        enddo
!        aamax=0.0_8
!        do  i=j,n
!          sum=a(i,j)
!          do  k=1,j-1
!            sum=sum-a(i,k)*a(k,j)
!          enddo
!          a(i,j)=sum
!          dum=vv(i)*abs(sum)
!          if (dum >= aamax) then
!            imax=i
!            aamax=dum
!          endif
!        enddo
!        if (j /= imax)then
!          do  k=1,n
!            dum=a(imax,k)
!            a(imax,k)=a(j,k)
!            a(j,k)=dum
!          enddo
!          d=-d
!          vv(imax)=vv(j)
!        endif
!        indx(j)=imax
!        if(a(j,j) == 0.0_8)a(j,j)=TINY
!        if(j /= n)then
!          dum=1.0_8/a(j,j)
!          do  i=j+1,n
!            a(i,j)=a(i,j)*dum
!          enddo
!        endif
!      enddo
!      return
!      END subroutine ludcmp

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

!      subroutine lubksb(a,n,np,indx,b)
!      implicit none
!      integer n,np,indx(n)
!      double precision a(np,np),b(n)
!      integer i,ii,j,ll
!      double precision sum

!!      print*,'lubksb - n,np', n,'; ',np

!      ii=0
!      do i=1,n
!        ll=indx(i)
!        sum=b(ll)
!        b(ll)=b(i)
!        if (ii /= 0)then
!          do j=ii,i-1
!            sum=sum-a(i,j)*b(j)
!          enddo
!        else if (sum /= 0.0_8) then
!          ii=i
!        endif
!        b(i)=sum
!      enddo
!      do i=n,1,-1
!        sum=b(i)
!        do j=i+1,n
!          sum=sum-a(i,j)*b(j)
!        enddo
!        b(i)=sum/a(i,i)
!      enddo
!      return
!      END subroutine lubksb


!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      subroutine funcvC(n,Pv,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq,fvecC)!, goutnm, foutnm, ginnm, finnm)
       implicit none
!       integer n
       double precision Pv(1:n),fvecC(1:n)
       integer nr,nsp,imu,l,lmx,nm, nc,nmrseq
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6)
       double precision  grseq(nr,2)

       integer n, nnod
       double precision mu,rm,sqru,u,up,upp,xk,xk1,En,normsave1,normsave2
!    .  c,csq               !Speed of light 274.072d0


       double precision  g1out1(nr),g2out1(nr),f1out1(nr),f2out1(nr),    &
         g1sout1(nr),g2sout1(nr),f1sout1(nr),f2sout1(nr),g1in1(nr),      &
         g2in1(nr),f1in1(nr),f2in1(nr),g1sin1(nr),g2sin1(nr),f1sin1(nr), &
         f2sin1(nr),bm(nr),vm(nr),                                       &
         g1out2(nr),g2out2(nr),f1out2(nr),f2out2(nr),                    &
         g1sout2(nr),g2sout2(nr),f1sout2(nr),f2sout2(nr),g1in2(nr),      &
         g2in2(nr),f1in2(nr),f2in2(nr),g1sin2(nr),g2sin2(nr),f1sin2(nr), &
         f2sin2(nr),                                                     &
         goutnm(4), foutnm(4), ginnm(4), finnm(4),                       &
         psiout1(nr,4,2,5),psiin1(nr,4,2,5),                             &
         psiout2(nr,4,2,5),psiin2(nr,4,2,5)

        call integrationC(1, Pv(1), xk, xk1, mu, up, u,           &
            upp, sqru, vm,                                        &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede, grseq,&
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,  &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,&
            f2sin1,nmrseq, nnod,normsave1)

        call integrationC(2, Pv(1), xk, xk1, mu, up, u,            &
            upp, sqru, vm,                                         &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede, grseq,&
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,   &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2, &
            f2sin2,nmrseq, nnod,normsave2)

        nm=nmrseq

        goutnm(1)=g1out1(nm); goutnm(2)=g2out1(nm)
        goutnm(3)=g1out2(nm); goutnm(4)=g2out2(nm)

        foutnm(1)=f1out1(nm); foutnm(2)=f2out1(nm)
        foutnm(3)=f1out2(nm); foutnm(4)=f2out2(nm)

        ginnm(1)=g1in1(nm); ginnm(2)=g2in1(nm)
        ginnm(3)=g1in2(nm); ginnm(4)=g2in2(nm)

        finnm(1)=f1in1(nm); finnm(2)=f2in1(nm)
        finnm(3)=f1in2(nm); finnm(4)=f2in2(nm)

        fvecC(1) =  (g1out1(nm) + Pv(3)*g1out2(nm)) - Pv(2)*(g1in1(nm) + Pv(4)*g1in2(nm))
        fvecC(2) =  (f1out1(nm) + Pv(3)*f1out2(nm)) - Pv(2)*(f1in1(nm) + Pv(4)*f1in2(nm))
        fvecC(3) =  (g2out1(nm) + Pv(3)*g2out2(nm)) - Pv(2)*(g2in1(nm) + Pv(4)*g2in2(nm))
        fvecC(4) =  (f2out1(nm) + Pv(3)*f2out2(nm)) - Pv(2)*(f2in1(nm) + Pv(4)*f2in2(nm))

!       if (l == 2) then
!          print*,'funcv - 4'
!          print*,'for debug'
!       endif


!        print*,' '
!        print*,'line 2235 f = ', fvecC

      end subroutine funcvC

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------


      subroutine fdjacC(n,Pv,df, xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq)
       use newtvC_common
     implicit none
       integer, intent(in) :: n
       integer nc, nm, nnod,nmrseq
       integer nr, nsp, imu, l, lmx
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6)
       double precision  grseq(nr,2)
       double precision mu, u, up, upp, xk, xk1, sqru
!        double precision fvecC,goutnm, foutnm, ginnm, finnm
!        parameter (NP=4)
!        common /newtv/ fvecC(NP),
!      . goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP)!,n
!       save /newtv/
       double precision Pv(1:n), df(n,n)
       double precision fvecp(1:n), fvecm(1:n)
       double precision deltaE,Eplus, Eminus
       double precision  g1out1m(nr),g2out1m(nr),f1out1m(nr),         &
         f2out1m(nr),g1out1p(nr),g2out1p(nr),f1out1p(nr),             &
         f2out1p(nr),                                                 &
         g1sout1(nr),g2sout1(nr),f1sout1(nr),f2sout1(nr),g1in1m(nr),  &
         g2in1m(nr),f1in1m(nr),f2in1m(nr),g1sin1(nr),g2sin1(nr),      &
         f1sin1(nr),g1in1p(nr),                                       &
         g2in1p(nr),f1in1p(nr),f2in1p(nr),                            &
         f2sin1(nr),bm(nr),vm(nr),                                    &
         g1out2m(nr),g2out2m(nr),f1out2m(nr),f2out2m(nr),             &
         g1out2p(nr),g2out2p(nr),f1out2p(nr),f2out2p(nr),             &
         g1sout2(nr),g2sout2(nr),f1sout2(nr),f2sout2(nr),g1in2m(nr),  &
         g2in2m(nr),f1in2m(nr),f2in2m(nr),g1sin2(nr),g2sin2(nr),      &
         f1sin2(nr),g1in2p(nr),                                       &
         g2in2p(nr),f1in2p(nr),f2in2p(nr),                            &
         f2sin2(nr),                                                  &
         psiout1(nr,4,2,5),psiin1(nr,4,2,5),                          &
         psiout2(nr,4,2,5),psiin2(nr,4,2,5)
!     .  goutnm(4), foutnm(4), ginnm(4), finnm(4),

        double precision  g1out1(nr),g2out1(nr),f1out1(nr), f2out1(nr),  &
       g1in1(nr), g2in1(nr),f1in1(nr),f2in1(nr),                         &
       g1out2(nr),g2out2(nr),f1out2(nr),f2out2(nr),                      &
       g1in2(nr), g2in2(nr),f1in2(nr),f2in2(nr),normsave1,normsave2


         deltaE = Pv(1)*1e-3
         Eplus  = Pv(1)+deltaE
         Eminus = Pv(1)-deltaE
         call integrationC(1, Eplus, xk, xk1, mu, up, u,                     &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,         &
            g1out1p,f1out1p,g2out1p,f2out1p,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1p,f1in1p,g2in1p,f2in1p,g1sin1,f1sin1,g2sin1,       &
            f2sin1,nmrseq, nnod,normsave1)
         call integrationC(2, Eplus, xk, xk1, mu, up, u,                     &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,         &
            g1out2p,f1out2p,g2out2p,f2out2p,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2p,f1in2p,g2in2p,f2in2p,g1sin2,f1sin2,g2sin2,       &
            f2sin2,nmrseq, nnod,normsave2)
         nm=nmrseq
         fvecp(1) =  (g1out1p(nm) + Pv(3)*g1out2p(nm)) - Pv(2)*(g1in1p(nm)   &
              + Pv(4)*g1in2p(nm))
         fvecp(2) =  (f1out1p(nm) + Pv(3)*f1out2p(nm)) - Pv(2)*(f1in1p(nm)   &
              + Pv(4)*f1in2p(nm))
         fvecp(3) =  (g2out1p(nm) + Pv(3)*g2out2p(nm)) - Pv(2)*(g2in1p(nm)   &
              + Pv(4)*g2in2p(nm))
         fvecp(4) =  (f2out1p(nm) + Pv(3)*f2out2p(nm)) - Pv(2)*(f2in1p(nm)   &
              + Pv(4)*f2in2p(nm))

         call integrationC(1, Eminus, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,         &
            g1out1m,f1out1m,g2out1m,f2out1m,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1m,f1in1m,g2in1m,f2in1m,g1sin1,f1sin1,g2sin1,       &
            f2sin1,nmrseq, nnod,normsave1)
         call integrationC(2, Eminus, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,         &
            g1out2m,f1out2m,g2out2m,f2out2m,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2m,f1in2m,g2in2m,f2in2m,g1sin2,f1sin2,g2sin2,       &
            f2sin2,nmrseq, nnod,normsave2)
         nm=nmrseq
         fvecm(1) =  (g1out1m(nm) + Pv(3)*g1out2m(nm)) - Pv(2)*(g1in1m(nm)   &
              + Pv(4)*g1in2m(nm))
         fvecm(2) =  (f1out1m(nm) + Pv(3)*f1out2m(nm)) - Pv(2)*(f1in1m(nm)   &
              + Pv(4)*f1in2m(nm))
         fvecm(3) =  (g2out1m(nm) + Pv(3)*g2out2m(nm)) - Pv(2)*(g2in1m(nm)   &
              + Pv(4)*g2in2m(nm))
         fvecm(4) =  (f2out1m(nm) + Pv(3)*f2out2m(nm)) - Pv(2)*(f2in1m(nm)   &
              + Pv(4)*f2in2m(nm))

         df(1,1) = (fvecp(1)-fvecm(1))/2/deltaE
         df(2,1) = (fvecp(2)-fvecm(2))/2/deltaE
         df(3,1) = (fvecp(3)-fvecm(3))/2/deltaE
         df(4,1) = (fvecp(4)-fvecm(4))/2/deltaE

         call integrationC(1, Pv(1), xk, xk1, mu, up, u,                     &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,         &
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,             &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,           &
            f2sin1,nmrseq, nnod,normsave1)
         call integrationC(2, Pv(1), xk, xk1, mu, up, u,                     &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,  grseq,         &
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,             &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,           &
            f2sin2,nmrseq, nnod,normsave2)

           nc=nmrseq
          df(1,2) = -g1in1(nc)-Pv(4)*g1in2(nc)
          df(1,3) = g1out2(nc)
          df(1,4) = -Pv(2)*g1in2(nc)

          df(2,2) = -f1in1(nc)-Pv(4)*f1in2(nc)
          df(2,3) = f1out2(nc)
          df(2,4) = -Pv(2)*f1in2(nc)

          df(3,2) = -g2in1(nc)-Pv(4)*g2in2(nc)
          df(3,3) = g2out2(nc)
          df(3,4) = -Pv(2)*g2in2(nc)

          df(4,2) = -f2in1(nc)-Pv(4)*f2in2(nc)
          df(4,3) = f2out2(nc)
          df(4,4) = -Pv(2)*f2in2(nc)

      end subroutine fdjacC

!_______end of F(Pv)=0________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------


!______Another equations system solver____________________________


      subroutine usrfun1C(Pv,n,NP,fvecC,fjac, xk, xk1, mu, up, u,  &
            upp, sqru, vm,                                       &
            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq)
        implicit none
       integer n, NP,nmrseq
       integer nr, nsp, imu, l, lmx
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6),Pv(1:n)
       double precision grseq(nr,2)
       double precision mu, u, up, upp, xk, xk1, sqru
!      parameter (NP=4)
       double precision fjac(n,n),fvecC(n),bm(nr),vm(nr)
       double precision goutnm(4), foutnm(4), ginnm(4), finnm(4)


!       print*,'In usrfun'

       call funcvC(n,Pv,xk, xk1, mu, up, u,                           &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,          &
            lmx, imu,scalede, grseq,fvecC)!, goutnm, foutnm, ginnm, finnm)

!       print*,'After funcvC'
!       print*,'Pv = ',Pv

       call fdjacC(n,Pv,fjac, xk, xk1, mu, up, u, &
            upp, sqru, vm,                       &
            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,grseq)

!       print*,'After fdjac'

      end subroutine usrfun1C

      subroutine mnewtC(ntrial,x,n,tolx,tolf, xk, xk1, mu, up, u, &
             upp, sqru, vm,                                      &
             bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq,fvecC)
      implicit none
       integer n, NP, ntrial, nmrseq
       double precision tolf, tolx, x(n)
       parameter (NP=4)
       double precision d,errf, errx,fjac(NP,NP),fvecC(NP),p(NP)
       integer i, k, indx(NP)
       integer nr, nsp, imu, l, lmx
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6)
       double precision grseq(nr,2)
       double precision mu, u, up, upp, xk, xk1, sqru
!      parameter (NP=4)
       double precision Pv(n),bm(nr),vm(nr)

!       print*,'In mnewtC'

       do k=1,ntrial
          call usrfun1C(x,n,NP,fvecC,fjac, xk, xk1, mu, up, u, &
             upp, sqru, vm,                                  &
             bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq)
          errf = 0
          do i=1,n
             errf = errf+abs(fvecC(i))
          enddo

!          write(3,*)' k = ',k
!          write(3,*)' E = ',x(1)
!          write(3,*)' F = ',fvecC(:)

          if (errf <= tolf) return
          do i=1,n
             p(i)=-fvecC(i)
          enddo
          call ludcmp(fjac,n,NP,indx,d)
          call lubksb(fjac,n,NP,indx,p)


          errx = 0
          do i=1,n
             errx=errx+abs(p(i))
             x(i)=x(i)+p(i)
          enddo


          if (errx <= tolx) return
       enddo
       return

      end subroutine mnewtC




!______Another equations system solver done________________________



!______Third monomization procedure________________________________

      subroutine broydnC(x,n,check, xk, xk1, mu, up, u, upp, sqru, vm, &
             bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq,fsave)
        use newtvC_common
      implicit none
!     x = Pv

      integer n,nmrseq
      integer, parameter :: MAXITS=200
      double precision x(n)
      real(8), parameter :: EPS=1.e-7_8,TOLF=1.e-9_8,TOLMIN=1.e-8_8, STPMX=100.0_8
      real(8), parameter :: TOLX=EPS
      logical check
!       common /newtv/ fvecC(NP),nn
!       save /newtv/
!CU    uses fdjacC,fminC,lnsrchC,qrdcmp,qrupdt,rsolv
      integer i,its,j,k
      double precision den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP), &
             g(NP),p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP) !,fminC
      logical restrt,sing,skip
!       external fminC

      integer nr, nsp, imu, l, lmx
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision  grseq(nr,2)
      double precision mu, u, up, upp, xk, xk1, sqru
      double precision Pv(n),bm(nr),vm(nr),fsave(n)

      procedure(real(8)) :: fminC

      nn=n
      f=fminC(x,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l, lmx, imu,scalede,grseq)
      test=0.0_8
      do i=1,n
        if(abs(fvecC(i)) > test)test=abs(fvecC(i))
      enddo
      if(test < .01_8*TOLF)then
        check=.false.
        return
      endif
      sum=0.0_8
      do i=1,n
        sum=sum+x(i)**2
      enddo
      stpmax=STPMX*max(sqrt(sum),float(n))
      restrt=.true.
      do its=1,MAXITS
        if(restrt)then
          call fdjacC(n,x,r,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede, grseq)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) call pause('singular Jacobian in broydnC')
          do i=1,n
            do j=1,n
              qt(i,j)=0.0_8
            enddo
            qt(i,i)=1.
          enddo
          do k=1,n-1
            if(c(k) /= 0.)then
              do j=1,n
                sum=0.0_8
                do i=k,n
                  sum=sum+r(i,k)*qt(i,j)
                enddo
                sum=sum/c(k)
                do i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
                enddo
              enddo
            endif
          enddo
          do i=1,n
            r(i,i)=d(i)
            do j=1,i-1
              r(i,j)=0.0_8
            enddo
          enddo
        else
          do i=1,n
            s(i)=x(i)-xold(i)
          enddo
          do i=1,n
            sum=0.0_8
            do  j=i,n
              sum=sum+r(i,j)*s(j)
            enddo
            t(i)=sum
          enddo
          skip=.true.
          do i=1,n
            sum=0.0_8
            do  j=1,n
              sum=sum+qt(j,i)*t(j)
            enddo
            w(i)=fvecC(i)-fvcold(i)-sum
            if(abs(w(i)) >= EPS*(abs(fvecC(i))+abs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.0_8
            endif
          enddo
          if(.not.skip)then
            do i=1,n
              sum=0.0_8
              do j=1,n
                sum=sum+qt(i,j)*w(j)
              enddo
              t(i)=sum
            enddo
            den=0.0_8
            do i=1,n
              den=den+s(i)**2
            enddo
            do i=1,n
              s(i)=s(i)/den
            enddo
            call qrupdt(r,qt,n,NP,t,s)
            do i=1,n
              if(r(i,i) == 0.0_8) call pause('r singular in broydnC')
              d(i)=r(i,i)
            enddo
          endif
        endif
        do i=1,n
          sum=0.0_8
          do  j=1,n
            sum=sum+qt(i,j)*fvecC(j)
          enddo
          p(i)=-sum
        enddo
        do i=n,1,-1
          sum=0.0_8
          do j=1,i
            sum=sum-r(j,i)*p(j)
          enddo
          g(i)=sum
        enddo
        do i=1,n
          xold(i)=x(i)
          fvcold(i)=fvecC(i)
        enddo
        fold=f
!        do  i=1,n
!          sum=0.
!          do j=1,n
!            sum=sum+qt(i,j)*fvecC(j)
!          enddo
!          p(i)=-sum
!        enddo
        call rsolv(r,n,NP,d,p)
        call lnsrchC(n,xold,fold,g,p,x,f,stpmax,check,fminC, &
            xk, xk1, mu, up, u,                            &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l, &
            lmx, imu,scalede, grseq)

        fsave(:) = fvecC(:)

        test=0.0_8
        do  i=1,n
          if(abs(fvecC(i)) > test)test=abs(fvecC(i))
        enddo
        if(test < TOLF)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.0_8
            den=max(f,.5_8*n)
            do i=1,n
              temp=abs(g(i))*max(abs(x(i)),1.)/den
              if(temp > test)test=temp
            enddo
            if(test < TOLMIN)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.0_8
          do  i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.0_8)
            if(temp > test)test=temp
          enddo
          if(test < TOLX)return
        endif
      enddo
      call pause('MAXITS exceeded in broydnC')
      END subroutine broydnC

!      subroutine qrdcmp(a,n,np,c,d,sing)
!      implicit none
!      integer n,np
!      double precision a(np,np),c(n),d(n)
!      logical sing
!      integer i,j,k
!      double precision scale,sigma,sum,tau

!      sing=.false.
!      do  k=1,n-1
!        scale=0.0_8
!        do  i=k,n
!          scale=max(scale,abs(a(i,k)))
!        enddo
!        if(scale == 0.0_8)then
!          sing=.true.
!          c(k)=0.0_8
!          d(k)=0.0_8
!        else
!          do  i=k,n
!            a(i,k)=a(i,k)/scale
!          enddo
!          sum=0.0_8
!          do  i=k,n
!            sum=sum+a(i,k)**2
!          enddo
!          sigma=sign(sqrt(sum),a(k,k))
!          a(k,k)=a(k,k)+sigma
!          c(k)=sigma*a(k,k)
!          d(k)=-scale*sigma
!          do  j=k+1,n
!           sum=0.0_8
!            do  i=k,n
!              sum=sum+a(i,k)*a(i,j)
!            enddo
!            tau=sum/c(k)
!            do i=k,n
!              a(i,j)=a(i,j)-tau*a(i,k)
!            enddo
!          enddo
!        endif
!      enddo
!      d(n)=a(n,n)
!      if(d(n) == 0.)sing=.true.
!      return
!      END subroutine qrdcmp

!      subroutine qrupdt(r,qt,n,np,u,v)
!      implicit none
!      integer n,np
!      double precision r(np,np),qt(np,np),u(np),v(np)
!!CU    USES rotate
!      integer i,j,k
!      do  k=n,1,-1
!        if(u(k) /= 0.0_8)goto 1
!      enddo
!      k=1
!1     do i=k-1,1,-1
!        call rotate(r,qt,n,np,i,u(i),-u(i+1))
!        if(u(i) == 0.0_8)then
!          u(i)=abs(u(i+1))
!        else if(abs(u(i)) > abs(u(i+1)))then
!          u(i)=abs(u(i))*sqrt(1.0_8+(u(i+1)/u(i))**2)
!        else
!          u(i)=abs(u(i+1))*sqrt(1.0_8+(u(i)/u(i+1))**2)
!        endif
!      enddo
!      do j=1,n
!        r(1,j)=r(1,j)+u(1)*v(j)
!      enddo
!      do  i=1,k-1
!        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
!      enddo
!      return
!      END subroutine qrupdt

!      subroutine rsolv(a,n,np,d,b)
!      implicit none
!      integer n,np
!      double precision a(np,np),b(n),d(n)
!      integer i,j
!      double precision sum
!      b(n)=b(n)/d(n)
!      do  i=n-1,1,-1
!        sum=0.0_8
!        do  j=i+1,n
!          sum=sum+a(i,j)*b(j)
!        enddo
!        b(i)=(b(i)-sum)/d(i)
!      enddo
!      return
!      END subroutine rsolv

!      subroutine rotate(r,qt,n,np,i,a,b)
!      implicit none
!      integer n,np,i
!      double precision a,b,r(np,np),qt(np,np)
!      integer j
!      double precision c,fact,s,w,y

!      if(a == 0.0_8)then
!        c=0.0_8
!        s=sign(1.0_8,b)
!      else if(abs(a) > abs(b))then
!        fact=b/a
!        c=sign(1.0_8/sqrt(1.0_8+fact**2),a)
!        s=fact*c
!      else
!        fact=a/b
!        s=sign(1.0_8/sqrt(1.0_8+fact**2),b)
!        c=fact*s
!      endif
!      do  j=i,n
!        y=r(i,j)
!        w=r(i+1,j)
!        r(i,j)=c*y-s*w
!        r(i+1,j)=s*y+c*w
!      enddo
!      do  j=1,n
!        y=qt(i,j)
!        w=qt(i+1,j)
!        qt(i,j)=c*y-s*w
!        qt(i+1,j)=s*y+c*w
!      enddo
!      return
!      END subroutine rotate


!______End of the third minimization procedure_____________________

      subroutine productn11(w1,w2,a,b,rofi,nr,n0,nf,normn)
!C-----Integrats wavefunction for normalization purposes-------
      integer n0, nf, n, i,nr
      double precision a,b,rofi(nr), normn
      double precision w1(nr), w2(nr)

      if ((n0==0).or.(n0==1)) n0=2
      normn = 0.0_8
      do  n = n0, nf-1
          f1 = a*(rofi(n  ) + b)/2.0_8
          f2 = a*(rofi(n+1) + b)/2.0_8
          normn = normn + f1*w1(n)*w2(n) + f2*w1(n+1)*w2(n+1)
      enddo
      end subroutine productn11




