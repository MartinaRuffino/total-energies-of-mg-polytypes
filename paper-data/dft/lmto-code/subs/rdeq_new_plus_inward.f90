      module newtv_common
      implicit none

      integer, parameter :: NP=4
      real(8) :: fvec(NP), goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP)
      integer :: nn

      end module newtv_common


!C      subroutine fdpp(enu1,enu2,ksop,shft,gmt,fmt,gmtde,fmtde,
!C     .  z,rmax,avw,l,lmx,imu,
!C     .  sr1,sr2,srav1,srav2,srdot1,srdot2,sravdot1,sravdot2,  ! not used now
!C     .  gsr1,gsr2,gsrav1,gsrav2,gsrdot1,gsrdot2,gsravdot1,gsravdot2, ! not used now
!C     .  pprel,gsmt)

!     module m_rdeq
!     implicit none
!     contains


      subroutine rdeq(enum1,enum2,ksop,z,v,Erseq,nmrseq,rofi,nr,nsp,a,b,l,&
           lmx,imu,scalede,avw,gmt,fmt,gmtde,fmtde,gsmt,pprel,Eout,gn,fn, &
           gdot,fdot,g2dot,f2dot)
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
!Ci   Erseq : energy output of rseq
!Ci   mnrseq: matching point passed from rseq
!Ci   avw   : Average WS radius
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
!Co   gn,fn :(2,2,nr) = (solution N,alpha,meshpoint) f and g on mesh
!Co   gdot,fdot:(2,2,nr) = (solution N,alpha,meshpoint) first energy
!Co      derivatives of f and g on mesh
!Co   g2dot,f2dot:(2,2,nr) = (solution N,alpha,meshpoint) second energy
!Co      derivatives of f and g on mesh
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
      integer j1, j2, ien, ntrial, nsh, nnod,nmrseq
      integer nnug, nnugzv, nclose
      double precision v(nr,nsp),rofi(nr),z,a,b,enum1,enum2,scalede
      double precision ksop(0:lmx,nsp,nsp,9)
      double precision pprel(4,0:lmx,2*(lmx+1),2,2)
      double precision gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),gsmt(2,2)
!C ... Local parameters
      logical check
      integer n,np,nrmx,ie,alpha,a1,ncut
!C     integer k1,k2,i,alp,ifi,fopn,ie
!C     character*(10) outs*80, fil*12
      double precision J, &   !Jacobian
       c,csq,            &   !Speed of light 274.072d0
       de,               &   !Finite energy difference for num diff.
       eav,enu1,enu2, gamma,hoc,mu,norm,pr,q,r,rm,  &  !,prodr,prodr2
       r1,rmax,smalpa(2,2),sqru,u,up,upp,w0,w1,w2,w3,x1,x2,xk,   &
       xk1,ff, e1, e2, re, En, enu1my, enu2my, E00, En00

      double precision  g1out1(nr),g2out1(nr),f1out1(nr),f2out1(nr),   &
       g1sout1(nr),g2sout1(nr),f1sout1(nr),f2sout1(nr),g1in1(nr),      &
       g2in1(nr),f1in1(nr),f2in1(nr),g1sin1(nr),g2sin1(nr),f1sin1(nr), &
       f2sin1(nr),nug1(nr),nug2(nr) ,nuf1(nr),nuf2(nr),bm(nr),         &
       vm(nr),psiout1(nr,4,2,5),psiin1(nr,4,2,5),psd(nr,4,2),Am, Aout, &
       Ain, psiout2(nr,4,2,5),psiin2(nr,4,2,5),psi(nr,4,2,5),psin(nr,4,2),&
       F(4),Pv(4),Erseq,                                               &
       g1out2(nr),g2out2(nr),f1out2(nr),f2out2(nr),                    &
       g1sout2(nr),g2sout2(nr),f1sout2(nr),f2sout2(nr),g1in2(nr),      &
       g2in2(nr),f1in2(nr),f2in2(nr),g1sin2(nr),g2sin2(nr),f1sin2(nr), &
       f2sin2(nr), g1final(nr),gfinal(nr),gfinalb(nr),ffinal(nr),      &
       ffinalb(nr),fvec(4),Pv0(4),fvec0(4),Pvpl(4),Pv02(4),fsave(4),   &
       g1finalE(nr),gfinalE(nr),gfinalbE(nr),ffinalE(nr),ffinalbE(nr), &
       g1finald(nr),gfinald(nr),gfinalbd(nr),ffinald(nr), ffinalbd(nr)

      double precision gn(2,2,nr), fn(2,2,nr),gn1(2,2,nr),&
       fn1(2,2,nr), gn2(2,2,nr), fn2(2,2,nr), Pvsh(4),gnE(2,2,nr),fnE(2,2,nr)

      double precision DeltaE, Eni(5), gp(2,2,5,nr), fp(2,2,5,nr),     &
       Eip(5), fdot(2,2,nr), gdot(2,2,nr), E0, sr, enu, avw, coeff,    &
       g2dot(2,2,nr), f2dot(2,2,nr)

      double precision gdotT(2,2), gnr(2,2), fnr(2,2), gnrT(2,2),      &
      Temp1(2,2), Temp2(2,2),  Anu(2,2), kappa(2,2), Temp3(2,2),       &
      Unit2(2,2),gammal(2,2), Cl(2,2), Wl(2,2), Dnu(2,2),One(2,2),Deltal(2,2),Fsq

      double precision tolx, tolf

      double precision Epl(1:61), Delpl, Npl

      double precision psifinal(nr,4,2),psifinalE(nr,4,2), normf, normfE

      procedure(real(8)) :: prodr,prodr2

      double precision  g1out10(nr),g2out10(nr),f1out10(nr),f2out10(nr),    &
       g1sout10(nr),g2sout10(nr),f1sout10(nr),f2sout10(nr),g1in10(nr),      &
       g2in10(nr),f1in10(nr),f2in10(nr),g1sin10(nr),g2sin10(nr),f1sin10(nr),&
       f2sin10(nr),psiout10(nr,4,2,5),psiin10(nr,4,2,5),psiout20(nr,4,2,5), &
       psiin20(nr,4,2,5),psi0(nr,4,2,5),psidot(nr,4,2),                      &
       g1out20(nr),g2out20(nr),f1out20(nr),f2out20(nr),                     &
       g1sout20(nr),g2sout20(nr),f1sout20(nr),f2sout20(nr),g1in20(nr),      &
       g2in20(nr),f1in20(nr),f2in20(nr),g1sin20(nr),g2sin20(nr),f1sin20(nr),&
       f2sin20(nr)

      double precision Pvsh2(4),E1sh, E2sh, Ash, dEsh, Efact, Vtot
      double precision Vnug(nr), Vnugzv(nr),dVnug(nr)
      double precision Eout((lmx+1),(2*lmx+2)), PercDif
      double precision Finout(4)
      double precision g1in0, f1in0, g2in0, f2in0
      double precision Pv10(4),Pv20(4),Pv100(4),Pv200(4),Pv1(4),Pv2(4)
      double precision fvec1(4), fvec2(4), fvec01(4), fvec02(4), fsave1(4), &
           fsave2(4), normsave1, normsave2, product


!C ... External calls
      external rx !, fctp0, fctp !dinv22,mul22,
!C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c


        eav = (enum1 + enum2)/2.0_8
        enu1 = -(l+1) * ksop(l,1,1,1)/2.0_8 + eav
        enu2 =   l    * ksop(l,1,1,1)/2.0_8 + eav

!c     enu1 = eav ; enu2 = eav         ! fixed enu for given l
        mu  = imu - l - 1.5d0

        if (z == 0) return

        open(1000,FILE='fvecs.dat')
        csq = c*c
        nrmx = nr
        hoc = (2*z/c)
        de = 0.03d0*scalede ;
        ff = 1

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
        e2 = enu2
        e1 = enu1

!        En = Erseq

        En = eav

        En00 = En
!!        print*,'Number of radial mesh points nr = ',nr
!!        print*,'Matching point (from rseq) nm = ',nmrseq
!!        print*,'Matching radius r(nm) = ',rofi(nmrseq)
!!        print*,' '


!C---First intergation to get P0=[g1(nr),f1(nr),g2(nr),f2(nr)]:

        g1in10 = 1e-20 !First step only gives P0, so inwards doesn't matter, we take outw(nr) as inw initial
        f1in10 = 1e-20
        g2in20 = 1e-20
        f2in20 = 1e-20

        call integration2(1, En, xk, xk1, mu, up, u,             &
           upp, sqru, vm,                                        &
           bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,      &
           g1in10,f1in10,g2in20,f2in20,                          &
           g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,  &
           f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,&
           f2sin1,psiout1, psiin1, nc, nnod,normsave1)

        call integration2(2, En, xk, xk1, mu, up, u,               &
           upp, sqru, vm,                                          &
           bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
           g1in10,f1in10,g2in20,f2in20,                            &
           g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
           f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
           f2sin2, psiout2, psiin2, nc, nnod,normsave2)

        nc = nmrseq

        g1in10 = g1out1(nr)*normsave1 !Passind outwards at nr to initial of inwards
        f1in10 = f1out1(nr)*normsave1
        g2in20 = g2out2(nr)*normsave2
        f2in20 = f2out2(nr)*normsave2

        call integration2(1, En, xk, xk1, mu, up, u,               &
           upp, sqru, vm,                                          &
           bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
           g1in10,f1in10,g2in20,f2in20,                            &
           g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,    &
           f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,  &
           f2sin1,psiout1, psiin1, nc, nnod,normsave1)

        call integration2(2, En, xk, xk1, mu, up, u,               &
           upp, sqru, vm,                                          &
           bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
           g1in10,f1in10,g2in20,f2in20,                            &
           g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
           f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
           f2sin2, psiout2, psiin2, nc, nnod,normsave2)

!        print*,'Initial gap at match F0 = ',(g1out1(nc)-g1in1(nc)), &
!       (f1out1(nc)-f1in1(nc)), (g2out2(nc)-g2in2(nc)), (f2out2(nc)-f2in2(nc))


!_________________Printing inwards, outwards before matching___________


!        ifi = 1000+100*l+imu
!        write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
!        write(ifi,*)'# l:',l,'imu:',imu
!        do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),g1out1(n),g1in1(n), &
!                   f1out1(n),f1in1(n),g2out2(n),g2in2(n),f2out2(n), &
!                   f2in2(n)
!        enddo
!        close(ifi)

!        ifi = 12000+100*l+imu
!        write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
!        write(ifi,*)'# l:',l,'imu:',imu
!        do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),g1sout1(n),g1sin1(n), &
!                   f1sout1(n),f1sin1(n),g2sout2(n),g2sin2(n),f2sout2(n), &
!                   f2sin2(n)
!        enddo
!        close(ifi)

!__________________End Print___________________________________________

!C--Building initial vector P(functions at nr by outwards integration)------

        nm = nmrseq

        Pv(1) =  g1out1(nr)*normsave1
        Pv(2) =  f1out1(nr)*normsave1
        Pv(3) =  g2out2(nr)*normsave2
        Pv(4) =  f2out2(nr)*normsave2

!!     print*,'Functions at boundary before matching Pv0 =',(Pv(1)/normsave1),&
!!             (Pv(2)/normsave1),(Pv(3)/normsave2),(Pv(4)/normsave2)
        Pv10(:) = Pv(:)
        Pv100(:) = Pv(:)

!________________________First F=0 solver_______________________
!C-- newt from Nimerical Recipes (Globaly convergent Newton).

 214    nF = 4;
        call newt(Pv,En,nF, xk, xk1, mu, up, u,                        &
             upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,     &
             lmx, imu,scalede,check,fvec)

!!        print*,' '
!!        print*,'GLOBAL NEWTON F = ',fvec


!!       print*,'Pv = ',(Pv(1)/normsave1),(Pv(2)/normsave1),(Pv(3)/normsave2),&
!!             (Pv(4)/normsave2)
        PercDif = abs(Pv10(1) - Pv(1))/abs(Pv10(1))*100
!!        print*,'rdeq l =',l,' imu =',imu,' % difference old-new g11(nr) =',PercDif, '%'

!____________________End First F=0 solver______________________


!__________Alternative F=0 solver__________________________________
!C--   mnewt from Numerical Recipes. Often finds solution too far   &
!C--   from the initial (scalar) energy value, though sometimes is  &
!C--   more effective.

!        ntrial = 100
!        tolf = 1e-12_8
!        tolx = 1e-12_8


!        call mnewt(ntrial,Pv10,En,nF,tolx,tolf, xk, xk1, mu, up, u,    &
!            upp, sqru, vm,                                             &
!            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,fvec01)


!        print*,' '
!        print*,'ALTERN METHOD F = ',fvec01

!__________End Alternative F=0 solver________________________________


!_________Third solver F=0 (broydn)__________________________________
!C-- broydn from Nimerical Recipes (Broydn metod).

!!      Pv02(:) = Pv(:)

!       call broydn(Pv100,En,nF,check, xk, xk1, mu, up, u,             &
!            upp, sqru, vm,                                            &
!            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,   &
!            fsave)


!        print*,' '
!        print*,'BROYDN F  = ',fsave

!        print*,'Is this a local minimum: ',check


!C-----Integrating with the energy found by newt: ---------------
!         Pv(:)=Pvsh(:)
!         Pv(:)=Pv100(:)
         Eout(l+1,imu) = eav
!         print*,'Eout =',Eout
         En = eav

         call integration2(1, En, xk, xk1, mu, up, u,             &
            upp, sqru, vm,                                        &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,      &
            Pv(1),Pv(2),Pv(3),Pv(4),                              &
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,  &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,&
            f2sin1,psiout1, psiin1, nc, nnod,normsave1)

         call integration2(2, En, xk, xk1, mu, up, u,             &
            upp, sqru, vm,                                          &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
            Pv(1),Pv(2),Pv(3),Pv(4),                                &
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
            f2sin2, psiout2, psiin2, nc, nnod,normsave2)


!_________________Print inwards and outwards after matching__________


!           ifi = 2000+100*l+imu
!           write(ifi,*)'r,g11o,g11i,g22o,g22i,f11o,f11i,f22o,f22i'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),g1out1(n),g1in1(n),  &
!                   g2out2(n),g2in2(n),f1out1(n),f1in1(n),f2out2(n),  &
!                   f2in2(n)
!           enddo
!           close(ifi)

!__________________End Print_________________________________________


!C---Building functions from Inwards and Outwars:--------------------

       nc=nmrseq


       gfinal(1:nc)  = g1out1(1:nc)
       gfinalb(1:nc) = g2out2(1:nc)
       ffinal(1:nc)  = f1out1(1:nc)
       ffinalb(1:nc) = f2out2(1:nc)

       gfinal((nc+1):nr) = g1in1((nc+1):nr)
       gfinalb((nc+1):nr)= g2in2((nc+1):nr)
       ffinal((nc+1):nr) = f1in1((nc+1):nr)
       ffinalb((nc+1):nr)= f2in2((nc+1):nr)

       gfinald(1:nc)  = g1sout1(1:nc)!/abs(g1sout1(nm))
       gfinalbd(1:nc) = g2sout2(1:nc)!/abs(g2sout2(nm))
       ffinald(1:nc)  = f1sout1(1:nc)!/abs(f1sout1(nm))
       ffinalbd(1:nc) = f2sout2(1:nc)!/abs(f2sout2(nm))

       gfinald((nc+1):nr) = g1sin1((nc+1):nr)!/abs(g1sin1(nm))
       gfinalbd((nc+1):nr)= g2sin2((nc+1):nr)!/abs(g2sin2(nm))
       ffinald((nc+1):nr) = f1sin1((nc+1):nr)!/abs(f1sin1(nm))
       ffinalbd((nc+1):nr)= f2sin2((nc+1):nr)!/abs(f2sin2(nm))

!_________________Print  final matched_________________________


!           ifi = 3000+100*l+imu
!           write(ifi,*)'r,gfin,ffin,gfin2,ffin2,same'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),gfinal(n),ffinal(n),     &
!                   gfinalb(n),ffinalb(n),gfinal(n),ffinal(n),gfinalb(n), &
!                   ffinalb(n)
!           enddo
!           close(ifi)

!__________________Print derivatives________________________

!           ifi = 4000+100*l+imu
!           write(ifi,*)'r,gdot11,fdot11,gdot22,fdot22,same'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),gfinald(n),ffinald(n),        &
!                   gfinalbd(n),ffinalbd(n),gfinald(n),ffinald(n),gfinalbd(n), &
!                   ffinalbd(n)
!           enddo
!           close(ifi)

!__________________End Print_________________________________________

!________________Sewing the inward/outward____________________________

!!!------ alpha-number always goes second-----------------------------

       nc = nmrseq

       gn(1,1,1:nc)=g1out1(1:nc)
       gn(2,2,1:nc)=g2out2(1:nc)
       gn(2,1,1:nc)=g2out1(1:nc)
       gn(1,2,1:nc)=g1out2(1:nc)
       fn(1,1,1:nc)=f1out1(1:nc)
       fn(2,2,1:nc)=f2out2(1:nc)
       fn(2,1,1:nc)=f2out1(1:nc)
       fn(1,2,1:nc)=f1out2(1:nc)

       gn(1,1,(nc+1):nr)=g1in1((nc+1):nr)
       gn(2,2,(nc+1):nr)=g2in2((nc+1):nr)
       gn(2,1,(nc+1):nr)=g2in1((nc+1):nr)
       gn(1,2,(nc+1):nr)=g1in2((nc+1):nr)
       fn(1,1,(nc+1):nr)=f1in1((nc+1):nr)
       fn(2,2,(nc+1):nr)=f2in2((nc+1):nr)
       fn(2,1,(nc+1):nr)=f2in1((nc+1):nr)
       fn(1,2,(nc+1):nr)=f1in2((nc+1):nr)

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

!___________________End Normalization, orthog.________________________


!_________________Print after normalization_________________________


!           ifi = 5000+100*l+imu
!           write(ifi,*)'r,gn11,fn11,gn22,fn22,etc after norm'
!           write(ifi,*)'# l:',l,'imu:',imu
!           do n=1,nr
!              write(ifi,'(9(x,g16.8))') rofi(n),gn(1,1,n),fn(1,1,n),     &
!                   gn(2,2,n),fn(2,2,n),gn(1,2,n),fn(1,2,n),gn(2,1,n),    &
!                   fn(2,1,n)
!           enddo
!           close(ifi)
!_______________________End print____________________________________

!          do j1=1,2
!             psin(:,1,j1) = gn(1,j1,:)
!             psin(:,2,j1) = fn(1,j1,:)
!             psin(:,3,j1) = gn(2,j1,:)
!             psin(:,4,j1) = fn(2,j1,:)
!          enddo


          gnr(1,1) = gn(1,1,nr); gnr(1,2) = gn(1,2,nr)
          gnr(2,1) = gn(2,1,nr); gnr(2,2) = gn(2,2,nr)
          fnr(1,1) = fn(1,1,nr); fnr(1,2) = fn(1,2,nr)
          fnr(2,1) = fn(2,1,nr); fnr(2,2) = fn(2,2,nr)

          gmt(:,:) = gnr(:,:)/rofi(nr)
          fmt(:,:) = fnr(:,:)/rofi(nr)

          if (imu == 1 .or. imu == 2*l+2) then
             nug1(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq
          else
             nug1(nr) = 1d0 + 2*z/rmax/csq + (enu1 - vm(nr) + up *bm(nr))/csq
          endif

          nug2(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq

          gsmt(1,:) = (c*nug1(nr)*psifinal(nr,2,:) - xk /rmax*psifinal(nr,1,:))
          gsmt(2,:) = (c*nug2(nr)*psifinal(nr,4,:) - xk1/rmax*psifinal(nr,3,:))

!_______________Energy derivative__________________________________

!C-------- Calculating gs, fs at 5 different energy points---------
!C-------- to find energy derivatives------------------------------


       DeltaE = 0.03d0*scalede !5e-4_8
!       E0 = Erseq
       En = eav
       E0 = En
       enu = E0

       Eni(1) = E0 - 2*DeltaE
       Eni(2) = E0 - DeltaE
       Eni(3) = E0
       Eni(4) = E0 + DeltaE
       Eni(5) = E0 + 2*DeltaE
!       Pv01(:)=Pv(:)

        do ien = 1,5
          call integration2(1, Eni(ien), xk, xk1, mu, up, u,       &
            upp, sqru, vm,                                         &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,       &
            Pv(1),Pv(2),Pv(3),Pv(4),                               &
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,   &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1, &
            f2sin1,psiout1, psiin1, nc, nnod,normsave1)
          call integration2(2, Eni(ien), xk, xk1, mu, up, u,       &
            upp, sqru, vm,                                         &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,       &
            Pv(1),Pv(2),Pv(3),Pv(4),                               &
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,   &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2, &
            f2sin2, psiout2, psiin2, nc, nnod,normsave2)

           nc=nmrseq
           Pv(1) =  g1out1(nr)*normsave1
           Pv(2) =  f1out1(nr)*normsave1
           Pv(3) =  g2out2(nr)*normsave2
           Pv(4) =  f2out2(nr)*normsave2

           nF = 4;

           call newt(Pv, Eni(ien), nF, xk, xk1, mu, up, u,         &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,    &
            lmx, imu,scalede,check,fvec)


           call integration2(1, Eni(ien), xk, xk1, mu, up, u,          &
            upp, sqru, vm,                                          &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
            Pv(1),Pv(2),Pv(3),Pv(4),                                &
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,    &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,  &
            f2sin1,psiout1, psiin1, nc, nnod,normsave1)

           call integration2(2, Eni(ien), xk, xk1, mu, up, u,          &
            upp, sqru, vm,                                          &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,        &
            Pv(1),Pv(2),Pv(3),Pv(4),                                &
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,    &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,  &
            f2sin2, psiout2, psiin2, nc, nnod,normsave2)

!_______________________Normalization of the whole___________________

           nc = nmrseq
           gfinalE(1:nc) =g1out1(1:nc)
           gfinalbE(1:nc)=g2out2(1:nc)
           ffinalE(1:nc) =f1out1(1:nc)
           ffinalbE(1:nc)=f2out2(1:nc)

           gfinalE((nc+1):nr) =g1in1((nc+1):nr)
           gfinalbE((nc+1):nr)=g2in2((nc+1):nr)
           ffinalE((nc+1):nr) =f1in1((nc+1):nr)
           ffinalbE((nc+1):nr)=f2in2((nc+1):nr)
!           psifinalE(:,1)=gfinalE(:)
!           psifinalE(:,2)=ffinalE(:)
!           psifinalE(:,3)=gfinalbE(:)
!           psifinalE(:,4)=ffinalbE(:)
!           call productn4(psifinalE(:,:),psifinalE(:,:),a,b,rofi,   &
!                nr,normfE)

!___________________End Normalization of the whole___________________

           gnE(1,1,1:nc)=g1out1(1:nc)
           gnE(2,2,1:nc)=g2out2(1:nc)
           gnE(2,1,1:nc)=g2out1(1:nc)
           gnE(1,2,1:nc)=g1out2(1:nc)
           fnE(1,1,1:nc)=f1out1(1:nc)
           fnE(2,2,1:nc)=f2out2(1:nc)
           fnE(2,1,1:nc)=f2out1(1:nc)
           fnE(1,2,1:nc)=f1out2(1:nc)

           gnE(1,1,(nc+1):nr)=g1in1((nc+1):nr)
           gnE(2,2,(nc+1):nr)=g2in2((nc+1):nr)
           gnE(2,1,(nc+1):nr)=g2in1((nc+1):nr)
           gnE(1,2,(nc+1):nr)=g1in2((nc+1):nr)
           fnE(1,1,(nc+1):nr)=f1in1((nc+1):nr)
           fnE(2,2,(nc+1):nr)=f2in2((nc+1):nr)
           fnE(2,1,(nc+1):nr)=f2in1((nc+1):nr)
           fnE(1,2,(nc+1):nr)=f1in2((nc+1):nr)
!_______________________Normalization of the whole___________________

           psifinalE(1:nr,1,1)=gnE(1,1,1:nr)
           psifinalE(1:nr,2,1)=fnE(1,1,1:nr)
           psifinalE(1:nr,3,1)=gnE(2,1,1:nr)
           psifinalE(1:nr,4,1)=fnE(2,1,1:nr)
           psifinalE(1:nr,1,2)=gnE(1,2,1:nr)
           psifinalE(1:nr,2,2)=fnE(1,2,1:nr)
           psifinalE(1:nr,3,2)=gnE(2,2,1:nr)
           psifinalE(1:nr,4,2)=fnE(2,2,1:nr)

           do j1=1,2
              call productn4(psifinalE(:,:,j1),psifinalE(:,:,j1),a,b,rofi,nr,normf)
              do j2=1,4
                 psifinalE(:,j2,j1)=psifinalE(:,j2,j1)/dsqrt(normf)
              enddo
           enddo

!--------------------Orthoganalization--------------------------------

           call productn4(psifinalE(:,:,1),psifinalE(:,:,2),a,b,rofi,nr,product)
           do j1=2,nr
              psifinalE(j1,1,2) = psifinalE(j1,1,2) - product*psifinalE(j1,1,1)
              psifinalE(j1,2,2) = psifinalE(j1,2,2) - product*psifinalE(j1,2,1)
              psifinalE(j1,3,2) = psifinalE(j1,3,2) - product*psifinalE(j1,3,1)
              psifinalE(j1,4,2) = psifinalE(j1,4,2) - product*psifinalE(j1,4,1)
           enddo

!----------------------Multiply the second solution (alpha=2) by
!---------------------- a phase factor (-1)

           if ((psifinalE(nr,1,1)*psifinalE(nr,3,2)) < 0) then
              psifinalE(:,:,2)=-psifinalE(:,:,2)
           endif


           gnE(1,1,1:nr) = psifinalE(1:nr,1,1)        ! gn(#,alpha,1:nr)
           fnE(1,1,1:nr) = psifinalE(1:nr,2,1)
           gnE(2,1,1:nr) = psifinalE(1:nr,3,1)
           fnE(2,1,1:nr) = psifinalE(1:nr,4,1)
           gnE(1,2,1:nr) = psifinalE(1:nr,1,2)
           fnE(1,2,1:nr) = psifinalE(1:nr,2,2)
           gnE(2,2,1:nr) = psifinalE(1:nr,3,2)
           fnE(2,2,1:nr) = psifinalE(1:nr,4,2)


!           do j1=1,2
!              psi(:,1,j1,ie) = gn(1,j1,:)
!              psi(:,2,j1,ie) = fn(1,j1,:)
!              psi(:,3,j1,ie) = gn(2,j1,:)
!              psi(:,4,j1,ie) = fn(2,j1,:)
!           enddo

           gp(1,1,ien,:) = gnE(1,1,:); gp(1,2,ien,:) = gnE(1,2,:)
           gp(2,1,ien,:) = gnE(2,1,:); gp(2,2,ien,:) = gnE(2,2,:)
           fp(1,1,ien,:) = fnE(1,1,:); fp(1,2,ien,:) = fnE(1,2,:)
           fp(2,1,ien,:) = fnE(2,1,:); fp(2,2,ien,:) = fnE(2,2,:)
           Eip(ien) = Eni(ien)

        enddo

!C____Energy Derivatives of g and f at nr__________________________

        do j1 = 1,2
           do j2 = 1,2

              gdot(j1,j2,:) = gp(j1,j2,1,:)*(Eip(3)-Eip(2))                   &
     *(Eip(3)-Eip(4))*                                                        &
     (Eip(3)-Eip(5))/(Eip(1)-Eip(2))/(Eip(1)-Eip(3))/(Eip(1)-Eip(4))/         &
     (Eip(1)-Eip(5)) +                                                        &
                           gp(j1,j2,2,:)*(Eip(3)-Eip(1))*                     &
      (Eip(3)-Eip(4))*                                                        &
     (Eip(3)-Eip(5))/(Eip(2)-Eip(1))/(Eip(2)-Eip(3))/(Eip(2)-Eip(4))/         &
     (Eip(2)-Eip(5))+                                                         &
                           gp(j1,j2,3,:)*((Eip(3)-Eip(2))*                    &
     (Eip(3)-Eip(4))*                                                         &
     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(4))*                       &
     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
     (Eip(3)-Eip(4)))                                                         &
     /(Eip(3)-Eip(1))/(Eip(3)-Eip(2))/(Eip(3)-Eip(4))/                        &
     (Eip(3)-Eip(5))+                                                         &
                           gp(j1,j2,4,:)*(Eip(3)-Eip(1))                      &
      *(Eip(3)-Eip(2))*                                                       &
     (Eip(3)-Eip(5))/(Eip(4)-Eip(1))/(Eip(4)-Eip(2))/(Eip(4)-Eip(3))/         &
     (Eip(4)-Eip(5))+                                                         &
                           gp(j1,j2,5,:)*(Eip(3)-Eip(1))                      &
      *(Eip(3)-Eip(2))*                                                       &
     (Eip(3)-Eip(4))/(Eip(5)-Eip(1))/(Eip(5)-Eip(2))/(Eip(5)-Eip(3))/         &
     (Eip(5)-Eip(4))


             fdot(j1,j2,:) = fp(j1,j2,1,:)*(Eip(3)-Eip(2))                    &
      *(Eip(3)-Eip(4))*                                                       &
     (Eip(3)-Eip(5))/(Eip(1)-Eip(2))/(Eip(1)-Eip(3))/(Eip(1)-Eip(4))/         &
     (Eip(1)-Eip(5)) +                                                        &
                           fp(j1,j2,2,:)*(Eip(3)-Eip(1))                      &
      *(Eip(3)-Eip(4))*                                                       &
     (Eip(3)-Eip(5))/(Eip(2)-Eip(1))/(Eip(2)-Eip(3))/(Eip(2)-Eip(4))/         &
     (Eip(2)-Eip(5))+                                                         &
                           fp(j1,j2,3,:)*((Eip(3)-Eip(2))                     &
      *(Eip(3)-Eip(4))*                                                       &
     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(4))*                       &
     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
     (Eip(3)-Eip(5)) + (Eip(3)-Eip(1))*(Eip(3)-Eip(2))*                       &
     (Eip(3)-Eip(4)))                                                         &
     /(Eip(3)-Eip(1))/(Eip(3)-Eip(2))/(Eip(3)-Eip(4))/                        &
     (Eip(3)-Eip(5))+                                                         &
                           fp(j1,j2,4,:)*(Eip(3)-Eip(1))                      &
      *(Eip(3)-Eip(2))*                                                       &
     (Eip(3)-Eip(5))/(Eip(4)-Eip(1))/(Eip(4)-Eip(2))/(Eip(4)-Eip(3))/         &
     (Eip(4)-Eip(5))+                                                         &
                           fp(j1,j2,5,:)*(Eip(3)-Eip(1))                      &
     *(Eip(3)-Eip(2))*                                                        &
     (Eip(3)-Eip(4))/(Eip(5)-Eip(1))/(Eip(5)-Eip(2))/(Eip(5)-Eip(3))/         &
     (Eip(5)-Eip(4))

             g2dot(j1,j2,:) = (-gp(j1,j2,1,:) + 16*gp(j1,j2,2,:) -            &
              30*gp(j1,j2,3,:) + 16*gp(j1,j2,4,:) - gp(j1,j2,5,:))/12/DeltaE**2

             f2dot(j1,j2,:) = (-fp(j1,j2,1,:) + 16*fp(j1,j2,2,:) -            &
              30*fp(j1,j2,3,:) + 16*fp(j1,j2,4,:) - fp(j1,j2,5,:))/12/DeltaE**2

           enddo
        enddo

!C__________ Orthogonalize phidot to phi_____________________________________

         psidot(1:nr,1,1)=gdot(1,1,1:nr)
         psidot(1:nr,2,1)=fdot(1,1,1:nr)
         psidot(1:nr,3,1)=gdot(2,1,1:nr)
         psidot(1:nr,4,1)=fdot(2,1,1:nr)
         psidot(1:nr,1,2)=gdot(1,2,1:nr)
         psidot(1:nr,2,2)=fdot(1,2,1:nr)
         psidot(1:nr,3,2)=gdot(2,2,1:nr)
         psidot(1:nr,4,2)=fdot(2,2,1:nr)

         do j2=1,2
           call productn4(psifinal(:,:,j2),psidot(:,:,j2),a,b,rofi,nr,product)
           do j1=2,nr
              gdot(1,j2,j1) = gdot(1,j2,j1) - product*gn(1,j2,j1)
              fdot(1,j2,j1) = fdot(1,j2,j1) - product*fn(1,j2,j1)
              gdot(2,j2,j1) = gdot(2,j2,j1) - product*gn(2,j2,j1)
              fdot(2,j2,j1) = fdot(2,j2,j1) - product*fn(2,j2,j1)
           enddo
         enddo

         call productn4(psifinal(:,:,1),psidot(:,:,2),a,b,rofi,nr,product)
         do j1=2,nr
           gdot(1,2,j1) = gdot(1,2,j1) - product*gn(1,1,j1)
           fdot(1,2,j1) = fdot(1,2,j1) - product*fn(1,1,j1)
           gdot(2,2,j1) = gdot(2,2,j1) - product*gn(2,1,j1)
           fdot(2,2,j1) = fdot(2,2,j1) - product*fn(2,1,j1)
         enddo
         call productn4(psifinal(:,:,2),psidot(:,:,1),a,b,rofi,nr,product)
         do j1=2,nr
           gdot(1,1,j1) = gdot(1,1,j1) - product*gn(1,2,j1)
           fdot(1,1,j1) = fdot(1,1,j1) - product*fn(1,2,j1)
           gdot(2,1,j1) = gdot(2,1,j1) - product*gn(2,2,j1)
           fdot(2,1,j1) = fdot(2,1,j1) - product*fn(2,2,j1)
         enddo

        gmtde(:,:) = gdot(:,:,nr)/rofi(nr)
        fmtde(:,:) = fdot(:,:,nr)/rofi(nr)

!_______________Energy derivative end_______________________________

!__________________Potemtial parameters_____________________________



        gdotT(1,1) = gdot(1,1,nr); gdotT(2,2) = gdot(2,2,nr) ! Deriv g transposed  nr
        gdotT(2,1) = gdot(1,2,nr); gdotT(1,2) = gdot(2,1,nr)

        gnrT(1,1) = gnr(1,1); gnrT(2,2) = gnr(2,2) ! g transposed at nr
        gnrT(2,1) = gnr(1,2); gnrT(1,2) = gnr(2,1)

!C--- Parameters are made by fdpp routine (all but smalpa (pprel(4,...)))

        call mul22(gnr, gdotT, Temp1)            ! Temp1 = g*gdot^T
        sr = rofi(nr)
        call matrconst(Temp1,sr,Temp2)           ! Temp2 = g*gdot^T * s
        call inver(Temp2, Temp1)                 ! Temp1 =(g*gdot^T*s)^-1

        Anu(:,:) = -Temp1(:,:)

        call inver(gnr, Temp1)                   ! Temp1 = g^-1
        call mul22(fnr, Temp1, Temp2)            ! Temp2 = f(nr)*g^-1
        call matrconst(Temp2, sr, Temp1)         ! Temp1 = s*f(nr)*g^-1

        Unit2(1,2) = 0.0_8; Unit2(2,1) = 0.0_8
        Unit2(1,1) = 1.0_8; Unit2(2,2) = 1.0_8

        One(1,1) = 1.0_8; One(1,2) = 1.0_8
        One(2,1) = 1.0_8; One(2,2) = 1.0_8

        Dnu(:,:) = Temp1(:,:) - kappa(:,:) - Unit2(:,:)


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

!        Deltal(1,1) = Wl(1,1)*Wl(1,1)+Wl(2,1)*Wl(2,1)
!        Deltal(1,2) = Wl(1,1)*Wl(1,2)+Wl(2,1)*Wl(2,2)
!        Deltal(1,1) = Wl(1,2)*Wl(1,1)+Wl(2,1)*Wl(2,2)
!        Deltal(1,1) = Wl(1,2)*Wl(1,2)+Wl(2,2)*Wl(2,2)

!        print*,'Parameters : '
!        print*,'Gamma = '
!        print*,gammal
!        print*,'W, C = ',Wl,Cl


!Co      pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),1:2,1:2)

!        pprel(1,l,imu,1:2,1:2) = Cl(:,:)
!        pprel(2,l,imu,1:2,1:2) = gammal(:,:)
!        pprel(3,l,imu,1:2,1:2) = Deltal(:,:)

!___________________Calculaing small p=<gdot|gdot>______________________

!        gdot1(:,:) = gdot(1,j2,:)
!        gdot2(:,:) = gdot(2,j2,:)

        do j1=1,2
           psd(:,1,j1) = gdot(1,j1,:)
           psd(:,2,j1) = fdot(1,j1,:)
           psd(:,3,j1) = gdot(2,j1,:)
           psd(:,4,j1) = fdot(2,j1,:)
        enddo

!C --- Small parameter p=<gdot|gdot> by trapezium rule ---
       smalpa = 0.0_8
       nrmx = nr
       do  n = 2, nrmx-1
         np = n + 1 ; r = rofi(n) ; r1 = rofi(np)
         w0 = (1 + xk *(xk+1) /((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
         w1 = (1 + xk *(xk+1) /((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
         w2 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
         w3 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
         x1 = a*(r+b) ; x2 = a*(r1+b)
         do  a1 = 1, 2
           smalpa(a1,:) = smalpa(a1,:)                                   &
            + x1*(w0*psd(n ,1,a1)*psd(n ,1,:)+psd(n ,2,a1)*psd(n ,2,:))  &
            + x2*(w1*psd(np,1,a1)*psd(np,1,:)+psd(np,2,a1)*psd(np,2,:))  &
            + x1*(w2*psd(n ,3,a1)*psd(n ,3,:)+psd(n ,4,a1)*psd(n ,4,:))  &
            + x2*(w3*psd(np,3,a1)*psd(np,3,:)+psd(np,4,a1)*psd(np,4,:))
         enddo
         smalpa = smalpa/2.0_8
       enddo

       pprel(4,l,imu,:,:) = smalpa(:,:)


!________________________Small p end_______________________________

!__________________Potential parameters end_________________________




 9999 continue


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
      close(1)

      end subroutine rdeq


      subroutine fdpp(enu1,enu2,ksop,shft,gmt,fmt,gmtde,z,rmax,avw,l, lmx,imu,pprel)
!C- Fully relativistic potential parameters
!C------------------------------------------------------------------------------
!Ci Inputs
!Ci   enu1 :linearization energy, spin 1 (called enum1 in rdeq.f)
!Ci   enu2 :linearization energy, spin 2 (called enum2 in rdeq.f)
!Ci   ksop  :spin orbit coupling parameters, used to adjust enu1,enu2
!Ci   shft  ?
!Ci   gmt   :2x2 normalized wave function g x r, at the MT surface
!Ci   fmt   :2x2 normalized wave function f x r, at the MT surface
!Ci   gmtde :gmtde :(2x2 matrix) gdot at rmax
!Ci   z     :nuclear charge
!Ci   rmax  :augmentation radius, in a.u.
!Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
!Ci         :used in the definition of potential parameters
!Ci   l     :l quantum number
!Ci   lmx   :dimensions ksop and pprel
!Ci   imu   :l + mu + 3/2, mu = quantum number
!Ci         :imu has range 1:2(l+1) and  mu has range (-l-1/2:l+1/2)
!Co Outputs
!Co   pprel :Relativistic potential parameters in kappa-mu representation
!Co         :pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),1:2,1:2)
!Co         :Made here:
!Co         :pprel(1,:,:,1:2,1:2) = C
!Co         :pprel(2,:,:,1:2,1:2) = gamma
!Co         :pprel(3,:,:,1:2,1:2) = delta
!Co         :Not made here:
!Co         :pprel(4,:,:,1:2,1:2) = small parameter p
!Cs Command-line switches
!Cl Local variables
!Cl         :
!Cr Remarks
!Cr
!Cu Updates
!Cu   18 Oct 13 (Belashchenko) cleanup
!C ----------------------------------------------------------------------
      implicit none
!C ... Passed parameters
      integer l,lmx,imu
      double precision gmt(2,2), fmt(2,2), gmtde(2,2), mu, rmax, avw, z, enu1, enu2, shft
!C      double precision sr1, sr2, srav1, fmtde(2,2), gsmt(2,2)
!C     .                 srav2, srdot1, srdot2, sravdot1, sravdot2,
!C     .                 gsrav1, gsrav2, gsrdot1, gsrdot2, gsravdot1,
!C     .                 gsravdot2, gsr1, gsr2
!C     Local variables
!C     integer  i, j
!C     double precision ftmt(2,2),fmtdet(2,2),temp(2,2)
      double precision pprel(4,0:lmx,2*(lmx+1),2,2)

!C ... Local parameters
      double precision gtmt(2,2), gmt1(2,2),                       &
                      fg(2,2),                                     &
                      xk, xk1, kapa(2,2), D(2,2),                  &
                      delta1(2,2),unit(2,2), q(2,2), q1(2,2),      &
                      gam(2,2), gam1(2,2), aa(2,2), aa1(2,2),      &
                      one(2,2), two(2,2), gmtdet(2,2), gamt(2,2),  &
                      yv(2,2), yv1(2,2), vm(2,2), enum(2,2),       &
                      gamma(2,2), gammat(2,2), delt(2,2),          &
                      deltt(2,2), cm(2,2), cm1(2,2),               &
                      savw, ksop(0:lmx,2,2,6),                     &
                      xx(2,2),c
!C     double precision clebsh1(2,2), clebsh1t(2,2),clebsh(2,2), u,
!C     double precision Ddot(2,2) Dtdot(2,2), gmtde1(2,2),fgde(2,2),delta(2,2),deltat(2,2),deltat1(2,2),fg1(2,2)

!C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      if (z == 0) goto 210

!C --- Setup ---
      mu = dble(imu-l) - 1.5d0
!c     cc = c
!c     cc = 274.072d0
!c     cc = 274.072d0
      savw = (avw/rmax)**(2*l+1)

!C --- Test: Eliminate Spin-Orbit for non-mag. materials
!C     mu = dble(l)+0.5d0

      if (imu == 1 .or. imu == 2*l+2) then
        xk1 = dble(-l) ; xk = xk1
      else
        xk1 = dble(-l) ; xk = dble(l+1)
      endif

      unit(1,1) = 1.0_8 ; unit(2,2) = 1.0_8 ;
      unit(1,2) = 0.0_8 ; unit(2,1) = 0.0_8

      enum = 0.0_8
      enum(1,1) = -(l+1) * ksop(l,1,1,1)/2.0_8 + (enu1 + enu2)/2.0_8
      enum(2,2) =   l    * ksop(l,1,1,1)/2.0_8 + (enu1 + enu2)/2.0_8

!c     enum(1,1) = (enu1 + enu2)/2 ; enum(2,2) = (enu1 + enu2)/2

      if (imu == 1 .or. imu == 2*l+2) enum(1,1) = enum(2,2)

      kapa = 0.0_8 ; kapa(1,1) = xk ; kapa(2,2) = xk1

!C --- Invert g and g_dot----------------------------
!C      The g and f do not commute. In Phys Rev. B
!C      43, 14414, 1991 it is D = c*S*g^(-1)*f - k - I
!C      In Sov. Phys. Solid State 31(8), 1285 1989 is
!C      D = c*S*f*g^(-1) - k - I. However only the second
!C      formula gives symmetric D. The order will depend
!C      on the order in g and f. If alpha numerates columns then
!C      the second formula works, if alpha numerates strokes then
!C      the first formula should be applied.
!C----------------------------------------------------
      call dinv22(gmt,gmt1) ; fg = matmul(fmt,gmt1)
!c     call dinv22(gmtde,gmtde1)
!c     fgde = matmul(fmtde,gmtde1)

!c     call mul22(gsmt,gmt1,fg1)

!C --- Logarithmic derivative matrices-----------------
      D(:,:) = rmax*c*fg(:,:) - kapa(:,:)
!c     D(:,:) = fg1(:,:) - unit(:,:)
!c     Ddot(:,:) = rmax*c*fgde(:,:) - kapa(:,:)

!C --- Symmetrize D (it is symmetric, but for numerical)
      D(1,2) = (D(1,2) + D(2,1))/2 ; D(2,1) = D(1,2)

!C --- Make corrections for diagonal D, D_dot--------

!cc       if (mu < 0) then
!cc       D(1,1) = D(1,1) - srav1 + sr1
!cc       D(2,2) = D(2,2) - srav1 + sr1
!cc       Ddot(1,1) = Ddot(1,1) - sravdot1 + srdot1
!cc       Ddot(2,2) = Ddot(2,2) - sravdot1 + srdot1
!cc       else
!cc       D(1,1) = D(1,1) - srav2 + sr2
!cc       D(2,2) = D(2,2) - srav2 + sr2
!cc       Ddot(1,1) = Ddot(1,1) - sravdot2 + srdot2
!cc       Ddot(2,2) = Ddot(2,2) - sravdot2 + srdot2
!cc       endif

!C------Make corrections for diagonal g, g_dot--------

!cc       if (mu < 0) then
!cc       gmt(1,1) = gmt(1,1) - gsrav1 + gsr1
!cc       gmt(2,2) = gmt(2,2) - gsrav1 + gsr1
!cc       gmtde(1,1) = gmtde(1,1) - gsravdot1 + gsrdot1
!cc       gmtde(2,2) = gmtde(2,2) - gsravdot1 + gsrdot1
!cc       else
!cc       gmt(1,1) = gmt(1,1) - gsrav2 + gsr2
!cc       gmt(2,2) = gmt(2,2) - gsrav2 + gsr2
!cc       gmtde(1,1) = gmtde(1,1) - gsravdot2 + gsrdot2
!cc       gmtde(2,2) = gmtde(2,2) - gsravdot2 + gsrdot2
!cc      endif

!C--Transpose gmt and Ddot -------------------------
      gtmt   = transpose(gmt) ; gmtdet = transpose(gmtde)
!c     Dtdot = transpose(Ddot)

!C --- Make unscreened potential parameters-------------
!c     delta(:,:) = D(:,:) - Ddot(:,:)
!c     deltat(:,:) = D(:,:) - Dtdot(:,:)

!c     call dinv22(delta, delta1)
!c     call dinv22(deltat, deltat1)

!cc       call mul22(gmtdet,gmt,one)
!cc       call mul22(gtmt,gmtde,two)

      call mul22(gmt,gmtdet,one) ! one = g x gdot^T
      call mul22(gmtde,gtmt,two) ! two = gdot x g^T

!c     aa(:,:) = -(delta1(:,:)+deltat1(:,:))/2  ! Eq. (19) in Solovyev
      aa(:,:) = -rmax*(one(:,:) + two(:,:))/2  ! This symmetrization is superfluous
      call dinv22(aa,aa1)                      ! aa1 : A_L in Shick

      xx(:,:) = D(:,:) - l * unit(:,:)
      call dinv22(xx,delta1)                   !delta1 = (D-l)^-1
      q(:,:) = xx(:,:) + aa1(:,:)              !q = D + A - l
      call dinv22(q,q1)                        !q1 = (D + A -l)^-1
      q(:,:) = 2*(2*l+1) * savw * ((2*l+1)*q1(:,:)+unit(:,:))
                                     ! 2(2l+1) savw (D+A+l+1)(D+A-l)^-1
                                     ! q = gamma_L^-1 (inverse gamma from Shick)

      xx = matmul(aa,xx) + unit      ! xx = A^-1 * (D + A - l)

      call dinv22(xx,gam1)           ! gam1 = (D + A - l)^-1 * A
!cc       call mul22(gtmt,gam1,gam)
      call mul22(gam1,gmt,xx)        ! xx = (D + A - l)^-1 * A * g

      gam(:,:) = dsqrt(2d0*rmax*savw)*(2*l+1)*xx(:,:)
      gamt = transpose(gam)
      yv = delta1 + aa        !yv = (D-l)^-1 + A^-1

      call dinv22(yv,yv1)     !yv1 = ((D-l)^-1+A^-1)^-1
      call mul22(gtmt,yv1,yv) !yv = g^T * (..)
      call mul22(yv,gmt,vm)   !vm = g^T * A [D + A - l]^-1 (D-l) * g

      vm(:,:) = rmax*vm(:,:) + enum(:,:)

!C --- Potential parameters in the nearly orthonormalized-----
!C     representation
!C------------------------------------------------------------

      call dinv22(q,gamma)  ! this is gamma_L from Shick

      gammat = transpose(gamma) ! gamma should be symmetric anyway

!cc       call mul22(gam,gamma,delt)
!cc       call mul22(gammat,gamt,deltt)

      call mul22(gamma,gam,delt) ! delt = (..) * g
      call mul22(gamt,gamma,deltt) ! deltt = g^T * (..)

!cc      call mul22(gam,gamma,cm)
!cc     call mul22(cm,gamt,cm1)

!cc       call mul22(delt,q,cm)
!cc       call mul22(cm,deltt,cm1)

      call mul22(deltt,q,cm)  ! cm  = g^T * gam_L^-1
      call mul22(cm,delt,cm1) ! cm1 = g^T * (gam_L)^-1 * g

      cm(:,:) = vm(:,:) + cm1(:,:) + shft * unit(:,:)

      pprel(1,l,imu,:,:) = cm(:,:)
      pprel(2,l,imu,:,:) = gamma(:,:)
      pprel(3,l,imu,:,:) = delt(:,:)

!C       open(190,file='cm.dat',status='unknown')
!C       open(191,file='gamma.dat',status='unknown')
!C       open(192,file='delt.dat',status='unknown')

!c     do  i = 1, 2
!c       do  j = 1, 2
!C          write(190,*) cm(i,j)
!C          write(191,*) gamma(i,j)
!C          write(192,*) delt(i,j)
!c       enddo
!c     enddo

!C --- Rotate cm---------------------------------------

!cc       u = mu/(dble(l)+0.5d0)

!cc       clebsh(1,1) = dsqrt(1d0+u)/dsqrt(2d0)
!cc       clebsh(2,2) = clebsh(1,1)
!cc       clebsh(2,1) = dsqrt(1d0-u)/dsqrt(2d0)
!cc       clebsh(1,2) = -clebsh(2,1)

!cc       call dinv22(clebsh,clebsh1)

!cc      call mul22(clebsh1t,cm,cm1)
!cc      call mul22(cm1,clebsh1,cm)

!cc      if (dabs(mu) == (l+0.5d0)) then
!cc      write(77,*) mu, xk-1, xk1-1
!cc      write(77,*) '--------------------------------------------'
!cc      write(77,*) cm(1,1), cm(1,2)
!cc      write(77,*) cm(2,1), cm(2,2)
!cc      write(77,*) '--------------------------------------------'
!cc      endif


!C-------This is a failed way to built the pot. pars-----------
!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             delta(i,j) = D(i,j) - Ddot(i,j)
!cc          enddo
!cc       enddo

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             deltat(i,j) = delta(j,i)
!cc          enddo
!cc       enddo


!cc       call dinv22(delta, delta1)
!cc       call dinv22(deltat, deltat1)

!cc       call mul22(gmtdet,gmt,one)
!cc       call mul22(gtmt,gmtde,two)

!cc       do  i = 1, 2
!cc        do  j = 1, 2
!cc        aa(i,j) = -0.5d0*(delta1(i,j) + deltat1(i,j))
!cc         aa(i,j) = -0.5d0*rmax*(one(i,j) + two(i,j))
!cc        enddo
!cc       enddo

!cc       call dinv22(aa,aa1)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             delta(i,j) = D(i,j) - lmatrix(i,j)
!cc           enddo
!cc      enddo

!cc       call dinv22(delta,delta1)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             q(i,j) = delta(i,j) + aa1(i,j)
!cc          enddo
!cc       enddo

!cc      write(77,*) '--------------------------------------------'
!cc      write(77,*) z, l, mu, xk, xk1
!cc      write(77,*) '--------------------------------------------'
!cc      write(77,*) delta(1,1), delta(1,2)
!cc      write(77,*) delta(2,1), delta(2,2)
!cc      write(77,*) '--------------------------------------------'
!cc      write(77,*)  aa1(1,1), aa1(1,2)
!cc      write(77,*)  aa1(2,1), aa1(2,2)

!cc       call dinv22(q,q1)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             q(i,j) = 2*(2*l+1)*savw*((2*l+1)*q1(i,j)+unit(i,j))
!cc          enddo
!cc       enddo

!cc       if (dabs(mu) == (l+0.5d0)) then
!cc          write(77,*) l,q(1,1)
!cc       endif

!cc       call mul22(aa,delta,gam)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             gam(i,j) = gam(i,j) + unit(i,j)
!cc          enddo
!cc       enddo

!cc       call dinv22(gam,gam1)
!cc       call mul22(gam1,gmt,gam)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             gam(i,j) = dsqrt(2d0*rmax*savw)*(2*l+1)*gam(i,j)
!cc          enddo
!cc       enddo

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             gamt(i,j) = gam(j,i)
!cc          enddo
!cc       enddo

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc          yv(i,j) = delta1(i,j) + aa(i,j)
!cc          enddo
!cc       enddo

!cc       call dinv22(yv,yv1)
!cc       call mul22(gtmt,yv1,yv)
!cc       call mul22(yv,gmt,vm)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc            vm(i,j) = rmax*vm(i,j) + enum(i,j)
!cc          enddo
!cc       enddo

!C --- Potential parameters in the nearly orthonormalized representation
!C------------------------------------------------------------

!cc       call dinv22(q,gamma)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc             gammat(i,j) = gamma(j,i)
!cc          enddo
!cc       enddo

!cc       call mul22(gam,gamma,delt)

!cc       call mul22(gam,gamma,cm)
!cc       call mul22(cm,gamt,cm1)

!cc       do  i = 1, 2
!cc          do  j = 1, 2
!cc            cm(i,j) = vm(i,j) + cm1(i,j)
!cc          enddo
!cc       enddo

!cc       open(190,file='cm.dat',status='unknown')
!cc       open(191,file='gamma.dat',status='unknown')
!cc       open(192,file='delt.dat',status='unknown')

!cc       do  i = 1, 2
!cc        do  j = 1, 2
!cc           pprel(1,l,imu,i,j) = cm(i,j)
!cc           write(190,*) cm(i,j)
!cc           pprel(2,l,imu,i,j) = gamma(i,j)
!cc           write(191,*) gamma(i,j)
!cc           pprel(3,l,imu,i,j) = delt(i,j)
!cc           write(192,*) delt(i,j)
!cc        enddo
!cc       enddo

!C-------End of the failed way--------------------------------------

  210 return

      end subroutine fdpp

      subroutine mul22(a,b,c)
!C-    Multiplies 2x2 matrices a*b
!C     a is from the left be carefull for non commuting matrices
      implicit none
      double precision a(2,2), b(2,2), c(2,2)

      c(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
      c(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
      c(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
      c(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

      end subroutine mul22

      subroutine productn4(w1,w2,a,b,rofi,nr,normn)
!C-----Integrats wavefunction for normalization purposes-------
      integer nr,n,i
      double precision a,b,rofi(nr), normn
      double precision w1(nr,4), w2(nr,4)

      normn = 0.0_8
      do i = 1, 4
        do  n = 2, nr-1
          f1 = a*(rofi(n  ) + b)/2.0_8
          f2 = a*(rofi(n+1) + b)/2.0_8
          normn = normn + f1*w1(n,i)*w2(n,i) + f2*w1(n+1,i)*w2(n+1,i)
        enddo
      enddo
      end subroutine productn4

      subroutine productn5(w1,w2,a,b,rofi,n0,nr,normn)
!C-----Integrats wavefunction for normalization purposes-------
      integer n0, nr, n, i
      double precision a,b,rofi(nr), normn
      double precision w1(nr,4), w2(nr,4)

      if ((n0==0).or.(n0==1)) n0=2
      normn = 0.0_8
      do i = 1, 4
        do  n = n0, nr-1
          f1 = a*(rofi(n  ) + b)/2.0_8
          f2 = a*(rofi(n+1) + b)/2.0_8
          normn = normn + f1*w1(n,i)*w2(n,i) + f2*w1(n+1,i)*w2(n+1,i)
        enddo
      enddo
      end subroutine productn5



      double precision function prodr(mode,w1,w2,a,b,rofi,nr)
!C-----Integrats wavefunction for normalization purposes-------
      implicit none
      integer mode,nr
      double precision w1(nr,*),w2(nr,*),a,b,rofi(nr)
      integer i,imax,n
      double precision f1,f2

      if (mode == 0) then
        imax = 4
      elseif (mode == 1) then
        imax = 1
      else
        call rx('prodr: unknown mode')
      endif

      prodr = 0.0_8
      do i = 1, imax
        do  n = 2, nr-1
          f1 = a * (rofi(n  ) + b)/2.0_8
          f2 = a * (rofi(n+1) + b)/2.0_8
          prodr = prodr + f1*w1(n,i)*w2(n,i) + f2*w1(n+1,i)*w2(n+1,i)
        enddo
      enddo
      end function prodr

      double precision function prodr2(w1,w2,a,b,rofi,nr)
      integer mode,nr
      double precision w1(nr),w2(nr),a,b,rofi(nr)
      integer i,imax,n
      double precision f1,f2

        do  n = 2, nr-1
          f1 = a * (rofi(n)  + b)/2 ; f2 = a * (rofi(n+1) + b)/2.0_8
          prodr2 = prodr2 + f1*w1(n)*w2(n) + f2*w1(n+1)*w2(n+1)
        enddo

      end function prodr2


      subroutine matrconst(a,b,c)
!C------Multiplies 2x2 matrix a  by number b
      implicit none
      double precision a(2,2), c(2,2)
      double precision b
        c(1,1) = a(1,1)*b
        c(1,2) = a(1,2)*b
        c(2,1) = a(2,1)*b
        c(2,2) = a(2,2)*b
      end subroutine matrconst


      subroutine inver(a,c)
!C-------Inverses 2x2 matrix
      implicit none
      double precision a(2,2), c(2,2)
      double precision Det

        Det = a(1,1)*a(2,2) - a(2,1)*a(1,2)
        c(1,1) = a(2,2)/Det
        c(2,2) = a(1,1)/Det
        c(2,1) = -a(2,1)/Det
        c(1,2) = -a(1,2)/Det

      end subroutine inver

!AV----------------------- Integration inward-outward_______________

      subroutine integration2(alpha, En, xk, xk1, mu, up, u,           &
          upp, sqru, vm,                                               &
          bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,             &
          g1in0,f1in0,g2in0,f2in0,                                     &
          g1out,f1out,g2out,f2out,g1sout,f1sout,g2sout,f2sout,         &
          g1in,f1in,g2in,f2in,g1sin,f1sin,g2sin,f2sin, psiout, psiin,  &
          nc,nnod,normsave)
!C----------Integrates inwards and outwards-----------------------------
      implicit none
!C ... Passed parameters
      integer nr,nsp,imu,l,lmx,ie, jalp,j1,j2,nfr,ncsave,nnod, ncut
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
!C ... Local parameters
      integer n,np,alpha,nrk,it,nrmx, n0, nprint, nit, nm, nc
!C     integer k1,k2,i,alp,ifi,fopn, ie
!C     character*(10) outs*80, fil*12
      double precision J, &   !Jacobian
       c,csq,                                                         & !Speed of light 274.072d0 !     .  de,                  !Finite energy difference for num diff.
       bmint,df11,df12,df13,df14,df21,df22,df23,df24,dg11,dg12,dg13,  &
       dg14,dg21,dg22,dg23,dg24,dx,e1,e2,En,f1c,f1p,f1sn,             &
       f1snf1,f1snm1,f1snm2,f1snm3,f1snm4,f2c,f2p,f2sn,f2snf2,f2snm1, &
       f2snm2,f2snm3,f2snm4,ff,g1c,g1p,g1sn,g1sng1,g1snm1,g1snm2,     &
       g1snm3,g1snm4,g2c,g2p,g2sn,g2sng2,g2snm1,g2snm2,g2snm3,g2snm4, &
       gamma,hoc,mu,norm,nuf1int,nuf2int,nug1int,nug2int,pr,prodr,q,r,&
       r1,rmax,smalpa(2,2),sqru,u,up,upp,vmint,w0,w1,w2,w3,x1,x2,xk,  &
       xk1, eps1, eps2, eps3, eps4, Delt, r0, mu1, q01, mu2, q02,rm

      double precision g1out(nr),g2out(nr),f1out(nr),f2out(nr),       &
       g1sout(nr),g2sout(nr),f1sout(nr),f2sout(nr),nug1(nr),          &
       nug2(nr),nuf1(nr),nuf2(nr),bm(nr),vm(nr),psiout(nr,4,2,5),     &
       psiin(nr,4,2,5),psd(nr,4,2), g1in(nr),g2in(nr),                &
       f1in(nr),f2in(nr),g1sin(nr),g2sin(nr),                         &
       f1sin(nr),f2sin(nr), psiinn(nr,4), psioutn(nr,4)
      double precision fr, dg1new, dg1old, g01,g02,normn,norm1,       &
       norm2

      double precision lamb1, lamb2, sq11, sq12, sq21, sq22, khi,     &
                       rin, sigma1, sigma2, b11, b12, a11, a12,       &
                       g01i, g02i, f01i, f02i
      double precision g1inmax, g1outmax
      double precision g1in0,f1in0,g2in0,f2in0,normsave

!C ... External calls
!C     external dinv22,mul22,rx
!C     Speed of light, or infinity in nonrelativistic case

      common /cc/ c


!______________________Zeros_________

         g1out(:)=0d0; f1out(:)=0d0; g2out(:)=0d0; f2out(:)=0d0
         g1sout(:)=0d0; f1sout(:)=0; g2sout(:)=0; f2sout(:)=0
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

      e1 = -(l+1) * ksop(l,1,1,1)/2.0_8 + En
      e2 =   l    * ksop(l,1,1,1)/2.0_8 + En

!      e1 = En
!      e2 = En

!      e2 = En
!      e1 = En - (2*l+1)* ksop(l,1,1,1)/2


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


!      nrmx = nr
!$      nrk = 6 ; nrmx = nr
      nrk = nr-1 ; nrmx = nr
      g1outmax = 0
      g1inmax = 0
      g1out(1) = 0.0_8 ; f1out(1) = 0.0_8
      g2out(1) = 0.0_8 ; f2out(1) = 0.0_8

!C --- Initial conditions at r-->0: V = V0 - 2Z/r, B = const
!C     The wave function is a polynomial g0i = r^(gi-1)*Sum_over_v(g1v*r^v)

      if (alpha == 1) then
        gamma = dsqrt(xk*xk - hoc*hoc)
        g1out(2) = 1d0 ; f1out(2) = (xk + gamma)/hoc
        g2out(2) = 0.0_8 ; f2out(2) = 0.0_8
      else
        gamma = dsqrt(xk1*xk1 - hoc*hoc)
        g1out(2) = 0.0_8 ; f1out(2) = 0.0_8
        g2out(2) = 1d0 ; f2out(2) = (xk1 + gamma)/hoc
      endif

!      e1 = En;
!      e2 = En;
!      e2 = enu2 + (ie - 3)*de ; e1 = enu1 + (ie - 3)*de


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


        r = rofi(nr)
        if (alpha == 1) then
          g1in(nr) = g1in0
          f1in(nr) = f1in0
          g2in(nr) = 0
          f2in(nr) = 0
        else
          g1in(nr) = 0
          f1in(nr) = 0
          g2in(nr) = g2in0
          f2in(nr) = f2in0
        endif

        J = a*(rofi(nr) + b) ; q = dexp(-a/2.0_8) ; r = rofi(nr)


        g1sin(nr) = J*(c*nug1(nr)*f1in(nr) - xk/r*g1in(nr))
        f1sin(nr) = J*( (nuf1(nr)*g1in(nr) - sqru*bm(nr)*g2in(nr))/c + xk/r*f1in(nr))
        g2sin(nr) = J*(c*nug2(nr)*f2in(nr) - xk1/r*g2in(nr))
        f2sin(nr) = J*( (nuf2(nr)*g2in(nr) - sqru*bm(nr)*g1in(nr))/c + xk1/r*f2in(nr))

!C   --- Runge-Kutta for the first NKR points ---
!$        do  n = nr, nr-5, -1
        do  n = nr, 3, -1

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
!$        n = nr-5
        n = 3

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

        fr=1  ! Random but must be postive. For the first check in the following loop
!$        n = nr-5
        n = 3
        nfr = 0
        if (alpha==1) nnod = 0

!________________alpha=1____________________________________
        if (alpha == 1) then

!$        do n = nr-5, 2, -1
        do n = 3, 2, -1
          J = a*(rofi(n) + b); r = rofi(n)
!C          J = J*q*q
!C          r = J/a - b

          g1p = g1in(n+5) - (3*dx/10.0_8)*(11*g1sn - 14*g1snm1 + 26*g1snm2 - 14*g1snm3 + 11*g1snm4)
          f1p = f1in(n+5) - (3*dx/10.0_8)*(11*f1sn - 14*f1snm1 + 26*f1snm2 - 14*f1snm3 + 11*f1snm4)
          g2p = g2in(n+5) - (3*dx/10.0_8)*(11*g2sn - 14*g2snm1 + 26*g2snm2 - 14*g2snm3 + 11*g2snm4)
          f2p = f2in(n+5) - (3*dx/10.0_8)*(11*f2sn - 14*f2snm1 + 26*f2snm2 - 14*f2snm3 + 11*f2snm4)

          eps1 = dabs(g1p)
          eps2 = dabs(f1p)
          eps3 = dabs(g2p)
          eps4 = dabs(f2p)

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

        if ((g1inmax <= dabs(g1c)).and.(n>1300)) g1inmax = dabs(g1c)
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

        dg1new = g1in(n-1) - g1in(n)
        dg1old = g1in(n) - g1in(n+1)
        fr = dg1new/dg1old

        if ((fr<0).and.(nfr==0)) then
           nfr=nfr+1
           ncsave = n
        endif

        if ((alpha==1).and.(g1in(n) /= 0)) then
           if ((g1in(n+1)/g1in(n))<=0) nnod = nnod + 1
        endif

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

      if (nfr==0) ncsave=1450

       nc = 1480 !1360 !1030 !ncsave
!      nrk = 6 ; nrmx = nc
       nrk = 6 ; nrmx = nr

      else
!________________alpha=2____________________________________-

!&        do n = nr-5, 2, -1
       do n = 3, 2, -1
          J = a*(rofi(n) + b); r = rofi(n)
!C          J = J*q*q
!C          r = J/a - b

          g1p = g1in(n+5) - (3*dx/10.0_8)*(11*g1sn - 14*g1snm1 + 26*g1snm2 - 14*g1snm3 + 11*g1snm4)
          f1p = f1in(n+5) - (3*dx/10.0_8)*(11*f1sn - 14*f1snm1 + 26*f1snm2 - 14*f1snm3 + 11*f1snm4)
          g2p = g2in(n+5) - (3*dx/10.0_8)*(11*g2sn - 14*g2snm1 + 26*g2snm2 - 14*g2snm3 + 11*g2snm4)
          f2p = f2in(n+5) - (3*dx/10.0_8)*(11*f2sn - 14*f2snm1 + 26*f2snm2 - 14*f2snm3 + 11*f2snm4)

          eps1 = dabs(g1p)
          eps2 = dabs(f1p)
          eps3 = dabs(g2p)
          eps4 = dabs(f2p)

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

          if ((g1inmax <= dabs(g1c)).and.(n>1300)) g1inmax = dabs(g1c)
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

          dg1new = g1in(n-1) - g1in(n)
          dg1old = g1in(n) - g1in(n+1)
          fr = dg1new/dg1old

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
      do ie=1,5
        psiin(:,1,alpha,ie) = g1in(:)
        psiin(:,2,alpha,ie) = f1in(:)
        psiin(:,3,alpha,ie) = g2in(:)
        psiin(:,4,alpha,ie) = f2in(:)
      enddo


!___________________End of Remnant of the initial rdeq__________

!___________________For normalization___________________________

      psiinn(:,1)=g1in(:)
      psiinn(:,2)=f1in(:)
      psiinn(:,3)=g2in(:)
      psiinn(:,4)=f2in(:)

!__________________End _For normalization________________________


!__________________Normalization_________________________________

      ncut = 1100

      call productn4(psioutn(:,:),psioutn(:,:),a,b,rofi,nr,normn)
      do j1=1,4
         psioutn(:,j1)=psioutn(:,j1)/dsqrt(normn)
      enddo

      normsave = dsqrt(normn)

!      call productn5(psiinn(:,:),psiinn(:,:),a,b,rofi,ncut,nr,normn)!1100
!      normn = 1
      do j1=1,4
         psiinn(:,j1)=psiinn(:,j1)/dsqrt(normn)
      enddo


!      ncut = 1250
!      call productn5(psiinn(:,:),psiinn(:,:),a,b,rofi,ncut,nr,normn)!1100
!      normn = normn + rofi(ncut)*abs(g1in(ncut))/2
!      do j1=1,4
!         psiinn(:,j1)=psiinn(:,j1)/dsqrt(normn)
!      enddo



      g1out(:)=psioutn(:,1)
      f1out(:)=psioutn(:,2)
      g2out(:)=psioutn(:,3)
      f2out(:)=psioutn(:,4)
      g1in(:)=psiinn(:,1)
      f1in(:)=psiinn(:,2)
      g2in(:)=psiinn(:,3)
      f2in(:)=psiinn(:,4)





!__________________End Normalization_______________________________




      if (alpha == 1) call printdata2('out_OutInt', 'rofi g1 12 12 2013', nr-1, rofi(2:), g1out(2:))

      do ie=1,5
        psiout(:,1,alpha,ie) = g1out(:)
        psiout(:,2,alpha,ie) = f1out(:)
        psiout(:,3,alpha,ie) = g2out(:)
        psiout(:,4,alpha,ie) = f2out(:)
        psiin(:,1,alpha,ie) = g1in(:)
        psiin(:,2,alpha,ie) = f1in(:)
        psiin(:,3,alpha,ie) = g2in(:)
        psiin(:,4,alpha,ie) = f2in(:)
      enddo

      end subroutine integration2


!___________F(P)=0______________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      subroutine newt(Pv,En,nF,xk, xk1, mu, up, u, upp, sqru, vm, bm, &
           ksop,z, v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,check,fvector)

      use newtv_common
      implicit none
!C____________________________________________________________________________
!C Given an initial guess for Pv(E,A,Aout,Ain) (see Ebert,J.Phys.:Cond.M.1(1989)
!C  9111-9116) finds the root of F = 0 (Ebert) by a globally convergent
!C Newton's method. Taken from 'Numerical Recipes in fortran 77, p 379.
!C----------------------------------------------------------------------------

!      integer nF,MAXITS
      integer nF,MAXITS
      logical check
      double precision Pv(nF),TOLF,TOLMIN,TOLX,STPMX

      parameter (MAXITS=300,TOLF=1.e-12_8,TOLMIN=1.e-6_8, TOLX=1.e-7_8, STPMX=100._8)
!       common /newtv/ fvec(NP),
!      . goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP),nn
!       save /newtv/
!CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
      integer i,its,j,indx(NP),nmrseq
      double precision d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),xold(NP), fvecold(NP) !,fmin
      integer nr,nsp,imu,l,lmx,nm
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision mu,rm,sqru,u,up,upp,xk,xk1
      double precision  bm(nr),vm(nr)

      double precision deltaf(NP), fvector(NP), deltaE, Eold, En
      procedure(real(8)) :: fmin
!       EXTERNAL fmin


      nn=nF
      f=fmin(Pv,En,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,  &
           v,rofi,nr,nmrseq,nsp,a,b,l, lmx, imu,scalede)

!      f = fmin(Pv)
!      call funcv(nF,Pv,fvec)
      test=0.0_8
      do i=1,nF
        if(abs(fvec(i)) > test)test=abs(fvec(i))
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

!        print*,'fvec= ',fvec
        call fdjac(nF,Pv,En,fjac, xk, xk1, mu, up, u, upp, sqru, vm, &
             bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede)
!        print*,'Jacobian = ',fjac
        do i=1,nF
          sum=0.0_8
          do j=1,nF
            sum=sum+fjac(j,i)*fvec(j)
          enddo
          g(i)=sum
        enddo
        do i=1,nF
          xold(i)=Pv(i)
          fvecold(i)=fvec(i)
        enddo
        Eold = Pv(1)
        fold=f
        do i=1,nF
          p(i)=-fvec(i)
        enddo

        call ludcmp(fjac,nF,NP,indx,d)
        call lubksb(fjac,nF,NP,indx,p)
        call lnsrch(nF,En,xold,fold,g,p,Pv,f,stpmax,check,fmin,&
            xk, xk1, mu, up, u,                                      &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,    &
            lmx, imu,scalede)
!        print*,'fvec = ',fvec


        fvector(:) = fvec(:)
      if ((l==1).and.(imu==2)) then
!      print*,'? '
!      print*,'Final values'
!      print*, 'newt - end','Energy=',Pv(1),' A = ',Pv(2), ' Aout = ', Pv(3),' Ain = ', Pv(4)
!      print*,'ln 2273 newt F = ',fvec
      endif

        deltaE = Pv(1)- Eold
        do i=1,nF
          deltaf(i) = fvec(i)-fvecold(i)
        enddo

        test=0.0_8
        do i=1,nF
          if(abs(fvec(i)) > test)test=abs(fvec(i))
        enddo
        if(test < TOLF)then
          write(1000,*)fvec(:),', test = ',test,'TOL=',TOLF
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
            write(1000,*)fvec(:),', test = ',test,'TOLM=',TOLMIN
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



      call pause('MAXITS exceeded in newt')
      end subroutine newt


      subroutine lnsrch(n,En,xold,fold,g,p,Pv,f,stpmax, &
        check,func, xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v,rofi,nr,  &
        nmrseq,nsp,ar,br,l, lmx, imu,scalede)

!C--------Taken from 'Numerical Recipes in fortran 77'---------------

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
      double precision v(nr,nsp),rofi(nr),z,ar,br,scalede,En
      double precision ksop(0:lmx,nsp,nsp,6)

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
      if(slope >= 0.0_8) call pause('roundoff problem in lnsrch')
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

        f=func(Pv,En,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v, &
             rofi,nr,nmrseq,nsp,ar,br,l,lmx, imu,scalede)

        if(alam < alamin)then
          do i=1,n
            Pv(i)=xold(i)
          enddo
          check=.true.
!           print*,' '
!           print*,'Final values'
!           print*, 'Energy =',Pv(1),' A = ',Pv(2),
!     .     ' Aout = ', Pv(3),' Ain = ', Pv(4)

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

      end subroutine lnsrch



      function fmin(Pv,En,xk, xk1, mu, up, u,upp, sqru, vm, bm, ksop,z, &
           v, rofi,nr, nmrseq,nsp,a,b,l,lmx, imu,scalede)

!C--------Taken from 'Numerical Recipes in fortran 77'------------

      use newtv_common
      implicit none
!       integer n
      double precision fmin,Pv(*)
!       ,fvec,
!      . goutnm, foutnm, ginnm, finnm
!       common /newtv/fvec(NP),
!      . goutnm(NP), foutnm(NP), ginnm(NP), finnm(NP),n
!       save /newtv/
      integer i
      double precision sum
      integer nr,nsp,imu,l,lmx,nm,nmrseq
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision mu,rm,sqru,u,up,upp,xk,xk1
      double precision  bm(nr),vm(nr),En


      call funcv(nn,Pv,En, xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,&
           z,v, rofi,nr,nmrseq,nsp,a,b,l, lmx, imu,scalede,fvec)


       sum=0.0_8
       do i=1,nn
         sum = sum + fvec(i)**2
       enddo
       fmin=0.5_8*sum

       return
      end function fmin

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      subroutine ludcmp(a,n,np,indx,d)

!C--------Taken from 'Numerical Recipes in fortran 77'---------------

      implicit none
      integer n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      parameter (NMAX=500,TINY=1.0e-20_8)
      integer i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)

!      print*,'ludcmp - n,np', n,'; ',np
      d=1.0_8
      do i=1,n
        aamax=0.0_8
        do j=1,n
          if (abs(a(i,j)) > aamax) aamax=abs(a(i,j))
        enddo
        if (aamax == 0.0_8) call pause('singular matrix in ludcmp')
        vv(i)=1./aamax
      enddo
      do  j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.0_8
        do  i=j,n
          sum=a(i,j)
          do  k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum >= aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j /= imax)then
          do  k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j) == 0.0_8)a(j,j)=TINY
        if(j /= n)then
          dum=1.0_8/a(j,j)
          do  i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      return
      END subroutine ludcmp

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      subroutine lubksb(a,n,np,indx,b)

!C--------Taken from 'Numerical Recipes in fortran 77'---------------

      implicit none
      integer n,np,indx(n)
      double precision a(np,np),b(n)
      integer i,ii,j,ll
      double precision sum

!      print*,'lubksb - n,np', n,'; ',np

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii /= 0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum /= 0.0_8) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
      return
      END subroutine lubksb


!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------

      subroutine funcv(n,Pv,En,xk, xk1, mu, up, u, upp, sqru, vm, bm, &
           ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,fvec)

!C--------Taken from 'Numerical Recipes in fortran 77'---------------

       implicit none
!       integer n
       double precision Pv(1:n),fvec(1:n)
       integer nr,nsp,imu,l,lmx,nm, nc,nmrseq
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6)

       integer n, nnod
       double precision mu,rm,sqru,u,up,upp,xk,xk1
       double precision g1in0, f1in0, g2in0, f2in0,normsave1,normsave2
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
         psiout2(nr,4,2,5),psiin2(nr,4,2,5),En


        call integration2(1, En, xk, xk1, mu, up, u,          &
            upp, sqru, vm,                                        &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,      &
            Pv(1),Pv(2), Pv(3),Pv(4),                             &
            g1out1,f1out1,g2out1,f2out1,g1sout1,f1sout1,g2sout1,  &
            f2sout1, g1in1,f1in1,g2in1,f2in1,g1sin1,f1sin1,g2sin1,&
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
        call integration2(2, En, xk, xk1, mu, up, u,          &
            upp, sqru, vm,                                        &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,      &
            Pv(1),Pv(2), Pv(3),Pv(4),                             &
            g1out2,f1out2,g2out2,f2out2,g1sout2,f1sout2,g2sout2,  &
            f2sout2, g1in2,f1in2,g2in2,f2in2,g1sin2,f1sin2,g2sin2,&
            f2sin2,psiout2, psiin2,nc, nnod,normsave2)


        nm=nmrseq

        fvec(1) =  g1out1(nm) - g1in1(nm)
        fvec(2) =  f1out1(nm) - f1in1(nm)
        fvec(3) =  g2out2(nm) - g2in2(nm)
        fvec(4) =  f2out2(nm) - f2in2(nm)


!        print*,' '
!        print*,'line 2235 f = ', fvec

      end subroutine funcv

!_______________________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------


      subroutine fdjac(n,Pv,En,df, xk, xk1, mu, up, u, upp, sqru, vm,&
           bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede)
!C---------------Provides Jacobian for minimization purposes---------
       use newtv_common
     implicit none
       integer, intent(in) :: n
       integer nc, nm, nnod
       integer nr, nsp, imu, l, lmx,nmrseq,i
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6)
       double precision mu, u, up, upp, xk, xk1, sqru
!        double precision fvec,goutnm, foutnm, ginnm, finnm
!        parameter (NP=4)
!        common /newtv/ fvec(NP),
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
         g1in1(nr), g2in1(nr),f1in1(nr),f2in1(nr),                       &
         g1out2(nr),g2out2(nr),f1out2(nr),f2out2(nr),                    &
         g1in2(nr), g2in2(nr),f1in2(nr),f2in2(nr),En
        double precision deltag1,deltaf1,deltag2,deltaf2,g1plus,f1plus,  &
             g2plus,f2plus,g1minus,f1minus,g2minus,f2minus,delta,        &
             normsave1,normsave2

         delta = 1e-4
         deltag1 = abs(Pv(1)*delta)
         deltaf1 = abs(Pv(2)*delta)
         deltag2 = abs(Pv(3)*delta)
         deltaf2 = abs(Pv(4)*delta)
!         print*,'Pv in fdjac =',Pv
         g1plus  = Pv(1) + deltag1; g1minus = Pv(1) - deltag1
         f1plus  = Pv(2) + deltaf1; f1minus = Pv(2) - deltaf1
         g2plus  = Pv(3) + deltag2; g2minus = Pv(3) - deltag2
         f2plus  = Pv(4) + deltaf2; f2minus = Pv(4) - deltaf2

         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            g1plus,Pv(2),Pv(3),Pv(4),                                        &
            g1out1p,f1out1p,g2out1p,f2out1p,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1p,f1in1p,g2in1p,f2in1p,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            g1plus,Pv(2),Pv(3),Pv(4),                                        &
            g1out2p,f1out2p,g2out2p,f2out2p,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2p,f1in2p,g2in2p,f2in2p,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin2,nc, nnod,normsave2)

         nm=nmrseq
!         print*,'nm = ',nm
         fvecp(1) =  g1out1p(nm) - g1in1p(nm)
         fvecp(2) =  f1out1p(nm) - f1in1p(nm)
         fvecp(3) =  g2out2p(nm) - g2in2p(nm)
         fvecp(4) =  f2out2p(nm) - f2in2p(nm)
!         print*,'fvecp = ',fvecp

         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            g1minus,Pv(2),Pv(3),Pv(4),                                       &
            g1out1m,f1out1m,g2out1m,f2out1m,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1m,f1in1m,g2in1m,f2in1m,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            g1minus,Pv(2),Pv(3),Pv(4),                                       &
            g1out2m,f1out2m,g2out2m,f2out2m,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2m,f1in2m,g2in2m,f2in2m,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin1,nc, nnod,normsave2)

         nm=nmrseq
         fvecm(1) =  g1out1m(nm) - g1in1m(nm)
         fvecm(2) =  f1out1m(nm) - f1in1m(nm)
         fvecm(3) =  g2out2m(nm) - g2in2m(nm)
         fvecm(4) =  f2out2m(nm) - f2in2m(nm)
!         print*,'fvecm = ',fvecm

         df(1,1) = (fvecp(1)-fvecm(1))/2/deltag1
         df(2,1) = (fvecp(2)-fvecm(2))/2/deltag1
         df(3,1) = (fvecp(3)-fvecm(3))/2/deltag1
         df(4,1) = (fvecp(4)-fvecm(4))/2/deltag1
!_______1_______________________________________________________________
         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),f1plus,Pv(3),Pv(4),                                        &
            g1out1p,f1out1p,g2out1p,f2out1p,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1p,f1in1p,g2in1p,f2in1p,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),f1plus,Pv(3),Pv(4),                                        &
            g1out2p,f1out2p,g2out2p,f2out2p,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2p,f1in2p,g2in2p,f2in2p,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin2,nc, nnod,normsave2)

         nm=nmrseq
         fvecp(1) =  g1out1p(nm) - g1in1p(nm)
         fvecp(2) =  f1out1p(nm) - f1in1p(nm)
         fvecp(3) =  g2out2p(nm) - g2in2p(nm)
         fvecp(4) =  f2out2p(nm) - f2in2p(nm)

         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),f1minus,Pv(3),Pv(4),                                       &
            g1out1m,f1out1m,g2out1m,f2out1m,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1m,f1in1m,g2in1m,f2in1m,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),f1minus,Pv(3),Pv(4),                                       &
            g1out2m,f1out2m,g2out2m,f2out2m,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2m,f1in2m,g2in2m,f2in2m,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin1,nc, nnod,normsave2)

         nm=nmrseq
         fvecm(1) =  g1out1m(nm) - g1in1m(nm)
         fvecm(2) =  f1out1m(nm) - f1in1m(nm)
         fvecm(3) =  g2out2m(nm) - g2in2m(nm)
         fvecm(4) =  f2out2m(nm) - f2in2m(nm)

         df(1,2) = (fvecp(1)-fvecm(1))/2/deltaf1
         df(2,2) = (fvecp(2)-fvecm(2))/2/deltaf1
         df(3,2) = (fvecp(3)-fvecm(3))/2/deltaf1
         df(4,2) = (fvecp(4)-fvecm(4))/2/deltaf1
!_______2_______________________________________________________________
         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),g2plus,Pv(4),                                        &
            g1out1p,f1out1p,g2out1p,f2out1p,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1p,f1in1p,g2in1p,f2in1p,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),g2plus,Pv(4),                                        &
            g1out2p,f1out2p,g2out2p,f2out2p,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2p,f1in2p,g2in2p,f2in2p,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin2,nc, nnod,normsave2)

         nm=nmrseq
         fvecp(1) =  g1out1p(nm) - g1in1p(nm)
         fvecp(2) =  f1out1p(nm) - f1in1p(nm)
         fvecp(3) =  g2out2p(nm) - g2in2p(nm)
         fvecp(4) =  f2out2p(nm) - f2in2p(nm)

         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),g2minus,Pv(4),                                       &
            g1out1m,f1out1m,g2out1m,f2out1m,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1m,f1in1m,g2in1m,f2in1m,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),g2minus,Pv(4),                                       &
            g1out2m,f1out2m,g2out2m,f2out2m,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2m,f1in2m,g2in2m,f2in2m,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin1,nc, nnod,normsave2)

         nm=nmrseq
         fvecm(1) =  g1out1m(nm) - g1in1m(nm)
         fvecm(2) =  f1out1m(nm) - f1in1m(nm)
         fvecm(3) =  g2out2m(nm) - g2in2m(nm)
         fvecm(4) =  f2out2m(nm) - f2in2m(nm)

         df(1,3) = (fvecp(1)-fvecm(1))/2/deltag2
         df(2,3) = (fvecp(2)-fvecm(2))/2/deltag2
         df(3,3) = (fvecp(3)-fvecm(3))/2/deltag2
         df(4,3) = (fvecp(4)-fvecm(4))/2/deltag2
!_______3_______________________________________________________________

         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),Pv(3),f2plus,                                        &
            g1out1p,f1out1p,g2out1p,f2out1p,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1p,f1in1p,g2in1p,f2in1p,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),Pv(3),f2plus,                                        &
            g1out2p,f1out2p,g2out2p,f2out2p,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2p,f1in2p,g2in2p,f2in2p,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin2,nc, nnod,normsave2)

         nm=nmrseq
         fvecp(1) =  g1out1p(nm) - g1in1p(nm)
         fvecp(2) =  f1out1p(nm) - f1in1p(nm)
         fvecp(3) =  g2out2p(nm) - g2in2p(nm)
         fvecp(4) =  f2out2p(nm) - f2in2p(nm)

         call integration2(1, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),Pv(3),f2minus,                                       &
            g1out1m,f1out1m,g2out1m,f2out1m,g1sout1,f1sout1,g2sout1,         &
            f2sout1, g1in1m,f1in1m,g2in1m,f2in1m,g1sin1,f1sin1,g2sin1,       &
            f2sin1,psiout1, psiin1,nc, nnod,normsave1)
         call integration2(2, En, xk, xk1, mu, up, u,                    &
            upp, sqru, vm,                                                   &
            bm, ksop,z,v,rofi,nr,nsp,a,b,l,lmx, imu,scalede,                 &
            Pv(1),Pv(2),Pv(3),f2minus,                                       &
            g1out2m,f1out2m,g2out2m,f2out2m,g1sout2,f1sout2,g2sout2,         &
            f2sout2, g1in2m,f1in2m,g2in2m,f2in2m,g1sin2,f1sin2,g2sin2,       &
            f2sin2,psiout2, psiin1,nc, nnod,normsave2)

         nm=nmrseq
         fvecm(1) =  g1out1m(nm) - g1in1m(nm)
         fvecm(2) =  f1out1m(nm) - f1in1m(nm)
         fvecm(3) =  g2out2m(nm) - g2in2m(nm)
         fvecm(4) =  f2out2m(nm) - f2in2m(nm)

         df(1,4) = (fvecp(1)-fvecm(1))/2/deltaf2
         df(2,4) = (fvecp(2)-fvecm(2))/2/deltaf2
         df(3,4) = (fvecp(3)-fvecm(3))/2/deltaf2
         df(4,4) = (fvecp(4)-fvecm(4))/2/deltaf2
!_______4_______________________________________________________________


      end subroutine fdjac

!_______end of F(Pv)=0________________________________________________________________________________________
!-------------------------------------------------------------------------------------------------------


!______Another equations system solver____________________________


      subroutine usrfun1(Pv,En,n,NP,fvec,fjac, xk, xk1, mu, up, u,  &
            upp, sqru, vm,                                       &
            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede)
        implicit none
       integer n, NP
       integer nr, nsp, imu, l, lmx,nmrseq
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6),Pv(1:n)
       double precision mu, u, up, upp, xk, xk1, sqru
!      parameter (NP=4)
       double precision fjac(n,n),fvec(n),bm(nr),vm(nr)
       double precision goutnm(4), foutnm(4), ginnm(4), finnm(4),En


       call funcv(n,Pv,En,xk, xk1, mu, up, u,                  &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,    &
            lmx, imu,scalede,fvec)

!       print*,'Pv in usrfun1 = ',Pv

       call fdjac(n,Pv,En,fjac, xk, xk1, mu, up, u,            &
            upp, sqru, vm,                                           &
            bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede)
!       print*,'Jac in usrfun1 = ',fjac

      end subroutine usrfun1



      subroutine mnewt(ntrial,x,En,n,tolx,tolf, xk, xk1, mu, up, u,&
             upp, sqru, vm,                                             &
             bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,fvec)
      implicit none
       integer n, NP, ntrial
       double precision tolf, tolx, x(n)
       parameter (NP=4)
       double precision d,errf, errx,fjac(NP,NP),fvec(NP),p(NP)
       integer i, k, indx(NP)
       integer nr, nsp, imu, l, lmx, nmrseq
       double precision v(nr,nsp),rofi(nr),z,a,b,scalede
       double precision ksop(0:lmx,nsp,nsp,6)
       double precision mu, u, up, upp, xk, xk1, sqru
!      parameter (NP=4)
       double precision Pv(n),bm(nr),vm(nr),En

!       print*,'In mnewt'

       do k=1,ntrial
!          print*,'k=',k
          call usrfun1(x,En,n,NP,fvec,fjac, xk, xk1, mu, up, u, &
             upp, sqru, vm,bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede)
          errf = 0
          do i=1,n
             errf = errf+abs(fvec(i))
          enddo

!          write(3,*)' k = ',k
!          write(3,*)' E = ',x(1)
!          write(3,*)' F = ',fvec(:)

          if (errf <= tolf) return
          do i=1,n
             p(i)=-fvec(i)
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

      end subroutine mnewt




!______Another equations system solver done________________________



!______Third minimization procedure________________________________

      subroutine broydn(x,En,n,check, xk, xk1, mu, up, u, upp, sqru, &
           vm,bm,ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede,fsave)

!C--------Taken from 'Numerical Recipes in fortran 77'---------------

        use newtv_common
      implicit none
!     x = Pv

      integer n,nmrseq
      integer, parameter :: MAXITS=200
      double precision x(n)
      real(8), parameter :: EPS=1.e-7_8,TOLF=1.e-6_8,TOLMIN=1.e-6_8, STPMX=100.0_8
      real(8), parameter :: TOLX=EPS
      logical check
!       common /newtv/ fvec(NP),nn
!       save /newtv/
!CU    uses fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
      integer i,its,j,k
      double precision den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP), &
             g(NP),p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP) !,fmin
      logical restrt,sing,skip
!       external fmin

      integer nr, nsp, imu, l, lmx
      double precision v(nr,nsp),rofi(nr),z,a,b,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision mu, u, up, upp, xk, xk1, sqru
      double precision Pv(n),bm(nr),vm(nr),fsave(n),En

      procedure(real(8)) :: fmin

      nn=n
      f=fmin(x,En,xk, xk1, mu, up, u, upp, sqru, vm, bm, ksop,z,v, &
           rofi,nr,nmrseq,nsp,a,b,l, lmx, imu,scalede)
      test=0.0_8
      do i=1,n
        if(abs(fvec(i)) > test)test=abs(fvec(i))
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
          call fdjac(n,x,En,r,xk, xk1, mu, up, u, upp, sqru, vm, bm, &
               ksop,z,v, rofi,nr,nmrseq,nsp,a,b,l,lmx, imu,scalede)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) call pause('singular Jacobian in broydn')
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
            w(i)=fvec(i)-fvcold(i)-sum
            if(abs(w(i)) >= EPS*(abs(fvec(i))+abs(fvcold(i))))then
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
              if(r(i,i) == 0.0_8) call pause('r singular in broydn')
              d(i)=r(i,i)
            enddo
          endif
        endif
        do i=1,n
          sum=0.0_8
          do  j=1,n
            sum=sum+qt(i,j)*fvec(j)
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
          fvcold(i)=fvec(i)
        enddo
        fold=f
!        do  i=1,n
!          sum=0.
!          do j=1,n
!            sum=sum+qt(i,j)*fvec(j)
!          enddo
!          p(i)=-sum
!        enddo
        call rsolv(r,n,NP,d,p)
        call lnsrch(n,En,xold,fold,g,p,x,f,stpmax,check,fmin,     &
            xk, xk1, mu, up, u,                                         &
            upp, sqru, vm, bm, ksop,z,v,rofi,nr,nmrseq,nsp,a,b,l,       &
            lmx, imu,scalede)

        fsave(:) = fvec(:)

        test=0.0_8
        do  i=1,n
          if(abs(fvec(i)) > test)test=abs(fvec(i))
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
      call pause('MAXITS exceeded in broydn')
      END subroutine broydn

      subroutine qrdcmp(a,n,np,c,d,sing)
      implicit none
      integer n,np
      double precision a(np,np),c(n),d(n)
      logical sing
      integer i,j,k
      double precision scale,sigma,sum,tau

      sing=.false.
      do  k=1,n-1
        scale=0.0_8
        do  i=k,n
          scale=max(scale,abs(a(i,k)))
        enddo
        if(scale == 0.0_8)then
          sing=.true.
          c(k)=0.0_8
          d(k)=0.0_8
        else
          do  i=k,n
            a(i,k)=a(i,k)/scale
          enddo
          sum=0.0_8
          do  i=k,n
            sum=sum+a(i,k)**2
          enddo
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do  j=k+1,n
            sum=0.0_8
            do  i=k,n
              sum=sum+a(i,k)*a(i,j)
            enddo
            tau=sum/c(k)
            do i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
            enddo
          enddo
        endif
      enddo
      d(n)=a(n,n)
      if(d(n) == 0.)sing=.true.
      return
      END subroutine qrdcmp

      subroutine qrupdt(r,qt,n,np,u,v)
      implicit none
      integer n,np
      double precision r(np,np),qt(np,np),u(np),v(np)
!CU    USES rotate
      integer i,j,k
      do  k=n,1,-1
        if(u(k) /= 0.0_8)goto 1
      enddo
      k=1
1     do i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i) == 0.0_8)then
          u(i)=abs(u(i+1))
        else if(abs(u(i)) > abs(u(i+1)))then
          u(i)=abs(u(i))*sqrt(1.0_8+(u(i+1)/u(i))**2)
        else
          u(i)=abs(u(i+1))*sqrt(1.0_8+(u(i)/u(i+1))**2)
        endif
      enddo
      do j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
      enddo
      do  i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
      enddo
      return
      END subroutine qrupdt

      subroutine rsolv(a,n,np,d,b)
      implicit none
      integer n,np
      double precision a(np,np),b(n),d(n)
      integer i,j
      double precision sum
      b(n)=b(n)/d(n)
      do  i=n-1,1,-1
        sum=0.0_8
        do  j=i+1,n
          sum=sum+a(i,j)*b(j)
        enddo
        b(i)=(b(i)-sum)/d(i)
      enddo
      return
      END subroutine rsolv

      subroutine rotate(r,qt,n,np,i,a,b)
      implicit none
      integer n,np,i
      double precision a,b,r(np,np),qt(np,np)
      integer j
      double precision c,fact,s,w,y

      if(a == 0.0_8)then
        c=0.0_8
        s=sign(1.0_8,b)
      else if(abs(a) > abs(b))then
        fact=b/a
        c=sign(1.0_8/sqrt(1.0_8+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1.0_8/sqrt(1.0_8+fact**2),b)
        c=fact*s
      endif
      do  j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
      enddo
      do  j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
      enddo
      return
      END subroutine rotate




!______End of the third minimization procedure_____________________

!C      subroutine printdata2(fname, desc, nini, nend, d1, d2)
      subroutine printdata2(fname, desc, n, d1, d2)
        implicit none
        character(len=*), intent(in) :: fname, desc
        integer, intent(in) :: n
        real(8), intent(in) :: d1(1:n), d2(1:n)

        integer :: i, ifi
        procedure(integer) :: fopna

        ifi = fopna(fname,-1,0)

        write(ifi, '("# ",a)') desc
        do i = 1,n
           write(ifi, *) d1(i), d2(i)
        end do

        call fclose(ifi)
      end subroutine printdata2

      subroutine pause(m)
          implicit none
          character(len=*) :: m
          character :: inp
          print *, m
          read(*,*) inp
      end subroutine pause
!     end module m_rdeq


      subroutine clebsh(l,imu,a,a1)
      implicit none
      double precision a(2,2),a1(2,2),cl(2,2),mu,u
!C     double precision cl1(2,2)
      integer ms1,ms2,l,imu


      mu= dble(imu - l) - 1.5d0
      u = mu/(dble(l)+0.5d0)

      cl(1,1) = dsqrt(0.5d0*(1d0+u))
      cl(2,2) = cl(1,1)
      cl(2,1) = -dsqrt(0.5d0*(1d0-u))
      cl(1,2) = -cl(2,1)

!C         call dinv22(cl,cl1)

      do  ms1 = 1, 2
        do  ms2 = 1, 2
          a1(ms1,ms2) = cl(1,ms1)*a(1,1)*cl(1,ms2)   &
                      + cl(1,ms1)*a(1,2)*cl(2,ms2)   &
                      + cl(2,ms1)*a(2,1)*cl(1,ms2)   &
                      + cl(2,ms1)*a(2,2)*cl(2,ms2)
        enddo
      enddo

      end




