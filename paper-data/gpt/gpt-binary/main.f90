program mdgpt

  use iso_fortran_env, only: ip => int32, dp => real64
  !use mpi
  implicit none

  type agpt
     character(len=128) :: dnfile
     integer(ip) :: mesh, nvol, nbas
     logical :: do, derivative
     integer(ip) :: vptr
     integer(ip), allocatable :: species(:) 
     real(dp) :: ra, omega, plat(3,3), step, rws, nval0, zval, &
          nvala, dnvala
     real(dp), allocatable :: dns(:,:,:), volumes(:), d2(:,:), &
          basis(:,:), nvals(:)
  end type agpt
  
  type cell
     character(len=128), allocatable :: filename(:,:)
     integer(ip) :: nbas, nspec
     integer(ip), allocatable :: species(:)
     real(dp) :: alat, plat(3,3), volume, omega, rws, rmax, rmaxrws, &
          dims(3), cent(3), iplat(3,3)
     real(dp), allocatable :: basis(:,:), xbasis(:,:), mass(:)
  end type cell

  type io
     integer(ip) :: nza, mesh, nvol, vptr
     real(dp) :: mass, p1, al, va, vb, p2, p4, r0, rp, pn, r0rws, &
          step, vol, rws, evol, vstep, evold, evoldd
     real(dp), allocatable :: numpot(:,:,:), volumes(:), evolry(:), &
          v2(:,:)
  end type io

  type md
     integer(ip) :: ensemble, nsteps, wsteps, ncol(2)
     integer(ip), allocatable :: colours(:)
     logical :: do, cmv, ne
     real(dp) :: itemp, temp, vsum(3), ke, pe, ham, tstep, e0, epair, &
          etot, ptot, pst(3,3), time, fsum, s, sdot, colourf(-1:1,0:3), &
          taut, friction, engint, conq, msd, pext, enthalpy, eta, piston, &
          taup
     real(dp), allocatable :: vels(:,:), force(:,:), pos0(:,:)
  end type md
  
  type nhbr
     integer(ip) :: ms(3)
     integer(ip), allocatable :: tot(:), ptr(:), head(:), list(:)
     real(dp) :: cent(3)
  end type nhbr
  
  type relax
     integer(ip) :: term, nmin
     logical :: do, elastic, screen, interpolate, linkl
     logical, allocatable :: flags(:,:)
     real(dp) :: fbreak, e0, epair, etot, tstep, tmax, f(3), mass, &
          acoef0, acoef, sq(3), pst(3,3), ptot
     real(dp), allocatable :: force(:,:), vels(:,:)
  end type relax

  ! derived types
  type(agpt) :: s_agpt
  type(cell) :: s_cell
  type(io), allocatable :: s_io(:,:)
  type(md) :: s_md
  type(nhbr) :: s_nhbr, s_nhbragpt
  type(relax) :: s_relax
  ! parameters
  real(dp), parameter :: third = 1/3.0_dp, pi = 4*atan(1.0_dp), &
       amass = 1.09716e-3_dp, fsec = 0.04837768653_dp, &
       rydbohrtogpa = 14710.516405757675_dp, &
       boltzmann = 6.333616e-06_dp
  ! allocatables
  integer(ip), allocatable :: basptr(:)
  real(dp), allocatable :: sites(:,:)
  ! locals
  integer(ip) :: i, j, k, l, icnt, site_count, nmax, ibsp, jbsp, ptr, &
       rcut(3), mu, nu, unit, a, b, c, d, idxs(2), jdxs(2), rptr, imd, &
       xyzu, clr, vptr, vptri, rlxu, foru, smax, sext(3), icell, ierr, &
       nsum, ksum
  logical :: rst
  real(dp) :: vdff(2), vol, vol1, r, rrws, rsqr, flr, fl2, v2a, v2b, &
       fscr, dq1, dp1, arg, arg1, arg2, arg12, f, arg13, fscr1, fscr2, &
       pa1, pa2, pb1, pb2, vv, vf, ff, ra(3), crad, ctan, cprime, &
       epp, fmax, cv, cf, qr, sf, cij(6,6), bmod, v2d, v2dd, dv2d, fv, &
       fkn(3,3), st, cijkl(3,3,3,3), rws, fij(3), fijk(3), sijab(3,3), &
       sijkab(3,3), omegai, nisi, chni, f0, sitex(3), wallt(2), inbas, &
       ivol, kbt, dt2, pvd, mass, etadt, sts(3,3)

  !call mpi_init(ierr)
  call cpu_time(wallt(1))
  
  call init_random_seed()
  call make_cell(s_agpt,s_cell,s_md,s_relax)

  if (s_agpt%do .and. s_cell%nspec > 1) then
     print '(a)', 'aGPT does not yet support alloys'
     call exit(1)
  end if

  s_cell%volume = s_cell%alat**3.0_dp*det(s_cell%plat)
  s_cell%omega = s_cell%volume/s_cell%nbas
  s_cell%rws = (0.75_dp*s_cell%omega/pi)**third
  s_cell%rmax = s_cell%rmaxrws*s_cell%rws/s_cell%alat
  s_cell%cent = (/0.5_dp,0.5_dp,0.5_dp/)
  s_cell%iplat = inv3(s_cell%plat)
  inbas = 1.0_dp/s_cell%nbas
  ivol = 1.0_dp/s_cell%volume
  
  allocate(s_cell%xbasis(s_cell%nbas,3))

  do i = 1, s_cell%nbas
     s_cell%xbasis(i,:) = matmul(s_cell%basis(i,:),s_cell%iplat(:,:))
  end do
  
  if (s_agpt%do) then
     print '(a)', 'Molecular Statics and Dynamics (aGPT)'
  else
     print '(a)', 'Molecular Statics and Dynamics (GPT)'
  end if

  allocate(s_io(s_cell%nspec,s_cell%nspec))

  do i = 1, s_cell%nspec
     do j = i, s_cell%nspec
        call read_potential(s_cell%filename(i,j),s_io(i,j))
        ! find vptr
        vdff(1) = 200.0_dp
        do k = 1, s_io(i,j)%nvol
           vdff(2) = abs(s_io(i,j)%volumes(k)-s_cell%omega)
           if (vdff(2) > vdff(1)) cycle
           s_io(i,j)%vptr = k
           vdff(1) = vdff(2)
        end do
        ! find vstep
        vol = s_io(i,j)%volumes(s_io(i,j)%vptr)
        vol1 = s_io(i,j)%volumes(s_io(i,j)%vptr-1)
        s_io(i,j)%vstep = (vol1/vol)**third-1.0_dp
        ! apply screening
        if (s_relax%screen) then
           do k = 1, s_io(i,j)%nvol
              rws = (0.75_dp*s_io(i,j)%volumes(k)/pi)**(third)
              s_io(i,j)%r0 = s_io(i,j)%r0rws*rws
              s_io(i,j)%rp = 1.8_dp*rws
              do l = 1, s_io(i,j)%mesh
                 r = s_io(i,j)%numpot(k,l,1)
                 rrws = r/rws
                 rsqr = r*r
                 flr = fl(r,s_io(i,j))
                 fl2 = flr*flr
                 v2a = s_io(i,j)%va*fl2*fl2
                 v2b = s_io(i,j)%vb*fl2
                 if (rrws < s_io(i,j)%r0rws) then
                    fscr = 1.0_dp
                    dq1 = 0.0_dp
                    dp1 = 0.0_dp
                 else
                    arg = rrws/s_io(i,j)%r0rws
                    arg1 = arg-1.0_dp
                    arg2 = arg*arg
                    arg12 = arg1*arg1
                    f = hgauss(arg,s_io(i,j))
                    arg13 = arg1*arg12
                    dp1 = 2.0_dp*(s_io(i,j)%al**2)*arg*&
                         arg13/(1.0_dp+s_io(i,j)%al*arg12)
                    dq1 = 2.0_dp*(s_io(i,j)%al**2)*arg2*arg12* &
                         (2.0_dp*s_io(i,j)%al*arg12-3.0_dp)/(1.0_dp+s_io(i,j)%al*arg12)
                    fscr = f*f
                 end if
                 fscr1 = -2.0_dp*dp1*fscr
                 fscr2 = 2.0_dp*(dp1*dp1+dq1)*fscr
                 pa1 = s_io(i,j)%p4 + 4.0_dp*dp1
                 pb1 = s_io(i,j)%p2 + 2.0_dp*dp1
                 pa2 = s_io(i,j)%p4*(s_io(i,j)%p4+1.0_dp) &
                      +4.0_dp*(2.0_dp*s_io(i,j)%p4*dp1+3.0_dp*dp1*dp1+dq1)
                 pb2 = s_io(i,j)%p2*(s_io(i,j)%p2+1.0_dp)+ &
                      2.0_dp*(s_io(i,j)%p4*dp1+dp1*dp1+dq1)
                 s_io(i,j)%numpot(k,l,4) = fscr*s_io(i,j)%numpot(k,l,4) &
                      + 2.0_dp*fscr1*s_io(i,j)%numpot(k,l,3) &
                      + fscr2*s_io(i,j)%numpot(k,l,2)/rsqr &
                      + (pa2*v2a - pb2*v2b)/rsqr
                 s_io(i,j)%numpot(k,l,3) = fscr*s_io(i,j)%numpot(k,l,3) &
                      + fscr1*s_io(i,j)%numpot(k,l,2)/rsqr &
                      + (-pa1*v2a+pb1*v2b)/rsqr
                 s_io(i,j)%numpot(k,l,2) = fscr*s_io(i,j)%numpot(k,l,2)+v2a-v2b
              end do
           end do
        end if
        ! get step
        s_io(i,j)%step = (s_io(i,j)%numpot(1,s_io(i,j)%mesh,1) &
             -s_io(i,j)%numpot(1,1,1))/(s_io(i,j)%mesh-1)
        if (s_agpt%do) then
           call read_density(s_agpt)
           ! find vptr
           vdff(1) = 200.0_dp
           do k = 1, s_agpt%nvol
              vdff(2) = abs(s_agpt%volumes(k)-s_cell%omega)
              if (vdff(2) > vdff(1)) cycle
              s_agpt%vptr = k
              vdff(1) = vdff(2)
           end do
           allocate(s_agpt%d2(s_agpt%mesh,5))
           if (s_relax%interpolate) then
              ! interpolate d2
              s_agpt%omega = s_cell%omega
              s_agpt%rws = s_cell%rws
           else
              ! interpolate d2
              s_agpt%omega = s_agpt%volumes(s_agpt%vptr)
              s_agpt%rws = (0.75_dp*s_agpt%omega/pi)**(third)
           end if
           ! calculate d2
           do k = 1, s_agpt%mesh
              s_agpt%d2(k,1) = s_agpt%dns(s_agpt%vptr,k,1)
              s_agpt%d2(k,2) = lagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,2), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
              s_agpt%d2(k,3) = lagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,4), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
              s_agpt%d2(k,4) = lagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,5), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
              s_agpt%d2(k,5) = dlagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,2), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
           end do
           s_agpt%nval0 = s_agpt%zval/s_cell%omega
           s_agpt%step = (s_agpt%d2(s_agpt%mesh,1)-s_agpt%d2(1,1))/(s_agpt%mesh-1)
        end if
        ! calculate volume derivatives of v2
        do k = 1, s_io(i,j)%nvol
           vol = s_io(i,j)%volumes(k)
           if (k < 3) then
              vptr = k+2
           else if (k > s_io(i,j)%nvol-1) then
              vptr = k-1
           else
              vptr = k
           end if
           do l = 1, s_io(i,j)%mesh
              s_io(i,j)%numpot(k,l,5) = dlagrange(vol, &
                   s_io(i,j)%numpot(vptr-2:vptr+1,l,2), &
                   s_io(i,j)%volumes(vptr-2:vptr+1),4)
           end do
        end do
        ! interpolate v2
        if (.not. s_agpt%do) then
           if (s_relax%interpolate) then
              s_io(i,j)%vol = s_cell%omega
              s_io(i,j)%rws = s_cell%rws
           else
              s_io(i,j)%vol = s_io(i,j)%volumes(s_io(i,j)%vptr)
              s_io(i,j)%rws = (0.75_dp*s_io(i,j)%vol/pi)**(third)
           end if
           ! allocate v2 and derivatives
           allocate(s_io(i,j)%v2(s_io(i,j)%mesh,5))
           s_io(i,j)%evol = lagrange(s_io(i,j)%vol, &
                s_io(i,j)%evolry(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1), &
                s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
           s_io(i,j)%evold = dlagrange(s_io(i,j)%vol, &
                s_io(i,j)%evolry(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1), &
                s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
           do k = 1, s_io(i,j)%mesh
              s_io(i,j)%v2(k,1) = s_io(i,j)%numpot(s_io(i,j)%vptr,k,1)
              do l = 2, 5
                 s_io(i,j)%v2(k,l) = lagrange(s_io(i,j)%vol, &
                      s_io(i,j)%numpot(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1,k,l), &
                      s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
              end do
           end do
        end if
        ! make symmetric
        s_io(j,i) = s_io(i,j)
     end do
  end do

  ! do relaxation
  allocate(s_relax%vels(s_cell%nbas,3), &
       s_relax%force(s_cell%nbas,3))

  s_relax%tstep = 1/(s_cell%alat*fsec)
  s_relax%tmax = 10*s_relax%tstep
  s_relax%f = (/1.1_dp, 0.5_dp, 0.99_dp/)
  s_relax%nmin = 5
  s_relax%vels(:,:) = 0.0_dp
  s_relax%mass = 4/amass
  s_relax%acoef0 = 0.1_dp
  s_relax%acoef = s_relax%acoef0
  s_relax%term = 10000
  s_relax%sq = 2.0_dp*pi*(/1.0_dp,1.0_dp,0.0_dp/)

  ! create site list
  s_nhbr%cent(:) = 0.0_dp
  do i = 1, 3
     rcut(i) = int(s_cell%rmax/norm(s_cell%plat(i,:)),ip)+2
     sext(i) = int(2*rcut(i)+1,ip)
  end do
  s_nhbr%cent(:) = 0.5_dp*sext(:)-s_cell%cent(:)

  ! Linked List Block Size
  call llblocks(s_nhbr%ms,s_cell%plat,s_cell%rmax,sext)
  
  smax = s_cell%nbas*product(sext)
  nmax = int(650*s_cell%nbas,ip)
  allocate(sites(smax,3),basptr(smax), &
       s_nhbr%list(smax),s_nhbr%head(product(s_nhbr%ms)))
  
  if (s_agpt%do) then
     call agpt_on_site(s_agpt,s_cell,s_nhbr)
     print 568, s_agpt%omega, s_cell%omega, s_cell%alat, s_cell%rmax, &
          s_agpt%ra, s_cell%nbas, s_agpt%nval0, s_agpt%nvala
568  format(' Ps-atom volume   : ',f10.6,', Cell atomic volume : ',f10.6 / &
            ' Lattice constant : ',f10.6,', Neighbour cut off  : ',f10.6 / &
            ' Averaging radius : ',f10.6,', Number of atoms    : ',i10   / &
            ' Uniform density  : ',f10.6,', On-site density    : ',f10.6)
     allocate(s_agpt%nvals(s_cell%nbas))
  else
     print 567, s_io(1,1)%vol, s_cell%omega, s_cell%alat, s_cell%rmax, &
          s_io(1,1)%evol, s_cell%nbas, s_io(1,1)%evold, smax
567  format(' PP atomic volume : ',f10.6,', Cell atomic volume : ',f10.6 / &
            ' Lattice constant : ',f10.6,', Neighbour cut off  : ',f10.6 / &
            ' E_vol            : ',f10.6,', Number of atoms    : ',i10   / &
            ' E_vol (1st der.) : ',f10.6,', Ghost cell size    : ',i10)
  end if

  allocate(s_nhbr%ptr(nmax),s_nhbr%tot(s_cell%nbas))

  if (s_relax%linkl) print '(a,3(i6,x))', ' Block dimensions : ', s_nhbr%ms

  print '(a)', ' Relaxation (1 Step if no relax)'
  print '(a)', ' Step     F_RMS        E_vol      E_pair     E_tot      P_tot      S(q)'

  open(newunit=xyzu,file='relax.xyz',action='write')
  open(newunit=rlxu,file='relax.out',action='write')
  write(rlxu,'(a)') '# rms force, total energy and pressure'
  
976 format('Lattice="',9(x,f12.6),'" Properties=species:S:1:pos:R:3 Time=',f10.6)
  
  rlx: do icnt = 1, s_relax%term

     sites(:,:) = 0.0_dp; site_count = 1
     s_nhbr%head(:) = 0
     do i = 1, s_cell%nbas
        do j = -rcut(1), rcut(1)
           do k = -rcut(2), rcut(2)
              do l = -rcut(3), rcut(3)
                 sites(site_count,:) = s_cell%basis(i,:) &
                      + j*s_cell%plat(1,:) &
                      + k*s_cell%plat(2,:) &
                      + l*s_cell%plat(3,:)
                 basptr(site_count) = i
                 sitex(:) = s_cell%xbasis(i,:) + real((/j,k,l/),dp)
                 ! dislocated cell fudge
                 icell = find_cell(s_nhbr,sext,sitex)
                 s_nhbr%list(site_count) = s_nhbr%head(icell)
                 s_nhbr%head(icell) = site_count
                 site_count = site_count + 1
              end do
           end do
        end do
     end do
     site_count = site_count-1

     ! make neighbours
     if (s_relax%linkl) then
        call linkl_nhbr(s_cell,s_nhbr,s_nhbr%ms,sext,sites)
     else
        call nhbr_maker(s_cell%basis,s_nhbr,1,s_cell%nbas,site_count, &
             sites,s_cell%rmax)
     end if
        
     if (s_agpt%do) then
        call density(s_agpt,s_cell,s_nhbr)
        s_relax%e0 = 0.0_dp
        s_relax%epair = 0.0_dp
        s_relax%pst = 0.0_dp
     else
        s_relax%e0 = s_cell%nbas*s_io(1,1)%evol
        s_relax%epair = 0.0_dp; s_relax%pst = 0.0_dp
        forall (i = 1:3) s_relax%pst(i,i) = s_io(1,1)%evold
     end if

     vv = 0.0_dp; ff = 0.0_dp; vf = 0.0_dp; sf = 0.0_dp
     nsum = 0; ksum = 0
     do i = 1, s_cell%nbas
        s_relax%force(i,:) = 0.0_dp
        ibsp = s_cell%species(i)
        if (s_agpt%do) then                                                                    
           s_relax%e0 = s_relax%e0 + evol(s_io(1,1), &
                s_agpt%nvals(i),s_agpt%zval)
           nisi = s_agpt%nvals(i)
           omegai = s_agpt%zval/nisi
           chni = -omegai/nisi
           vptri = vpointer(omegai,s_io(ibsp,ibsp))
           f0 = chni*dlagrange(omegai, &                                                       
                s_io(ibsp,ibsp)%evolry(vptri-2:vptri+1), &                                     
                s_io(ibsp,ibsp)%volumes(vptri-2:vptri+1),4)                                    
           forall (j = 1:3) s_relax%pst(j,j) = s_relax%pst(j,j) &                              
                + f0*s_agpt%dnvala*inbas                                                 
        end if
        do j = 1, s_nhbr%tot(i)
           ptr = s_nhbr%ptr(nsum+j)
           jbsp = s_cell%species(basptr(ptr))
           ra = s_cell%alat*(sites(ptr,:)-s_cell%basis(i,:))                                
           r = norm(ra)
           if (s_agpt%do) then                                                              
              call agpt_force(s_agpt,s_cell,s_io(ibsp,jbsp),s_nhbr, &                       
                   i,ra,epp,fij,fijk,sijab,sijkab,ptr,sites,basptr, &
                   ksum)
              s_relax%force(i,:) = s_relax%force(i,:)+fij(:)+fijk(:)
              s_relax%epair = s_relax%epair + 0.5_dp*epp
              s_relax%pst(:,:) = s_relax%pst(:,:) + sijab(:,:) &                            
                   + sijkab(:,:)                                                            
           else
              call force_constants(s_io(ibsp,jbsp),r,crad,ctan,epp,v2d)
              s_relax%force(i,:) = s_relax%force(i,:) + ctan*ra(:)                          
              s_relax%epair = s_relax%epair + 0.5_dp*epp                                    
              forall (i = 1:3) s_relax%pst(i,i) = s_relax%pst(i,i) &                     
                   + 0.5_dp*v2d*inbas
              do mu = 1, 3                                                                  
                 do nu = 1, 3                                                               
                    s_relax%pst(mu,nu) = s_relax%pst(mu,nu) &                               
                         + 0.5_dp*ctan*ra(mu)*ra(nu)*ivol
                 end do
              end do
           end if
        end do
        do a = 1, 3                                                                            
           if (s_relax%flags(i,a)) then                                                        
              vv = vv + s_relax%vels(i,a)**2                                                   
              vf = vf + s_relax%force(i,a)*s_relax%vels(i,a)                                   
              ff = ff + s_relax%force(i,a)**2                                                  
           end if
        end do
        qr = dot_product(s_relax%sq,s_cell%basis(i,:))                                         
        sf = sf + cos(qr)*inbas
        nsum = nsum + s_nhbr%tot(i)
        ksum = ksum + s_nhbr%tot(i)
     end do
     
     s_relax%etot = s_relax%e0 + s_relax%epair
     s_relax%ptot = -third*trace(s_relax%pst,3)
     fmax = sqrt(third*ff*inbas)
     
     print '(x,i6,x,es11.4,6(x,f10.6))', icnt, fmax, &
          s_relax%e0*inbas, &
          s_relax%epair*inbas, &
          s_relax%etot*inbas, &
          s_relax%ptot*rydbohrtogpa, &
          sf

     ! write relax out
     write(rlxu,'(i6,x,es27.16,2(x,f24.16))') icnt, fmax, &
          s_relax%etot, s_relax%ptot*rydbohrtogpa
     
     ! change at some point - write xyz
     write(xyzu,'(i8)') s_cell%nbas
     write(xyzu,976) s_cell%alat*transpose(s_cell%plat), real(icnt,dp)
     do i = 1, s_cell%nbas
        if (any(s_relax%flags(i,:))) then
           if (s_cell%species(i)==1) then
              write(xyzu,'(a,2x,3(f12.6,2x))') 'Mg', &
                   s_cell%alat*s_cell%basis(i,:)
           else if (s_cell%species(i)==2) then
              write(xyzu,'(a,2x,3(f12.6,2x))') 'Ca', &
                   s_cell%alat*s_cell%basis(i,:)
           end if
        else
           write(xyzu,'(a,2x,3(f12.6,2x))') 'In', &
                s_cell%alat*s_cell%basis(i,:)
        end if
     end do
     
     if (fmax < s_relax%fbreak .or. .not. s_relax%do) exit rlx

     if (vf < 0.0_dp) then
        s_relax%vels(:,:) = 0.0_dp
        s_relax%tstep = s_relax%tstep*s_relax%f(2)
        s_relax%acoef = s_relax%acoef0
     else if (vf > 0.0_dp) then
        cf = s_relax%acoef*sqrt(vv/ff)
        cv = 1 - s_relax%acoef
        s_relax%vels(:,:) = cv*s_relax%vels(:,:) + cf*s_relax%force(:,:)
        s_relax%tstep = min(s_relax%tstep*s_relax%f(1),s_relax%tmax)
        s_relax%acoef = s_relax%acoef*s_relax%f(3)
     end if 
     
     do i = 1, s_cell%nbas
        do a = 1, 3
           if (s_relax%flags(i,a)) then
              s_relax%vels(i,a) = s_relax%vels(i,a) + &
                   (s_relax%tstep/s_relax%mass)*s_relax%force(i,a)
              s_cell%basis(i,a) = s_cell%basis(i,a) + &
                   s_relax%tstep*s_relax%vels(i,a)
           end if
        end do
     end do

  end do rlx

  close(rlxu)
  close(xyzu)
  
  ! write stress tensor
  open(newunit=unit,file='stress.out')
  do mu = 1, 3
     write(unit,'(3(f24.16,x))') s_relax%pst(mu,:)
  end do
  close(unit)

  ! dump atomic positions and forces
  open(newunit=unit,file='pos.out')
  open(newunit=foru,file='force.out')
  do i = 1, s_cell%nbas
     write(unit,'(i8,3(f24.16,x))') s_cell%species(i), &
          s_cell%basis(i,:)
     if (all(s_relax%flags(i,:))) then
        write(foru,'(i6,x,3(f24.16,x))') 1, s_relax%force(i,:)
     else
        write(foru,'(i6,x,3(f24.16,x))') 0, s_relax%force(i,:)
     end if
  end do
  close(unit)
  close(foru)
  
  ! calculate elastic constants
  if (s_relax%elastic .and. .not. s_agpt%do) then
  
     fv = 0.5_dp*ivol
     sites(:,:) = 0.0_dp; site_count = 1
     s_nhbr%head(:) = 0
     do i = 1, s_cell%nbas
        do j = -rcut(1), rcut(1)
           do k = -rcut(2), rcut(2)
              do l = -rcut(3), rcut(3)
                 sites(site_count,:) = s_cell%basis(i,:) &
                      + j*s_cell%plat(1,:) &
                      + k*s_cell%plat(2,:) &
                      + l*s_cell%plat(3,:)
                 basptr(site_count) = i
                 sitex(:) = s_cell%xbasis(i,:) + real((/j,k,l/),dp)
                 ! dislocated cell fudge
                 icell = find_cell(s_nhbr,sext,sitex)
                 s_nhbr%list(site_count) = s_nhbr%head(icell)
                 s_nhbr%head(icell) = site_count
                 site_count = site_count + 1
              end do
           end do
        end do
     end do
     site_count = site_count-1

     ! make neighbours
     if (s_relax%linkl) then
        call linkl_nhbr(s_cell,s_nhbr,s_nhbr%ms,sext,sites)
     else
        call nhbr_maker(s_cell%basis,s_nhbr,1,s_cell%nbas,site_count, &
             sites,s_cell%rmax)
     end if

     do a = 1, 3
        do b = 1, 3
           do c = 1, 3
              do d = 1, 3
                 cijkl(a,b,c,d) = 0.5_dp*(s_relax%pst(a,d)*kron(b,c) + &
                      s_relax%pst(b,d)*kron(a,c) + &
                      s_relax%pst(a,c)*kron(b,d) + &
                      s_relax%pst(b,c)*kron(a,d) - &
                      2.0_dp*s_relax%pst(a,b)*kron(c,d))
              end do
           end do
        end do
     end do

     nsum = 0
     do i = 1, s_cell%nbas
        ibsp = s_cell%species(i)
        do j = 1, s_nhbr%tot(i)
           ptr = s_nhbr%ptr(nsum+j)
           jbsp = s_cell%species(basptr(ptr))
           ra = s_cell%alat*(sites(ptr,:)-s_cell%basis(i,:))
           r = norm(ra)
           call force_constants(s_io(ibsp,jbsp),r,crad,ctan,epp,v2d)
           do a = 1, 3
              do b = 1, 3
                 do c = 1, 3
                    do d = 1, 3
                       cijkl(a,b,c,d) = cijkl(a,b,c,d) + &
                            fv*crad*ra(a)*ra(b)*ra(c)*ra(d)
                    end do
                 end do
              end do
           end do
        end do
        nsum = nsum + s_nhbr%tot(i)
     end do

     ! voigt and unit conversion
     print '(x,a)', 'Elastic moduli'
     do a = 1, 6
        do b = 1, 6
           idxs = voigt(a)
           jdxs = voigt(b)
           cij(a,b) = rydbohrtogpa*cijkl(idxs(1),idxs(2),jdxs(1),jdxs(2))
        end do
        print 454, cij(a,:)
     end do
     cprime = 0.5_dp*(cij(1,1)-cij(1,2))
     bmod = (cij(1,1) + cij(1,2) + cij(1,3) + &
             cij(2,1) + cij(2,2) + cij(2,3) + &
             cij(3,1) + cij(3,2) + cij(3,3)) / 9.0_dp
     print 455, bmod, cprime
  
454 format(2x,6(f10.6,x))
455 format(x,'Bulk modulus     : ',f10.6,', Shear modulus      : ', f10.6)
  
  end if

  ! molecular dynamics
  if (s_md%do) then

     if (s_md%ensemble == 1) then
        print '(a)', ' Molecular Dynamics (NVE)'
     else if (s_md%ensemble == 2 .and. .not. s_md%ne) then
        print '(a)', ' Molecular Dynamics (Gaussian Isokinetic Thermostat)'
     else if (s_md%ensemble == 2 .and. s_md%ne) then
        print '(a)', ' Molecular Dynamics (Non-Equilibrium Gaussian Isokinetic Thermostat)'
        ! num. colours and compensating force
        s_md%ncol = 0
        do i = 1, s_cell%nbas
           if (s_md%colours(i) > 0) then
              s_md%ncol(1) = s_md%ncol(1) + 1
           else
              s_md%ncol(2) = s_md%ncol(2) + 1
           end if
        end do
        s_md%colourf(-1,0) = s_md%ncol(1)*s_md%colourf(1,0)/s_md%ncol(2)
        s_md%colourf(-1,1:3) = -s_md%colourf(1,1:3)
     else if (s_md%ensemble == 3) then
        print '(a)', ' Molecular Dynamics (NVT - Langevin)'
     else if (s_md%ensemble == 4) then
        print '(a)', ' Molecular Dynamics (NPH - BZP)'
     else if (s_md%ensemble == 5) then
        print '(a)', ' Molecular Dynamics (NPT - Langevin + BZP)'
     end if

     s_md%cmv = .true.
     allocate(s_md%vels(s_cell%nbas,3), &
              s_md%force(s_cell%nbas,3), &
              s_md%pos0(s_cell%nbas,3))

     s_md%engint = 0.0_dp
     s_md%eta = 0.0_dp
     
     inquire(file='out.rst',exist=rst)
     if (rst) then
        print '(a)', ' Restarting MD...'
        call read_rst(s_cell,s_md)
        s_cell%omega = s_cell%volume/s_cell%nbas
        s_cell%rws = (0.75_dp*s_cell%omega/pi)**third
        s_cell%alat = (s_cell%volume/det(s_cell%plat))**(third)
        !!!!!! interpolate v2
        if (.not. s_agpt%do) then
           do i = 1, s_cell%nspec
              do j = i, s_cell%nspec
                 if (s_relax%interpolate) then
                    s_io(i,j)%vol = s_cell%omega
                    s_io(i,j)%rws = s_cell%rws
                 else
                    s_io(i,j)%vol = s_io(i,j)%volumes(s_io(i,j)%vptr)
                    s_io(i,j)%rws = (0.75_dp*s_io(i,j)%vol/pi)**(third)
                 end if
                 s_io(i,j)%evol = lagrange(s_io(i,j)%vol, &
                      s_io(i,j)%evolry(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1), &
                      s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
                 s_io(i,j)%evold = dlagrange(s_io(i,j)%vol, &
                      s_io(i,j)%evolry(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1), &
                      s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
                 do k = 1, s_io(i,j)%mesh
                    s_io(i,j)%v2(k,1) = s_io(i,j)%numpot(s_io(i,j)%vptr,k,1)
                    do l = 2, 5
                       s_io(i,j)%v2(k,l) = lagrange(s_io(i,j)%vol, &
                            s_io(i,j)%numpot(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1,k,l), &
                            s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
                    end do
                 end do
                 s_io(j,i) = s_io(i,j)
              end do
           end do
        else
           ! find vptr
           vdff(1) = 200.0_dp
           do k = 1, s_agpt%nvol
              vdff(2) = abs(s_agpt%volumes(k)-s_cell%omega)
              if (vdff(2) > vdff(1)) cycle
              s_agpt%vptr = k
              vdff(1) = vdff(2)
           end do
           if (s_relax%interpolate) then
              ! interpolate d2
              s_agpt%omega = s_cell%omega
              s_agpt%rws = s_cell%rws
           else
              ! interpolate d2
              s_agpt%omega = s_agpt%volumes(s_agpt%vptr)
              s_agpt%rws = (0.75_dp*s_agpt%omega/pi)**(third)
           end if
           ! calculate d2
           do k = 1, s_agpt%mesh
              s_agpt%d2(k,1) = s_agpt%dns(s_agpt%vptr,k,1)
              s_agpt%d2(k,2) = lagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,2), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
              s_agpt%d2(k,3) = lagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,4), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
              s_agpt%d2(k,4) = lagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,5), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
              s_agpt%d2(k,5) = dlagrange(s_agpt%omega, &
                   s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,2), &
                   s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
           end do
           s_agpt%nval0 = s_agpt%zval/s_cell%omega
           s_agpt%step = (s_agpt%d2(s_agpt%mesh,1)-s_agpt%d2(1,1))/(s_agpt%mesh-1)
           ! recalculate nvala
           deallocate(s_agpt%nvals)
           call agpt_on_site(s_agpt,s_cell,s_nhbragpt)
           allocate(s_agpt%nvals(s_cell%nbas))
        end if
        !!!!!!
     else
        call init_md(s_cell,s_md)
        s_md%time = 0.0_dp
     end if

     ! zercmv tracks ke change vis. Langevin
     ! converted to a.m.u
     s_md%taut = 200.0_dp/fsec
     s_md%taup = 500.0_dp/fsec

     s_md%pext = s_md%pext/rydbohrtogpa
     
     s_md%friction = 0.5_dp*s_md%taut**(-1.0_dp)
     kbt = boltzmann*s_md%itemp
     s_md%piston = (3*s_cell%nbas-3)*kbt*s_md%taup**2

     s_md%pos0(:,:) = s_cell%basis(:,:)
     
     call kinetic(s_cell,s_md)

     s_md%tstep = s_md%tstep/fsec
     
     print 234, s_md%temp, s_md%itemp, s_md%tstep*fsec, s_md%nsteps

     open(newunit=xyzu,file='XYZ',action='write')

     do imd = 0, s_md%nsteps

        if (imd > 0) then

           ! Thermostat
           if (s_md%ensemble == 2) then
              call gaussian_isokinetic(s_cell,s_md)
           else if (s_md%ensemble == 3 .or. s_md%ensemble == 5) then
              call langevin_thermostat(s_cell,s_md)
           end if

           if (s_md%ensemble == 4 .or. s_md%ensemble == 5) then
              ! recalculate ptot
              sts = s_md%pst
              call vel_stress(s_cell,s_md)
              s_md%ptot = -third*trace(s_md%pst,3)
              s_md%pst = sts
              !!!
              vf = 0.0_dp; ff = 0.0_dp
              do i = 1, s_cell%nbas
                 ibsp = s_cell%species(i)
                 mass = s_cell%mass(ibsp)
                 vf = vf + dot_product(s_md%vels(i,:),s_md%force(i,:))
                 ff = ff + dot_product(s_md%force(i,:),s_md%force(i,:))/mass
              end do
              vf = vf/s_md%piston
              ff = ff/s_md%piston
              dt2 = 0.5_dp*s_md%tstep
              pvd = s_cell%volume*(s_md%ptot-s_md%pext)
              s_md%eta = s_md%eta + 3.0_dp*dt2*(pvd+2.0_dp*kbt)/s_md%piston &
                   + vf*dt2**2 + third*ff*dt2**3
              ! for third step
              etadt = s_md%eta*s_md%tstep
           end if
              
           do i = 1, s_cell%nbas
              ibsp = s_cell%species(i)
              ! update vels (t -> t + Dt/2)
              if (s_md%ensemble == 1 .or. &
                  s_md%ensemble == 3 .or. &
                  s_md%ensemble == 4 .or. &
                  s_md%ensemble == 5) then
                 
                 s_md%vels(i,:) = s_md%vels(i,:) &
                      + 0.5_dp*s_md%tstep*s_md%force(i,:)/s_cell%mass(ibsp)

              else if (s_md%ensemble == 2) then
                 s_md%vels(i,:) = (s_md%vels(i,:) &
                      + s_md%s*s_md%force(i,:)/s_cell%mass(ibsp))/s_md%sdot
              end if

              ! update pos (t -> t + Dt)
              if (s_md%ensemble == 1 .or. &
                  s_md%ensemble == 2 .or. &
                  s_md%ensemble == 3) then
                 
                 s_cell%basis(i,:) = s_cell%basis(i,:) &
                      + s_md%tstep*s_md%vels(i,:)/s_cell%alat
              else if (s_md%ensemble == 4 .or. s_md%ensemble == 5) then
                 s_cell%basis(i,:) = exp(etadt/s_cell%alat)*s_cell%basis(i,:) &
                      + sinh(etadt/s_cell%alat)*s_md%vels(i,:)/s_md%eta
                 s_md%vels(i,:) = exp(-etadt)*s_md%vels(i,:)
              end if
           end do
           
           if (s_md%ensemble == 4 .or. s_md%ensemble == 5) then
              s_cell%volume = exp(3.0_dp*etadt)*s_cell%volume
              s_cell%alat = (s_cell%volume/det(s_cell%plat))**(third)
              s_cell%omega = s_cell%volume/s_cell%nbas
              s_cell%rws = (0.75_dp*s_cell%omega/pi)**third
              
              !!!!!! interpolate v2
              if (.not. s_agpt%do) then
                 do i = 1, s_cell%nspec
                    do j = i, s_cell%nspec
                       if (s_relax%interpolate) then
                          s_io(i,j)%vol = s_cell%omega
                          s_io(i,j)%rws = s_cell%rws
                       else
                          s_io(i,j)%vol = s_io(i,j)%volumes(s_io(i,j)%vptr)
                          s_io(i,j)%rws = (0.75_dp*s_io(i,j)%vol/pi)**(third)
                       end if
                       s_io(i,j)%evol = lagrange(s_io(i,j)%vol, &
                            s_io(i,j)%evolry(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1), &
                            s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
                       s_io(i,j)%evold = dlagrange(s_io(i,j)%vol, &
                            s_io(i,j)%evolry(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1), &
                            s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
                       do k = 1, s_io(i,j)%mesh
                          s_io(i,j)%v2(k,1) = s_io(i,j)%numpot(s_io(i,j)%vptr,k,1)
                          do l = 2, 5
                             s_io(i,j)%v2(k,l) = lagrange(s_io(i,j)%vol, &
                                  s_io(i,j)%numpot(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1,k,l), &
                            s_io(i,j)%volumes(s_io(i,j)%vptr-2:s_io(i,j)%vptr+1),4)
                          end do
                       end do
                       s_io(j,i) = s_io(i,j)
                    end do
                 end do
              else
                 ! find vptr
                 vdff(1) = 200.0_dp
                 do k = 1, s_agpt%nvol
                    vdff(2) = abs(s_agpt%volumes(k)-s_cell%omega)
                    if (vdff(2) > vdff(1)) cycle
                    s_agpt%vptr = k
                    vdff(1) = vdff(2)
                 end do
                 if (s_relax%interpolate) then
                    ! interpolate d2
                    s_agpt%omega = s_cell%omega
                    s_agpt%rws = s_cell%rws
                 else
                    ! interpolate d2
                    s_agpt%omega = s_agpt%volumes(s_agpt%vptr)
                    s_agpt%rws = (0.75_dp*s_agpt%omega/pi)**(third)
                 end if
                 ! calculate d2
                 do k = 1, s_agpt%mesh
                    s_agpt%d2(k,1) = s_agpt%dns(s_agpt%vptr,k,1)
                    s_agpt%d2(k,2) = lagrange(s_agpt%omega, &
                         s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,2), &
                         s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
                    s_agpt%d2(k,3) = lagrange(s_agpt%omega, &
                         s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,4), &
                         s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
                    s_agpt%d2(k,4) = lagrange(s_agpt%omega, &
                         s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,5), &
                         s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
                    s_agpt%d2(k,5) = dlagrange(s_agpt%omega, &
                         s_agpt%dns(s_agpt%vptr-2:s_agpt%vptr+1,k,2), &
                         s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
                 end do
                 s_agpt%nval0 = s_agpt%zval/s_cell%omega
                 s_agpt%step = (s_agpt%d2(s_agpt%mesh,1)-s_agpt%d2(1,1))/(s_agpt%mesh-1)
                 ! recalculate nvala
                 deallocate(s_agpt%nvals)
                 call agpt_on_site(s_agpt,s_cell,s_nhbragpt)
                 allocate(s_agpt%nvals(s_cell%nbas))
              end if
              !!!!!!
           end if
           
        end if
        
        sites(:,:) = 0.0_dp; site_count = 1
        s_nhbr%head(:) = 0
        do i = 1, s_cell%nbas
           do j = -rcut(1), rcut(1)
              do k = -rcut(2), rcut(2)
                 do l = -rcut(3), rcut(3)
                    sites(site_count,:) = s_cell%basis(i,:) &
                         + j*s_cell%plat(1,:) &
                         + k*s_cell%plat(2,:) &
                         + l*s_cell%plat(3,:)
                    basptr(site_count) = i
                    sitex(:) = s_cell%xbasis(i,:) + real((/j,k,l/),dp)
                    ! dislocated cell fudge
                    icell = find_cell(s_nhbr,sext,sitex)
                    s_nhbr%list(site_count) = s_nhbr%head(icell)
                    s_nhbr%head(icell) = site_count
                    site_count = site_count + 1
                 end do
              end do
           end do
        end do
        site_count = site_count-1

        ! make neighbours
        if (s_relax%linkl) then
           call linkl_nhbr(s_cell,s_nhbr,s_nhbr%ms,sext,sites)
        else
           call nhbr_maker(s_cell%basis,s_nhbr,1,s_cell%nbas,site_count, &
                sites,s_cell%rmax)
        end if
        
        if (s_agpt%do) then
           call density(s_agpt,s_cell,s_nhbr)
           s_md%e0 = 0.0_dp
           s_md%epair = 0.0_dp
           s_md%pst = 0.0_dp
        else
           s_md%e0 = s_cell%nbas*s_io(1,1)%evol
           s_md%epair = 0.0_dp; s_md%pst = 0.0_dp
           forall (i = 1:3) s_md%pst(i,i) = s_io(1,1)%evold
        end if
        
        sf = 0.0_dp; nsum = 0; ksum = 0
        s_md%fsum = 0.0_dp
        do i = 1, s_cell%nbas
           s_md%force(i,:) = 0.0_dp
           ibsp = s_cell%species(i)
           if (s_agpt%do) then
              s_md%e0 = s_md%e0 + evol(s_io(1,1), &
                   s_agpt%nvals(i),s_agpt%zval)
              nisi = s_agpt%nvals(i)
              omegai = s_agpt%zval/nisi
              chni = -omegai/nisi
              vptri = vpointer(omegai,s_io(ibsp,ibsp))
              f0 = chni*dlagrange(omegai, &
                   s_io(ibsp,ibsp)%evolry(vptri-2:vptri+1), &
                   s_io(ibsp,ibsp)%volumes(vptri-2:vptri+1),4)
              forall (j = 1:3) s_md%pst(j,j) = s_md%pst(j,j) &
                   + f0*s_agpt%dnvala*inbas                                                 
           end if
           do j = 1, s_nhbr%tot(i)
              ptr = s_nhbr%ptr(nsum+j)
              jbsp = s_cell%species(basptr(ptr))
              ra = s_cell%alat*(sites(ptr,:)-s_cell%basis(i,:))
              r = norm(ra)
              if (s_agpt%do) then
                 call agpt_force(s_agpt,s_cell,s_io(ibsp,jbsp),s_nhbr, &
                      i,ra,epp,fij,fijk,sijab,sijkab,ptr,sites,basptr, &
                      ksum)
                 s_md%force(i,:) = s_md%force(i,:)+fij(:)+fijk(:)                           
                 s_md%epair = s_md%epair + 0.5_dp*epp
                 s_md%pst(:,:) = s_md%pst(:,:)+sijab(:,:)+sijkab(:,:)
              else
                 call force_constants(s_io(ibsp,jbsp),r,crad,ctan,epp,v2d)
                 s_md%force(i,:) = s_md%force(i,:) + ctan*ra(:)
                 s_md%epair = s_md%epair + 0.5_dp*epp
                 forall (i = 1:3) s_md%pst(i,i) = s_md%pst(i,i) &                           
                      + 0.5_dp*v2d*inbas
                 do mu = 1, 3
                    do nu = 1, 3
                       s_md%pst(mu,nu) = s_md%pst(mu,nu) &                                     
                            + 0.5_dp*ctan*ra(mu)*ra(nu)/s_cell%volume
                    end do
                 end do
                 
              end if
           end do
           qr = dot_product(s_relax%sq,s_cell%basis(i,:))
           sf = sf + cos(qr)*inbas
           nsum = nsum + s_nhbr%tot(i)
           ksum = ksum + s_nhbr%tot(i)
           s_md%fsum = s_md%fsum + third*sum(s_md%force(i,:))*inbas
        end do
        s_md%etot = s_md%e0 + s_md%epair
        
        if (imd > 0) then
           if (s_md%ensemble == 2) call gaussian_isokinetic(s_cell,s_md)

           if (s_md%ensemble == 4 .or. s_md%ensemble == 5) then
              ! recalculate ptot
              sts = s_md%pst
              call vel_stress(s_cell,s_md)
              s_md%ptot = -third*trace(s_md%pst,3)
              s_md%pst = sts
              !!!!!
              vf = 0.0_dp; ff = 0.0_dp
              do i = 1, s_cell%nbas
                 ibsp = s_cell%species(i)
                 mass = s_cell%mass(ibsp)
                 vf = vf + dot_product(s_md%vels(i,:),s_md%force(i,:))
                 ff = ff + dot_product(s_md%force(i,:),s_md%force(i,:))/mass
              end do
              vf = vf/s_md%piston
              ff = ff/s_md%piston
              dt2 = 0.5_dp*s_md%tstep
              pvd = s_cell%volume*(s_md%ptot-s_md%pext)
              s_md%eta = s_md%eta + 3.0_dp*dt2*(pvd+2.0_dp*kbt)/s_md%piston &
                   + vf*dt2**2 + third*ff*dt2**3
           end if

           do i = 1, s_cell%nbas
              ibsp = s_cell%species(i)
              ! update vels (t -> t + Dt/2)
              if (s_md%ensemble == 1 .or. &
                  s_md%ensemble == 3 .or. &
                  s_md%ensemble == 4 .or. &
                  s_md%ensemble == 5) then

                 s_md%vels(i,:) = s_md%vels(i,:) &
                      + 0.5_dp*s_md%tstep*s_md%force(i,:)/s_cell%mass(ibsp)
                 
              else if (s_md%ensemble == 2) then
                 s_md%vels(i,:) = (s_md%vels(i,:) &
                      + s_md%s*s_md%force(i,:)/s_cell%mass(ibsp))/s_md%sdot
              end if
           end do
           
           if (s_md%ensemble == 3 .or. s_md%ensemble == 5) then
              call langevin_thermostat(s_cell,s_md)
           end if
           
           s_md%time = s_md%time + s_md%tstep*fsec
           if (s_md%cmv) call zercmv(s_cell,s_md)
        end if
        
        call kinetic(s_cell,s_md)
        
        ! recalculate ptot                                                                                                                                            
        sts = s_md%pst
        call vel_stress(s_cell,s_md)
        s_md%ptot = -third*trace(s_md%pst,3)
        s_md%pst = sts
        !!!!!
        
        s_md%ham = s_md%ke + s_md%etot

        s_md%enthalpy = s_md%ham &
             - 2.0_dp*kbt*log(s_cell%volume) &
             + s_md%pext*s_cell%volume &
             + 0.5_dp*s_md%piston*s_md%eta**2
        
        if (s_md%ensemble == 4 .or. s_md%ensemble == 5) then
           s_md%conq = s_md%enthalpy + s_md%engint
        else
           s_md%conq = s_md%ham + s_md%engint
        end if
                
        if (modulo(imd,s_md%wsteps) == 0) then
           print '(x,f12.3,x,es12.4,20(x,f12.6))', &
                s_md%time, &
                s_md%fsum, &
                s_md%etot*inbas, &
                s_md%ke*inbas, &
                s_md%ham*inbas, &
                s_md%enthalpy*inbas, &
                s_md%conq*inbas, &
                s_md%temp, &
                s_md%ptot*rydbohrtogpa, &
                s_cell%omega
                
           call write_rst(s_cell,s_md)
           write(xyzu,'(i6)') s_cell%nbas
           write(xyzu,'(a,f16.3)') '# time ', s_md%time
           do i = 1, s_cell%nbas
              ibsp = s_cell%species(i)
              write(xyzu,'(i6,3(f16.10,x))') ibsp, &
                   s_cell%alat*s_cell%basis(i,:)
           end do
        end if
           
     end do
     
     close(xyzu)

234  format(x,'Initial temp.    : ',f10.6,', Target temp. (NVT) : ',f10.6, / &
            x,'Time Step (fsec) : ',f10.6,', Num. time steps    : ',i10)
     
  end if

  !wallt(2) = mpi_wtime()
  call cpu_time(wallt(2))
  
  print '(a,f16.6,x,a)', ' md.gpt complete  :', wallt(2)-wallt(1), '(sec)'
  !call mpi_finalize(ierr)
  
contains

  pure subroutine llblocks(blks,plat,rmax,sext)

    implicit none
    ! passed variables
    integer(ip), intent(in) :: sext(3)
    real(dp), intent(in) :: plat(3,3), rmax
    integer(ip), intent(out) :: blks(3)
    ! local variables
    integer(ip) :: i
    real(dp) :: cxyz(3,3), wxyz(3)

    cxyz(1,:) = cross3(plat(2,:),plat(3,:))
    cxyz(2,:) = cross3(plat(3,:),plat(1,:))
    cxyz(3,:) = cross3(plat(1,:),plat(2,:))
    
    wxyz(1) = dot_product(plat(1,:),cxyz(1,:))/norm(cxyz(1,:))
    wxyz(2) = dot_product(plat(2,:),cxyz(2,:))/norm(cxyz(2,:))
    wxyz(3) = dot_product(plat(3,:),cxyz(3,:))/norm(cxyz(3,:))
    
    do i = 1, 3
       blks(i) = int(sext(i)*wxyz(i)/rmax,ip)
       if (blks(i) < 3) blks(i) = 3
    end do
    
  end subroutine llblocks
  
  pure subroutine density(s_agpt,s_cell,s_nhbr)
  
    implicit none
    ! passed variables
    type(agpt), intent(inout) :: s_agpt
    type(cell), intent(in) :: s_cell
    type(nhbr), intent(inout) :: s_nhbr
    ! local variables
    integer(ip) :: i, j, k, ibsp, jbsp, ptr, rptr, nsum
    real(dp) :: ra(3), r, npa
    
    ! calculate density
    nsum = 0
    s_agpt%nvals(:) = s_agpt%nvala
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       do j = 1, s_nhbr%tot(i)
          ptr = s_nhbr%ptr(j+nsum)
          jbsp = s_cell%species(basptr(ptr))
          ra = s_cell%alat*(sites(ptr,:)-s_cell%basis(i,:))
          r = norm(ra)
          rptr = int((r-s_agpt%d2(1,1))/s_agpt%step,ip)+1
          npa = lagrange(r,s_agpt%d2(rptr-1:rptr+2,2), &                                    
               s_agpt%d2(rptr-1:rptr+2,1),4)                                                
          s_agpt%nvals(i) = s_agpt%nvals(i) + npa
       end do
       nsum = nsum+s_nhbr%tot(i)
    end do
    
  end subroutine density

  subroutine agpt_force(s_agpt,s_cell,s_io,s_nhbr,idxi,rija,epp, &
       fij,fijk,sijab,sijkab,jptr,sites,basptr,ksum)
  
    implicit none
    ! passed variables
    type(agpt), intent(in) :: s_agpt
    type(cell), intent(in) :: s_cell
    type(io), intent(in) :: s_io
    type(nhbr), intent(in) :: s_nhbr
    integer(ip), intent(in) :: idxi, jptr, basptr(:), ksum
    real(dp), intent(in) :: rija(3), sites(:,:)
    real(dp), intent(out) :: epp, fij(3), fijk(3), sijab(3,3), &
         sijkab(3,3)
    ! local variables
    integer(ip) :: k, l, m, n, o, vptri, vptrj, vptrij, rptrij, rdptrij, &
         kptr, kbsp, vptrjk, rptrjk, rdptrik, rdptrjk, jbsp
    real(dp) :: rij, nisi, nisj, omegai, omegaj, chni, chnj, f02n, &
         nisij, omegaij, chnij, pots(4), dpots(4), npots(4), nbij, &
         wija(3), fv2r, fv2nij, srmax, rika(3), rjka(3), rik, &
         rjk, wika(3), wjka(3), nisk, nisjk, omegajk, chnjk, nbik, &
         nbjk, fv2njk, dnpaij, dnpaik, dnpajk
    
    rij = norm(rija)
    wija(:) = rija(:)/rij
    ! set up i and j
    nisi = s_agpt%nvals(idxi)
    jbsp = basptr(jptr)
    nisj = s_agpt%nvals(jbsp)
    omegai = s_agpt%zval/nisi
    omegaj = s_agpt%zval/nisj
    chni = -omegai/nisi
    chnj = -omegaj/nisj
    vptri = vpointer(omegai,s_io)
    vptrj = vpointer(omegaj,s_io)
    ! evol derivatives
    f02n = chni*dlagrange(omegai,s_io%evolry(vptri-2:vptri+1), &
         s_io%volumes(vptri-2:vptri+1),4) &
         + chnj*dlagrange(omegaj,s_io%evolry(vptrj-2:vptrj+1), &
         s_io%volumes(vptrj-2:vptrj+1),4)
    ! ij combination
    nisij = 0.5_dp*(nisi+nisj)
    omegaij = s_agpt%zval/nisij
    vptrij = vpointer(omegaij,s_io)
    chnij = -omegaij/nisij
    ! volume interpolation
    rptrij = int((rij-s_io%numpot(vptrij,1,1))/s_io%step,ip)+1
    do l = 1, 4
       pots(l) = lagrange(omegaij,s_io%numpot(vptrij-2:vptrij+1,rptrij+l-2,2), &
            s_io%volumes(vptrij-2:vptrij+1),4)
       dpots(l) = lagrange(omegaij,s_io%numpot(vptrij-2:vptrij+1,rptrij+l-2,3), &
            s_io%volumes(vptrij-2:vptrij+1),4)
       npots(l) = lagrange(omegaij,s_io%numpot(vptrij-2:vptrij+1,rptrij+l-2,5), &
            s_io%volumes(vptrij-2:vptrij+1),4)
    end do
    rdptrij = int((rij-s_agpt%d2(1,1))/s_agpt%step,ip)+1
    nbij = lagrange(rij,s_agpt%d2(rdptrij-1:rdptrij+2,3), &
         s_agpt%d2(rdptrij-1:rdptrij+2,1),4)
    epp = lagrange(rij,pots,s_io%numpot(vptrij,rptrij-1:rptrij+2,1),4)
    fv2r = lagrange(rij,dpots,s_io%numpot(vptrij,rptrij-1:rptrij+2,1),4)
    fv2nij = chnij*lagrange(rij,npots,s_io%numpot(vptrij,rptrij-1:rptrij+2,1),4)
    fij(:) = fv2r*rija(:) + (f02n + fv2nij)*nbij*wija(:)
    dnpaij = lagrange(rij, &
         s_agpt%d2(rdptrij-1:rdptrij+2,5), &
         s_agpt%d2(rdptrij-1:rdptrij+2,1),4)
    ! virial stress
    do l = 1, 3
       do m = 1, 3
          sijab(l,m) = 0.5_dp*fij(l)*rija(m)*ivol &
               + 0.5_dp*(f02n+fv2nij)*dnpaij*kron(l,m)*inbas &
               + 0.5_dp*(fv2nij)*s_agpt%dnvala*kron(l,m)*inbas
       end do
    end do
    ! three body
    fijk(:) = 0.0_dp; sijkab(:,:) = 0.0_dp
    srmax = s_cell%rmax*s_cell%alat
    do k = 1, s_nhbr%tot(idxi)
       kptr = s_nhbr%ptr(k+ksum)
       kbsp = basptr(kptr)
       rika = s_cell%alat*(sites(kptr,:)-s_cell%basis(i,:))
       rjka = s_cell%alat*(sites(kptr,:)-sites(jptr,:))
       rik = norm(rika)
       rjk = norm(rjka)
       if (rjk < 1e-4_dp .or. rjk > srmax) cycle
       wika(:) = rika(:)/rik
       wjka(:) = rjka(:)/rjk
       nisk = s_agpt%nvals(kbsp)
       nisjk = 0.5_dp*(nisj+nisk)
       omegajk = s_agpt%zval/nisjk
       vptrjk = vpointer(omegajk,s_io)
       chnjk = -omegajk/nisjk
       rptrjk = int((rjk-s_io%numpot(vptrjk,1,1))/s_io%step,ip)+1
       rdptrik = int((rik-s_agpt%d2(1,1))/s_agpt%step,ip)+1
       rdptrjk = int((rjk-s_agpt%d2(1,1))/s_agpt%step,ip)+1
       nbik = lagrange(rik,s_agpt%d2(rdptrik-1:rdptrik+2,3), &
            s_agpt%d2(rdptrik-1:rdptrik+2,1),4)
       nbjk = lagrange(rjk,s_agpt%d2(rdptrjk-1:rdptrjk+2,3), &
            s_agpt%d2(rdptrjk-1:rdptrjk+2,1),4)
       do n = 1, 4
          npots(n) = lagrange(omegajk, &
               s_io%numpot(vptrjk-2:vptrjk+1,rptrjk+n-2,5), &
               s_io%volumes(vptrjk-2:vptrjk+1),4)
       end do
       fv2njk = chnjk*lagrange(rjk,npots, &
            s_io%numpot(vptrjk,rptrjk-1:rptrjk+2,1),4)
       fijk(:) = fijk(:) + &
            0.5_dp*fv2nij*nbik*wika + &
            0.25_dp*fv2njk*(nbij*wija+nbik*wika)
       dnpaik = lagrange(rik, &
            s_agpt%d2(rdptrik-1:rdptrik+2,5), &
            s_agpt%d2(rdptrik-1:rdptrik+2,1),4)
       dnpajk = lagrange(rjk, &
            s_agpt%d2(rdptrjk-1:rdptrjk+2,5), &
            s_agpt%d2(rdptrjk-1:rdptrjk+2,1),4)
       do n = 1, 3
          do o = 1, 3
             sijkab(n,o) = sijkab(n,o) &
                  + 0.25_dp*ivol*fv2nij*(nbik*wika(n)*rika(o)+nbjk*wjka(n)*rjka(o)) &
                  + 0.25_dp*fv2nij*kron(n,o)*(dnpaik+dnpajk)*inbas
          end do
       end do
    end do
    
  end subroutine agpt_force

  subroutine agpt_on_site(s_agpt,s_cell,s_nhbr)
  
    implicit none
    ! passed variables
    type(agpt), intent(inout) :: s_agpt
    type(cell), intent(in) :: s_cell
    type(nhbr), intent(inout) :: s_nhbr
    ! local variables
    integer(ip) :: i, j, k, l, ptr, nmax, site_count, ibsp, jbsp, &
         rptr, n, nsum 
    integer(ip), allocatable :: basptr(:)
    real(dp) :: nval0, npa, ra(3), alat, omega, rws, r, step, &
         nvalas, iomega, alats(4), rcut(3)
    real(dp), allocatable :: sites(:,:), nvala(:)
  
    allocate(sites(smax,3),basptr(smax))
    nmax = int(650*s_agpt%nbas,ip)
    allocate(s_nhbr%ptr(nmax),s_nhbr%tot(s_agpt%nbas), &
         s_agpt%nvals(s_agpt%nbas))
    ! create site list
    do i = 1, 3
       rcut(i) = int(s_cell%rmax/norm(s_agpt%plat(i,:)),ip)+2
    end do

    sites(:,:) = 0.0_dp; site_count = 1
    do i = 1, s_agpt%nbas
       do j = -rcut(1), rcut(1)
          do k = -rcut(2), rcut(2)
             do l = -rcut(3), rcut(3)
                sites(site_count,:) = s_agpt%basis(i,:) &
                     + j*s_agpt%plat(1,:) &
                     + k*s_agpt%plat(2,:) &
                     + l*s_agpt%plat(3,:)
                basptr(site_count) = i
                site_count = site_count + 1
             end do
          end do
       end do
    end do
    site_count = site_count-1
    
    ! make neighbours
    call nhbr_maker(s_agpt%basis,s_nhbr,1,s_agpt%nbas,site_count, &
         sites,s_cell%rmax)

    if (.not. s_agpt%derivative) then
       s_agpt%nvals(:) = 0.0_dp
       nsum = 0
       do i = 1, s_agpt%nbas
          s_agpt%nvals(i) = 0.0_dp
          do j = 1, s_nhbr%tot(i)
             ptr = s_nhbr%ptr(nsum+j)
             ra = s_cell%alat*(sites(ptr,:)-s_agpt%basis(i,:))
             r = norm(ra)
             rptr = int((r-s_agpt%d2(1,1))/s_agpt%step,ip)+1
             npa = lagrange(r,s_agpt%d2(rptr-1:rptr+2,2), &
                  s_agpt%d2(rptr-1:rptr+2,1),4)
             s_agpt%nvals(i) = s_agpt%nvals(i) + npa
          end do
          nsum = nsum+s_nhbr%tot(i)
       end do
       s_agpt%nvala = sum(s_agpt%nval0-s_agpt%nvals)/s_agpt%nbas
       s_agpt%dnvala = 0.0_dp
    else
       allocate(nvala(4))
       do n = -2, 1
          omega = s_agpt%volumes(s_agpt%vptr+n)
          nval0 = s_agpt%zval/omega
          alats(3+n) = (s_cell%nbas*omega/det(s_cell%plat))**(third)
          nsum = 0
          do i = 1, s_agpt%nbas
             ibsp = s_agpt%species(i)
             s_agpt%nvals(i) = 0.0_dp
             do j = 1, s_nhbr%tot(i)
                ptr = s_nhbr%ptr(nsum+j)
                ra = alats(3+n)*(sites(ptr,:)-s_agpt%basis(i,:))
                r = norm(ra)
                rptr = int((r-s_agpt%dns(s_agpt%vptr+n,1,1))/s_agpt%step,ip)+1
                npa = lagrange(r, &
                     s_agpt%dns(s_agpt%vptr+n,rptr-1:rptr+2,2), &
                     s_agpt%dns(s_agpt%vptr+n,rptr-1:rptr+2,1),4)
                s_agpt%nvals(i) = s_agpt%nvals(i) + npa
             end do
             nsum = nsum+s_nhbr%tot(i)
          end do
          nvala(3+n) = sum(nval0-s_agpt%nvals(:))/s_agpt%nbas
       end do
       s_agpt%nvala = lagrange(s_agpt%omega,nvala, &
            s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
       s_agpt%dnvala = dlagrange(s_agpt%omega,nvala, &
            s_agpt%volumes(s_agpt%vptr-2:s_agpt%vptr+1),4)
       deallocate(nvala)
    end if
    
    deallocate(s_nhbr%ptr,s_nhbr%tot,s_agpt%nvals)
  
  end subroutine agpt_on_site

  subroutine write_rst(s_cell,s_md)

    implicit none
    ! passed variables
    type(cell), intent(inout) :: s_cell
    type(md), intent(inout) :: s_md
    ! local variables
    integer(ip) :: i, unit
    
    open(newunit=unit,file='out.rst',action='write')
    write(unit,'(f20.12)') s_md%time, s_md%engint, s_md%eta, s_cell%volume
    do i = 1, s_cell%nbas
       write(unit,'(6(f16.12,x))') s_cell%basis(i,:), s_md%vels(i,:)
    end do
    close(unit)
    
  end subroutine write_rst

  subroutine read_rst(s_cell,s_md)

    implicit none
    ! passed variables
    type(cell), intent(inout) :: s_cell
    type(md), intent(inout) :: s_md
    ! local variables
    integer(ip) :: i, unit
    
    open(newunit=unit,file='out.rst',action='read')
    read(unit,*) s_md%time, s_md%engint, s_md%eta, s_cell%volume
    do i = 1, s_cell%nbas
       read(unit,*) s_cell%basis(i,:), s_md%vels(i,:)
    end do
    close(unit)
    
  end subroutine read_rst

  subroutine read_density(s_agpt)

    implicit none
    ! passed variables
    type(agpt), intent(inout) :: s_agpt
    ! local variables
    integer(ip) :: unit, i, j
    real(dp) :: rws

    open(newunit=unit,file=trim(s_agpt%dnfile),action='read')
    read(unit,*) s_agpt%nvol, s_agpt%mesh
    allocate(s_agpt%dns(s_agpt%nvol,s_agpt%mesh,5), &
         s_agpt%volumes(s_agpt%nvol))
    do i = 1, s_agpt%nvol
       read(unit,*) s_agpt%volumes(i), s_agpt%ra
       rws = (0.75_dp*s_agpt%volumes(i)/pi)**(third)
       do j = 1, s_agpt%mesh
          read(unit,*) s_agpt%dns(i,j,:)
       end do
    end do
    close(unit)

  end subroutine read_density
    
  subroutine gaussian_isokinetic(s_cell,s_md)

    implicit none
    ! passed variables
    type(cell), intent(in) :: s_cell
    type(md), intent(inout) :: s_md
    ! local variables
    integer(ip) :: i, ibsp
    real(dp) :: alpha, fdp, fdf, mass, dt2, twoke, a, b, sqrtb

    fdp = 0.0_dp; fdf = 0.0_dp; twoke = 0.0_dp
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       mass = s_cell%mass(ibsp)
       fdp = fdp + dot_product(s_md%force(i,:),s_md%vels(i,:))
       fdf = fdf + dot_product(s_md%force(i,:),s_md%force(i,:))/mass
       twoke = twoke + mass*dot_product(s_md%vels(i,:),s_md%vels(i,:))
    end do
    a = fdp/twoke
    b = fdf/twoke
    dt2 = 0.5_dp*s_md%tstep
    sqrtb = sqrt(b)
    s_md%s = (a/b)*(cosh(dt2*sqrtb)-1.0_dp) + sinh(dt2*sqrtb)/sqrtb
    s_md%sdot = (a/sqrtb)*sinh(dt2*sqrtb) + cosh(dt2*sqrtb)

  end subroutine gaussian_isokinetic
  
  pure subroutine force_constants(s_io,r,crad,ctan,epp,v2d)

    implicit none
    ! Passed Variables
    type(io), intent(in) :: s_io
    real(dp), intent(in) :: r
    real(dp), intent(out) :: ctan, crad, epp, v2d
    ! Local Variables
    integer(ip) :: ptr
    real(dp) :: pdd

    ptr = int((r-s_io%v2(1,1))/s_io%step,ip)+1
    ctan = lagrange(r,s_io%v2(ptr-1:ptr+2,3),s_io%v2(ptr-1:ptr+2,1),4)
    pdd = lagrange(r,s_io%v2(ptr-1:ptr+2,4),s_io%v2(ptr-1:ptr+2,1),4)
    epp = lagrange(r,s_io%v2(ptr-1:ptr+2,2),s_io%v2(ptr-1:ptr+2,1),4)
    v2d = lagrange(r,s_io%v2(ptr-1:ptr+2,5),s_io%v2(ptr-1:ptr+2,1),4)
    !v2dd = lagrange(r,s_io%v2(ptr-1:ptr+2,6),s_io%v2(ptr-1:ptr+2,1),4)
    !dv2d = lagrange(r,s_io%v2(ptr-1:ptr+2,7),s_io%v2(ptr-1:ptr+2,1),4)
    crad = (pdd - ctan)/(r**2)
    
  end subroutine force_constants

  subroutine init_md(s_cell, s_md)

    implicit none
    ! Passed Variables
    type(cell), intent(in) :: s_cell
    type(md), intent(inout) :: s_md
    ! Local Variables
    integer(ip) :: i, j, ibsp
    real(dp) :: sigma
    
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       sigma = sqrt(boltzmann*s_md%itemp/s_cell%mass(ibsp))
       do j = 1, 3
          s_md%vels(i,j) = sigma*random_gauss()
       end do
    end do
    
    call zercmv(s_cell,s_md)

  end subroutine init_md

  subroutine kinetic(s_cell, s_md)

    implicit none
    ! Passed Variables
    type(cell), intent(in) :: s_cell
    type(md), intent(inout) :: s_md
    ! Local Variables
    integer(ip) :: i, ibsp
    real(dp) :: mass

    s_md%ke = 0.0_dp
    s_md%vsum = 0.0_dp
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       mass = s_cell%mass(ibsp)
       s_md%ke = s_md%ke + 0.5_dp*mass*dot_product(s_md%vels(i,:),s_md%vels(i,:))
       s_md%vsum = s_md%vsum + s_md%vels(i,:)
    end do
    
    if (s_md%cmv) then
       s_md%temp = 2.0_dp*s_md%ke/((3*s_cell%nbas-3)*boltzmann)
    else
       s_md%temp = 2.0_dp*third*s_md%ke/(s_cell%nbas*boltzmann)
    end if

  end subroutine kinetic

  subroutine vel_stress(s_cell,s_md)

    implicit none
    ! passed variables
    type(cell), intent(in) :: s_cell
    type(md), intent(inout) :: s_md
    ! local variables
    integer(ip) :: i, mu, nu
    real(dp) :: ivol, mass

    ivol = 1.0_dp/s_cell%volume
    
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       mass = s_cell%mass(ibsp)
       do mu = 1, 3
          do nu = 1, 3
             s_md%pst(mu,nu) = s_md%pst(mu,nu) + &
                  ivol*mass*s_md%vels(i,mu)*s_md%vels(i,nu)
          end do
       end do
    end do

  end subroutine vel_stress
  
  subroutine zercmv(s_cell, s_md)

    implicit none
    ! Passed Variables
    type(cell), intent(in) :: s_cell
    type(md), intent(inout) :: s_md
    ! Local Variable
    integer(ip) :: i, ibsp
    real(dp) :: sum_vels(3), v2

    sum_vels = sum(s_md%vels,dim=1)*inbas
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       v2 = dot_product(s_md%vels(i,:),s_md%vels(i,:))
       s_md%engint = s_md%engint + 0.5_dp*s_cell%mass(ibsp)*v2 
       s_md%vels(i,:) = s_md%vels(i,:) - sum_vels(:)
       v2 = dot_product(s_md%vels(i,:),s_md%vels(i,:))
       s_md%engint = s_md%engint - 0.5_dp*s_cell%mass(ibsp)*v2
    end do

  end subroutine zercmv

  pure subroutine linkl_nhbr(s_cell,s_nhbr,blks,sext,sites)

    implicit none
    ! passed variables
    type(cell), intent(in) :: s_cell
    type(nhbr), intent(inout) :: s_nhbr
    integer(ip), intent(in) :: blks(3), sext(3)
    real(dp), intent(in) :: sites(:,:)
    ! local variables
    integer(ip) :: i, icell, imx, imy, imz, j, jcell, nsum
    real(dp) :: bx(3), r

    nsum = 1 
    do i = 1, s_cell%nbas
       s_nhbr%tot(i) = 0
       bx = matmul(s_cell%basis(i,:),s_cell%iplat)
       icell = find_cell(s_nhbr,sext,bx)
       do imx = -1, 1
          do imy = -blks(1), blks(1), blks(1)
             do imz = -blks(1)*blks(2), blks(1)*blks(2), blks(1)*blks(2)
                jcell = icell+imx+imy+imz
                j = s_nhbr%head(jcell)
                do while (j /= 0)
                   r = norm(sites(j,:)-s_cell%basis(i,:))
                   if (r <= s_cell%rmax .and. r > 0.0_dp) then
                      s_nhbr%ptr(nsum) = j
                      s_nhbr%tot(i) = s_nhbr%tot(i)+1
                      nsum = nsum + 1
                   end if
                   j = s_nhbr%list(j)
                end do
             end do
          end do
       end do
    end do
    
  end subroutine linkl_nhbr
  
  pure subroutine nhbr_maker(basis,s_nhbr,start,end,site_count,sites,rmax)

    implicit none
    ! Passed Variables
    integer(ip), intent(in) :: site_count, start, end
    type(nhbr), intent(inout) :: s_nhbr
    real(dp), intent(in) :: basis(:,:), sites(:,:), rmax
    ! Local Variables
    integer(ip) :: i, j, nsum
    real(dp) :: rij, rija(3)

    nsum = 1
    do i = start, end
       s_nhbr%tot(i) = 0
       do j = 1, site_count
          rija(:) = sites(j,:)-basis(i,:)
          rij = norm(rija)
          if (rij <= rmax .and. rij > 0.0_dp) then
             s_nhbr%ptr(nsum) = j
             s_nhbr%tot(i) = s_nhbr%tot(i)+1
             nsum = nsum+1
          end if
       end do
    end do
    
  end subroutine nhbr_maker

  subroutine read_potential(filename, s_io)

    implicit none
    ! Passed Variables
    character(*), intent(in) :: filename
    type(io), intent(out) :: s_io
    ! Local Variables
    integer(ip) :: mode, i, n, unit
    real(dp) :: unitr

    open(newunit=unit,file=trim(filename),action='read')

    read(unit,*) s_io%nza, unitr, s_io%mesh, s_io%nvol, mode, s_io%pn

    allocate(s_io%numpot(s_io%nvol,s_io%mesh,5), &
         s_io%volumes(s_io%nvol), &
         s_io%evolry(s_io%nvol))

    do n = 1, s_io%nvol
       read(unit,*) s_io%volumes(n), s_io%evolry(n), s_io%mass
       read(unit,*) s_io%p1, s_io%r0rws, s_io%al, s_io%va, s_io%vb
       do i = 1, s_io%mesh
          read(unit,*) s_io%numpot(n,i,1:4)
       end do
    end do

    close(unit)
    
    s_io%p2 = 2.0_dp*s_io%p1
    s_io%p4 = 2.0_dp*s_io%p2

  end subroutine read_potential

  subroutine make_cell(s_agpt, s_cell, s_md, s_relax)

    implicit none
    ! passed variables
    type(agpt), intent(out) :: s_agpt
    type(cell), intent(out) :: s_cell
    type(md), intent(out) :: s_md
    type(relax), intent(out) :: s_relax
    ! local variables
    character(len=128) :: line, tmp_filename
    integer(ip) :: unit, sptot, i, j, k

    ! defaults
    s_agpt%do = .false.
    s_agpt%derivative = .false.
    s_relax%do = .false.
    s_relax%elastic = .false.
    s_relax%screen = .false.
    s_relax%interpolate = .true.
    s_relax%linkl = .false.
    !s_nhbr%ms = (/3,3,3/)
    s_md%do = .false.
    s_md%ne = .false.
    s_md%pext = 0.0_dp
    s_md%wsteps = 1
    
    open(newunit=unit,file='CONTROL',action='read')
10  continue
    read(unit,*,end=30) line
    if (scan(trim(line),'#') == 1) goto 10
    select case(trim(line))
    
    case('ALAT')
       read(unit,*) s_cell%alat
    case('AGPT')
       read(unit,*) s_agpt%do
    case('AGPT_BASIS')
       if (.not. allocated(s_agpt%basis) .or. &
            .not. allocated(s_agpt%species)) then
          print 500
          call exit(1)
       end if
       do i = 1, s_agpt%nbas
          read(unit,*) s_agpt%species(i), s_agpt%basis(i,:)
       end do
    case('AGPT_DENSITY')
       read(unit,*) s_agpt%dnfile
    case('AGPT_DERIVATIVE')
       read(unit,*) s_agpt%derivative
    case('AGPT_NBAS')
       read(unit,*) s_agpt%nbas
       allocate(s_agpt%basis(s_agpt%nbas,3), &
            s_agpt%species(s_agpt%nbas))
    case('AGPT_PLAT')
       do i = 1, 3
          read(unit,*) s_agpt%plat(i,:)
       end do
    case('AGPT_ZVAL')
       read(unit,*) s_agpt%zval
    case('BASIS')
       if (.not. allocated(s_cell%basis) .or. &
            .not. allocated(s_cell%species) .or. &
            .not. allocated(s_relax%flags) .or. &
            .not. allocated(s_md%colours)) then
          print 500
          call exit(1)
       end if
       do i = 1, s_cell%nbas
          if (s_md%ne) then
             read(unit,*) s_cell%species(i), s_cell%basis(i,:), &
                  s_relax%flags(i,:), s_md%colours(i)
          else
             read(unit,*) s_cell%species(i), s_cell%basis(i,:), &
                  s_relax%flags(i,:)
          end if
       end do
    case('ELASTIC')
       read(unit,*) s_relax%elastic
    case('MD_ENSEMBLE')
       read(unit,*) s_md%ensemble
    case('HARD_CUT')
       read(unit,*) s_cell%rmaxrws
    case('INTERPOLATE')
       read(unit,*) s_relax%interpolate
    case('LINKED_LIST')
       read(unit,*) s_relax%linkl
    !case('LL_BLOCKS')
    !   read(unit,*) s_nhbr%ms(:)
    case('MASSES')
       if (.not. allocated(s_cell%mass)) then
          print 500
          call exit(1)
       end if
       do i = 1, s_cell%nspec
          read(unit,*) s_cell%mass(i)
          s_cell%mass(i) = s_cell%mass(i)/amass
       end do
    case ('MD')
       read(unit,*) s_md%do
    case ('MD_COLOUR_FORCE')
       read(unit,*) s_md%colourf(1,:)
    case('MD_ITEMP')
       read(unit,*) s_md%itemp
    case('MD_NSTEPS')
       read(unit,*) s_md%nsteps
    case('MD_NE')
       read(unit,*) s_md%ne
    case('MD_PRESSURE')
       read(unit,*) s_md%pext
    case('MD_TSTEP')
       read(unit,*) s_md%tstep
    case('MD_WSTEPS')
       read(unit,*) s_md%wsteps
    case('NBAS')
       read(unit,*) s_cell%nbas
       allocate(s_cell%basis(s_cell%nbas,3), &
            s_cell%species(s_cell%nbas), &
            s_relax%flags(s_cell%nbas,3), &
            s_md%colours(s_cell%nbas))
       s_relax%flags(:,:) = .true.
    case('NSPEC')
       read(unit,*) s_cell%nspec
       allocate(s_cell%filename(s_cell%nspec,s_cell%nspec), &
            s_cell%mass(s_cell%nspec))
    case('POTENTIAL')
       if (.not. allocated(s_cell%filename)) then
          print 500
          call exit(1)
       end if
       sptot = triangle(s_cell%nspec)
       do i = 1, sptot
          read(unit,*) j, k, tmp_filename
          s_cell%filename(j,k) = tmp_filename
       end do
    case('PLAT')
       do i = 1, 3
          read(unit,*) s_cell%plat(i,:)
       end do
    case('RELAX')
       read(unit,*) s_relax%do
    case('RELAX_BREAKFORCE')
       read(unit,*) s_relax%fbreak
    case('SCREEN')
       read(unit,*) s_relax%screen

    end select
    goto 10
    
30  continue
    close(unit)

500 format('EXIT 1 : INCORRECT ORDER IN CONTROL')

  end subroutine make_cell

  subroutine init_random_seed()

    implicit none
    ! Local Variables
    integer(ip) :: i, n, clock
    integer(ip), allocatable :: seed(:)

    call random_seed(size = n)
    allocate(seed(n))
    
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)

  end subroutine init_random_seed

  subroutine langevin_thermostat(s_cell,s_md)

    implicit none
    ! passed variables
    type(cell), intent(in) :: s_cell
    type(md), intent(inout) :: s_md
    ! local variables
    integer(ip) :: i, j, ibsp
    real(dp) :: c(3), dt2, kbt, keij

    dt2 = 0.5_dp*s_md%tstep
    kbt = boltzmann*s_md%itemp
    c(1) = exp(-s_md%friction*dt2)
    do i = 1, s_cell%nbas
       ibsp = s_cell%species(i)
       c(2) = sqrt((1.0_dp-c(1)**2)*kbt/s_cell%mass(ibsp))
       do j = 1, 3
          keij = 0.5_dp*s_cell%mass(ibsp)*s_md%vels(i,j)**2
          s_md%engint = s_md%engint + keij
          s_md%vels(i,j) = c(1)*s_md%vels(i,j) + c(2)*random_gauss()
          keij = 0.5_dp*s_cell%mass(ibsp)*s_md%vels(i,j)**2
          s_md%engint = s_md%engint - keij
       end do
    end do
    
    if (s_md%ensemble == 5) then 
       c(3) = sqrt((1.0_dp-c(1)**2)*kbt/s_md%piston)
       keij = 0.5_dp*s_md%piston*s_md%eta**2
       s_md%engint = s_md%engint + keij
       s_md%eta = c(1)*s_md%eta + c(3)*random_gauss()
       keij = 0.5_dp*s_md%piston*s_md%eta**2
       s_md%engint = s_md%engint - keij
    end if
    
  end subroutine langevin_thermostat

  integer(ip) pure function find_cell(s_nhbr,sext,site)

    implicit none
    ! passed variables
    type(nhbr), intent(in) :: s_nhbr
    integer(ip), intent(in) :: sext(3)
    real(dp), intent(in) :: site(3)
    ! local variables
    integer(ip) :: i, ixs(3)

    do i = 1, 3
       ixs(i) = int((site(i)+s_nhbr%cent(i))/sext(i)*s_nhbr%ms(i))
       if (ixs(i) >= s_nhbr%ms(i)) ixs(i) = s_nhbr%ms(i)-1 
    end do
    find_cell = 1 + ixs(1) + s_nhbr%ms(1)*ixs(2) &
         + s_nhbr%ms(2)*s_nhbr%ms(1)*ixs(3)
        
  end function find_cell
    
  integer(ip) pure function triangle(trm) 
    
    implicit none
    ! Passed Variables
    integer(ip), intent(in) :: trm
    ! Local Variables
    integer(ip) :: i

    triangle = 0
    do i = 1, trm
       triangle = triangle + i
    end do
    
  end function triangle

  real(dp) function random_gauss()

    implicit none
    ! Local Variable
    integer(ip), save :: iset = 0
    real(dp) :: uniform_r(2), chi
    real(dp), save :: gset
    
    if (iset == 0) then
       call random_number(uniform_r)
       chi = sqrt(-2*log(uniform_r(1)))
       random_gauss = chi*cos(2*pi*uniform_r(2))
       gset = chi*sin(2*pi*uniform_r(2))
       iset = 1
    else
       random_gauss = gset
       iset = 0
    end if
    
  end function random_gauss
  
  real(dp) function evol(s_io,ni,zval)

    implicit none
    ! passed variables
    type(io), intent(in) :: s_io
    real(dp), intent(in) :: ni, zval
    ! local variables
    integer(ip) ::  ptr
    real(dp) :: omega
    
    omega = zval/ni
    ptr = vpointer(omega,s_io)
    evol = lagrange(omega,s_io%evolry(ptr-2:ptr+1), &
         s_io%volumes(ptr-2:ptr+1),4)
    
  end function evol

  integer(ip) pure function vpointer(omega,s_io)

    implicit none
    ! passed variables
    type(io), intent(in) :: s_io
    real(dp), intent(in) :: omega
    ! local variables
    integer(ip) :: vp, idx
    real(dp) :: x

    x = (omega/s_io%volumes(s_io%vptr))**third
    idx = nint((1.0_dp-x)/s_io%vstep)
    vpointer = s_io%vptr+idx
    
    if (vpointer <= 2) then
       vpointer = 3
    else if (vpointer > s_io%nvol-1) then
       vpointer = s_io%nvol-1
    end if
    
  end function vpointer

  pure function cross3(a,b) result (c)

    implicit none
    ! passed variables
    real(dp), intent(in) :: a(3), b(3)
    ! local variables
    real(dp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    
  end function cross3
  
  pure function inv3(matrix) result(imat)

    implicit none
    ! Passed Variables
    real(dp), intent(in) :: matrix(3,3) ! Must be 3x3 matrix
    ! Local Variables
    real(dp) :: detm, imat(3,3)

    detm = det(matrix)

    imat(1,1) = matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)
    imat(1,2) = matrix(1,3)*matrix(3,2)-matrix(1,2)*matrix(3,3)
    imat(1,3) = matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2)

    imat(2,1) = matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3)
    imat(2,2) = matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1)
    imat(2,3) = matrix(1,3)*matrix(2,1)-matrix(1,1)*matrix(2,3)
    
    imat(3,1) = matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1)
    imat(3,2) = matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2)
    imat(3,3) = matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)

    imat = imat/detm
    
  end function inv3
  
  real(dp) pure function det(matrix)

    implicit none
    ! Passed Variables
    real(dp), intent(in) :: matrix(3,3) ! Must be 3x3 matrix

    det = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
         - matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
         + matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
    
  end function det

  real(dp) pure function fl(r, s_io)

    implicit none
    ! Passed Variables
    type(io), intent(in) :: s_io
    real(dp), intent(in) :: r
    ! Local Variables
    real(dp) :: term, quan

    term = (s_io%rp/r)**(s_io%p1)

    if (r < s_io%r0) then
       fl = term
    else
       quan = s_io%al*(r/s_io%r0-1.0_dp)**2
       fl = term*(1.0_dp+quan)*exp(-quan)
    end if

  end function fl

  real(dp) pure function hgauss(r, s_io)

    implicit none
    ! Passed Variables
    type(io), intent(in) :: s_io
    real(dp), intent(in) :: r

    hgauss = (1.0_dp+s_io%al*(r-1.0_dp)**2)*exp(-s_io%al*(r-1.0_dp)**2)

  end function hgauss

  real(dp) pure function lagrange(x, f, x1, ord)

    implicit none
    ! Passed Variables
    integer(ip), intent(in) :: ord
    real(dp), intent(in) :: f(ord), x1(ord), x
    ! Local Variables
    integer(ip) :: i, j
    real(dp) :: f1, f2

    lagrange = 0.0_dp
    do i = 1, ord
       f1 = f(i)
       f2 = 1.0_dp
       do j = 1, ord
          if (i == j) cycle
          f1 = f1*(x - x1(j))
          f2 = f2*(x1(i) - x1(j))
       end do
       lagrange = lagrange + f1/f2
    end do

  end function lagrange

  real(dp) pure function dlagrange(x, f, x1, ord)

    implicit none
    ! Passed Variables
    integer(ip), intent(in) :: ord
    real(dp), intent(in) :: f(ord), x1(ord), x
    ! Local Variables
    integer(ip) :: i, j, m
    real(dp) :: lj, pp

    dlagrange = 0.0_dp
    do j = 1, ord
       lj = 0.0_dp
       do i = 1, ord
          if (i == j) cycle
          pp = 1.0_dp
          do m = 1, ord
             if ((m == i) .or. (m == j)) cycle
             pp = pp*(x - x1(m))/(x1(j)-x1(m))
          end do
          lj = lj + pp/(x1(j)-x1(i))
       end do
       dlagrange = dlagrange + f(j)*lj
    end do

  end function dlagrange

  real(dp) pure function ddlagrange(x, f, x1, ord)

    implicit none
    ! Passed Variables
    integer(ip), intent(in) :: ord
    real(dp), intent(in) :: f(ord), x1(ord), x
    ! Local Variables
    integer(ip) :: i, j, m, l
    real(dp) :: lj, inr, pp

    ddlagrange = 0.0_dp
    do j = 1, ord
       lj = 0.0_dp
       do i = 1, ord
          if (i == j) cycle
          inr = 0.0_dp
          do m = 1, ord
             if ((m == i) .or. (m == j)) cycle
             pp = 1.0_dp
             do l = 1, ord
                if ((l == i) .or. (l == j) .or. (l == m)) cycle
                pp = pp*(x - x1(l))/(x1(j)-x1(l))
             end do
             inr = inr + pp/(x1(j)-x1(m))
          end do
          lj = lj + inr/(x1(j)-x1(i))
       end do
       ddlagrange = ddlagrange + f(j)*lj
    end do

  end function ddlagrange

  real(dp) pure function norm(vec)

    implicit none
    ! passed variables
    real(dp), intent(in) :: vec(:)

    norm = sqrt(dot_product(vec,vec))

  end function norm

  real(dp) pure function kron(a,b)
    
    implicit none
    ! passed variables
    integer(ip), intent(in) :: a, b
    
    if (a == b) then
       kron = 1.0_dp
    else
       kron = 0.0_dp
    end if
  
  end function kron

  pure function voigt(a)

    implicit none
    ! passed variables
    integer(ip), intent(in) :: a
    integer(ip) :: voigt(2)
  
    voigt(:) = 0
    if (a < 4) then
       voigt(1:2) = a
    else if (a == 4) then
       voigt(1) = 2
       voigt(2) = 3
    else if (a == 5) then
       voigt(1) = 1
       voigt(2) = 3
    else if (a == 6) then
       voigt(1) = 1
       voigt(2) = 2
    end if
    
  end function voigt

  real(dp) pure function trace(mat,n)
    
    implicit none
    ! passed variables
    integer(ip), intent(in) :: n
    real(dp), intent(in) :: mat(n,n)
    ! local variables
    integer(ip) :: i
    
    trace = 0.0_dp
    do i = 1, n
       trace = trace + mat(i,i)
    end do
    
  end function trace

end program mdgpt
