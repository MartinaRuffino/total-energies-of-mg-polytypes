      subroutine lmrefdos(mode,npr,eminr,emaxr,npf,eminf,emaxf,nchan,
     .  efermi,eshft,wt,wg,refdos,fitdos,outdos,sigmrq)
C- Shift and scale reference dos in preparation to compare against fit dos
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1s digit generates outdos
Ci         :1 make outdos by spline interpolation reference dos
Ci         :2 same as 1, modified spline interpolation
Ci         :4 Also scale area of outdos to match fitdos
Ci            Requires 10s digit to be set
Ci         :10s digit
Ci         :1 generates "std deviation" sigmrq
Ci   npr   :number of energy points, reference dos
Ci   eminr :lower bound, reference dos
Ci   emaxr :lower bound, reference dos
Ci   npf   :number of energy points, fit and out dos
Ci   eminf :lower bound, fit and out dos
Ci   emaxf :lower bound, fit and out dos
Ci   nchan :number of dos channels
Ci   efermi:Fermi energy.  efermi(1), fit dos; efermi(2), ref dos
Ci   eshft :Fermi energy shift to align reference with fit dos:
Ci         :efermi(fit) - efermi(ref)
Ci   refdos:reference dos
Ci   fitdos:fit dos
Ci   wt    :Fitting weights for DOS
Ci         :wt(1) determines functional form of weights
Ci         :wt(2) is energy shift, denoted w2 below
Ci         :wt(3) is energy scale, denoted w3 below
Ci         :wt(1):
Ci         :0 unit weights
Ci         :1 Fermi function, wt = {1+exp[-(E-(w2+Ef))/w3]}^-1
Ci         :2 Fermi-like gaussian cutoff,
Ci            wt = {1+exp[-((E_kn-(w2+Ef))/w3)^2]}^-1
Ci         :3 Exponential weights, wt = exp[-|(E_kn-(w2+Ef))/w3|]
Ci         :4 Gaussian weights, wt = exp[-((E_kn-(w2+Ef))/w3)^2]
Ci   wg    :Reference DOS is broadened by a gaussian of width wg
Co Outputs
Co   outdos:modified reference dos.  Modified as follows:
Co         :(1) shifted by eshft
Co         :(2) spline modified to minimize curvature in gap regions
Co         :    (2's part of 1's digit mode nonzero)
Co         :(3) interpolated to fit mesh
Co         :(4) area scaled to match that of fit dos
Co         :    (4's part of 1's digit mode nonzero)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   16 Apr 10 Allow gaussian broadening of DOS
Cu   20 Mar 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,npr,npf,nchan
      double precision eminr,emaxr,eminf,emaxf,efermi,eshft,wt(3),wg
      double precision fitdos(npf,nchan),sigmrq(npf,nchan),
     .  refdos(npr,nchan),outdos(npf,nchan)
C ... Dynamically allocated arrays
      real(8),allocatable:: emeshr(:),emeshf(:),y2(:),fc(:)
C ... Local parameters
      integer i,ich,mode0,stdo,iprint
      double precision der,def,sig,xx,xx2,elo,ehi,ddot
C      double precision seminr,semaxr
      procedure(integer) :: nglob

      if (mode == 0 .or. nchan == 0) return

      mode0 = mod(mode,10)
      stdo = nglob('stdo')

C ... Energy mesh for fit data
      allocate(emeshf(npf))
      def = (emaxf-eminf)/(npf-1)
      do  i = 1, npf
        emeshf(i) = eminf + (i-1)*def
      enddo

C --- Make and interpolate outdos ---
      if (mod(mode,10) /= 0) then
C      seminr = eminr + eshft ! lower bound, ref
C      semaxr = emaxr + eshft ! upper bound, ref
      allocate(emeshr(npr),y2(npr))
      der = (emaxr-eminr)/(npr-1)
      do  i = 1, npr
        emeshr(i) = eminr + (i-1)*der + eshft
      enddo

      do  ich = 1, nchan
C       Get second derivatives for spline
        call cspline(emeshr,refdos(1,ich),npr,1d99,1d99,y2)
C       call prrmsh('dos',emeshr,refdos(1,ich),npr,npr,1)
C       call prrmsh('y2',emeshr,y2,npr,npr,1)
C       Zero out curvature in gap regions
        if (mod(mode0,4) == 2) then
        do  i  = 2, npr-1
          if (refdos(i-1,ich)*refdos(i,ich)*refdos(i+1,ich) == 0)
     .      y2(i) = 0
        enddo
C       call prrmsh('y2 after adjustment',emeshr,y2,npr,npr,1)
        endif
C       Interpolate refdos to fit mesh by cubic spline
        do  i = 1, npf
          xx2 = emeshf(i)
          call csplint(emeshr,refdos(1,ich),y2,npr,xx2,outdos(i,ich),xx)
        enddo
C        call prrmsh('outdos',emeshf,outdos(1,ich),npf,npf,1)
C        call prrmsh('fit dos',emeshf,fitdos(1,ich),npf,npf,1)

C        Convolve
         if (wg /= 0) then
           allocate(fc(npf))
           call dcopy(npf,outdos(1,ich),1,fc,1)
           call convolve(2,0d0,0d0,0d0,0d0,wg,npf,emeshf,fc,
     .       outdos(1,ich))
C          call prrmsh('after convolve',emeshf,outdos(1,ich),npf,npf,1)
           deallocate(fc)
         endif
      enddo
      deallocate(emeshr,y2)
      endif

C --- Generate sigmrq ---
      if (mode/10 /= 0) then
        do  i = 1, npf
        sig = 1
        if (nint(wt(1)) /= 0) then
          xx = (emeshf(i)-wt(2)-efermi)/wt(3)
          if (nint(wt(1)) == 1) then
            sig = dsqrt((exp(xx)+1))
          elseif (nint(wt(1)) == 2) then
            sig = dsqrt((exp(xx**2)+1))
          elseif (nint(wt(1)) == 3) then
            sig = dsqrt(dexp(dabs(xx)))
          elseif (nint(wt(1)) == 4) then
            sig = dsqrt(dexp(dabs(xx**2)/2))
          else
            call rx1('LMASA: fitting weights mode %d not '//
     .        'implemented',wt)
          endif
        endif
        if (i == 1 .or. i == npf) sig = 0
C       Points outside reference window get no weight
        if (emeshf(i) < eminr+eshft) sig = 0
        if (emeshf(i) > emaxr+eshft) sig = 0
C       print *, i, emeshf(i), sig
C       Assign the weights
        do  ich = 1, nchan
          sigmrq(i,ich) = sig
        enddo
      enddo

C ... Get effective range (interval where weights nonzero)
      do  i = 1, npf
        if (sigmrq(i,1) /= 0) then
          elo = emeshf(i)
          exit
        endif
      enddo
      do  i = npf, 1, -1
        if (sigmrq(i,1) /= 0) then
          ehi = emeshf(i)
          exit
        endif
      enddo


      endif

      call info8(30,0,0,'%N Adjust ref DOS, eshft=%d'//
     .  '%?#n#, sculpt##'//
     .  '%?#n#, scale##'//
     .  '%?#n#, wt>0 in (%;4d,%;4d)#%2j#',
     .  eshft,mod(mode0,4)/2,mode0/4,mode/10,elo,ehi,0,0)
      if (mode/10 == 0 .and. mode0 >= 4)
     .  call rx('LMREFDOS requires 10s digit mode to enable scaling')

C --- Scale outdos ---
      if (mod(mode,10) >= 4) then
C      seminr = eminr + eshft ! lower bound, ref
C      semaxr = emaxr + eshft ! upper bound, ref

      call info0(30,0,0,' chan    fit NOS     ref NOS     scale')

      allocate(y2(npf))
      do  i = 1, npf
        if (sigmrq(i,1) /= 0) then
          y2(i) = 1/dsqrt(sigmrq(i,1))
        else
          y2(i) = 0
        endif
      enddo

C      call prrmsh('outdos',emeshf,outdos,npf,npf,1)
C      call prrmsh('wt',emeshf,y2,npf,npf,1)

      do  ich = 1, nchan
        xx = ddot(npf,fitdos(1,ich),1,y2,1) -
     .    (fitdos(1,ich)*y2(1)-fitdos(npf,ich)*y2(npf))/2
        xx2 = ddot(npf,outdos(1,ich),1,y2,1) -
     .    (refdos(1,ich)*y2(1)-refdos(npf,ich)*y2(npf))/2
        if (iprint() >= 30)
     .    write(stdo,333) ich, def*xx, def*xx2,  xx/xx2
  333   format(i4,2f12.5,f10.4)
        call dscal(npf,xx/xx2,outdos(1,ich),1)
      enddo
      deallocate(y2)
      endif

      deallocate(emeshf)

      end
