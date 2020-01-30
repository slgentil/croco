#ifdef GO_MEDDY

! to be included in nr_inittracer_meddy.F90 7 fev. 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        icase = 2                      !!!! icase > 100 pour Nvar

        zmid = -h(1,1)/2           ! mid-depth
!        Nfreq2 = 5.0d-6
!        Nfreq = sqrt(Nfreq2)
        Nfreq = 2.23d-3
!       Nfreq = 2.d-3
        Nfreq2 = Nfreq**2
!       pi = 4.*atan(1.)
!       print *,'pi=', pi

!       x1=0
!       y1=0
        x1=xl/2.
        y1=el/2.
        linear_balance   = .false.      ! linear balance       
        random_init      = .true.      ! put to true to get random_seeds
        random_pert_urot = .false.
        Tcoef            = (4.0d0/25.0d0)
        Scoef            = (4.0d0/5.0d0)

!..........................................................................

      if (icase.eq.1) then   ! 'lorentzienne' shape, Ncst (papier QG: power = 1)
        stretch = 0.9d0
        Burger  = 0.15d0
        thick_y = 3.2e4
        power   = 1
        z1      = -1100

! ..... Background stratification
        if (mynode == 0)  write(STDOUT,*) 'N CONSTANT'
        dhh    = 4000
        shift  = 0.0
        deltaN = 0.0
        zth1   = 1.0
        zth2   = 1.0
        delta1 = 0.0
        delta2 = 0.0
        wid1   = 1.0
        wid2   = 1.0

! ..... T et S
        T0 = 12.0d0
        S0 = 36.6d0
        FRACTS = 0.5d0

      elseif (icase.eq.2) then              !  'gaussienne' shape, Ncst
!       stretch = 0.167d0
!       stretch = 0.3d0
!       stretch = 0.4d0
!       stretch = 0.98d0
!       stretch = 1.d0
!       stretch = 1.33d0
        stretch = 2.5d0
!       stretch = 3.d0
!       stretch = 6.d0

        Burger  = 0.3d0
!       Burger  = 0.2d0
!       Burger  = 0.15d0
        thick_y = 2.8e4
        power   = 1
        z1      = -1100
!       z1      = -1200

! ..... Background stratification
        if (mynode == 0)  write(STDOUT,*) 'N CONSTANT'
        dhh    = 4000
        shift  = 0.0
        deltaN = 0.0
        zth1   = 1.0
        zth2   = 1.0
        delta1 = 0.0
        delta2 = 0.0
        wid1   = 1.0
        wid2   = 1.0

! ..... T et S
        T0 = 12.0d0
        S0 = 36.6d0
        FRACTS = 0.5d0


!..........................................................................
!..........................................................................


      elseif (icase.eq.50) then              !  GO data
        thick_y = 2.8e4
        z1      = -1100
        stretch = 0.0d0
        Burger  = 0.0d0

! ..... Background stratification
        if (mynode == 0)  write(STDOUT,*) 'N GO'
!        dhh    = 4000
!        shift  = 0.0
!        deltaN = 0.0
!        zth1   = 1.0
!        zth2   = 1.0
!        delta1 = 0.0
!        delta2 = 0.0
!        wid1   = 1.0
!        wid2   = 1.0

! ..... T et S
        T0 = 9.6d0
        S0 = 36.6d0
        FRACTS = 0.45d0




!..........................................................................
!..........................................................................

	elseif (icase.eq.101) then                  !  'lorentzienne' shape, Nvar
        stretch = 0.88d0
        Burger  = 0.15d0
        thick_y = 3.5d4
        power   = 2
        z1      = -1070

! ..... Background stratification
        if (mynode == 0) 
     &     write(STDOUT,*) 'thermoclines set up at',zth1, 'and', zth2,' (m)'
        dhh     = 400
        shift  = 1./6.
        deltaN = 0.24
        zth1   = 624
        zth2   = 1332
        delta1 = 2.1d-3
        delta2 = 692d-6
        wid1   = 0.8/120
        wid2   = wid1*1.3


! ..... T et S
        T0 = 12.0d0
        S0 = 36.6d0
        FRACTS = 0.5d0

	elseif (icase.eq.102) then           !  'lorentzienne' shape, Nvar Lien
!        stretch = 0.88d0
!        Burger = 0.15d0
!        thick_y = 3.5d4
        stretch = 0.9d0
        Burger = 0.15d0
        thick_y = 3.5d4
        power  = 2
        z1     = -1000

! ..... Background stratification
        if (mynode == 0) 
     &   write(STDOUT,*) 'thermoclines set up at',zth1, 'and', zth2,' (m)'
        dhh    = 4000
        shift  = 0.0
        deltaN = 0.0

!        zunit  = 120.		!meters
        zth1   = 550
!        zth2   = 120
        zth2   = 1360
!        delta1 = 1.4*2./(0.9/120)*Nfreq2
!        delta2 = 2./7.*2./(1.1/120)*Nfreq2
        delta1 = 0.5*2./(0.9/120)*Nfreq2
        delta2 = 0.2*2./(1.1/120)*Nfreq2
        wid1   = 0.9/120
        wid2   = 1.1/120

! ..... T et S
        T0 = 12.0d0
        S0 = 36.6d0
        FRACTS = 0.5d0


!..........................................................................

      else
        if (mynode == 0) then
           print *,'================================'
           print *,'icase not found: job is stopped '
           print *,'================================'
        endif
        stop
      endif
        

        Rossby = -2*stretch*Burger
        thick_z = thick_y*(f(1,1)/Nfreq)*sqrt(Burger)

! ..... random perturbation
      if (random_init) then     
        nmdef  = 1
        energy = 0.1                    ! energy used to scale amplitude of perturbation
        facsolpert = 1.                 ! =0 for perturb only, =1  superposition (random pert)
      endif

				!.....................

#endif  /*  GO_MEDDY */
