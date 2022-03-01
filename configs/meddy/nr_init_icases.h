#ifdef GO_MEDDY

! to be included in nr_inittracer_meddy.F90 7 fev. 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zmid = -h(1,1)/2           ! mid-depth
        Nfreq = 2.23d-3
        Nfreq2 = Nfreq**2

        x1=xl/2.
        y1=el/2.
        linear_balance   = .false.      ! linear balance       
        random_init      = .true.      ! put to true to get random_seeds

!..........................................................................

        stretch = 0.98d0

        Burger  = 0.15d0
        thick_y = 2.8e4
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
        FRACTS = 0.5d0
        

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
