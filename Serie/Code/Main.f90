program Main

    ! modules
    use Constants_module
    use Initial_state_module
    use Pbc_module
    use Forces_module
    use Integrator_module
    use pair_distribution_function_module

    implicit none
    ! files
    character(len=80) :: infile, outfile1, outfile2, outfile3, outfile4
    ! parameters
    integer :: nparts
    character(len=3) :: geometry          ! geometry of th system: SC/FCC
    double precision :: in_rho            ! density of the system in g/cm^3
    double precision :: mass              ! molecular mass of the system in g/mol
    double precision :: LJ_sig            ! Lennard-Jones length parameter in A
    double precision :: LJ_eps            ! Lennard-Jones energy parameter in KJ/mol
    double precision :: cutoff            ! radius of cutoff divided by the length of the cell
    double precision :: in_temp           ! room temperature in K
    double precision :: in_temp_init      ! initial temperature of the system in K
    character(len=3) :: bimodal           ! if "Yes" the system initialize velocities using a bimodal distribution associated with 
    !                                       the initial temperature 
    character(len=3) :: disordered_system ! if "Yes" the system realizes 100000 steps at initial temperature
    character(len=3) :: thermostat        ! if "Yes" the system uses a thermostat at room temperature
    character(len=6) :: integrator        ! the integrators can be Euler or Velocity Verlet
    double precision :: in_dt             ! time between two steps in ps
    integer          :: init_steps        ! initial steps in case of disordered system
    integer          :: steps             ! steps of the simulation
    integer          :: measure_steps     ! steps between two measures
    integer          :: traj_steps        ! steps between two trajectories
    integer          :: boxes             ! number of boxes of the histogram
    ! parameters in reduced units
    double precision :: rho, temp, temp_init, dt
    ! observables
    double precision :: time              ! time evolution in reduced units
    double precision :: kin               ! kinetic energy per particle in reduced units
    double precision :: pot               ! potential energy per particle in reduced units
    double precision :: tot               ! total energy per particle in reduced units
    double precision :: tempi             ! temperature of the system in reduced units
    double precision :: pres              ! pressure in reduced units
    double precision :: dr2               ! square displacement in reduced units
    ! conversion factors
    double precision :: lconv, econv, rhoconv, tempconv, pconv, timeconv
    ! simulation box
    integer :: box_m
    double precision :: box_l, box_a
    ! dynamics
    double precision, allocatable :: newpos(:,:), oldpos(:,:), vel(:,:), force(:,:)
    ! pair distribution function
    integer :: measures 
    double precision :: dr
    double precision, allocatable :: gofr(:)
    ! others
    double precision :: dnparts
    double precision :: sigma, nu
    integer :: iter, l, therm_on, integrator_num

    ! read parameters
    infile = "parameters.txt"
    call Read_parameters(infile, nparts, geometry, in_rho, mass, LJ_sig, LJ_eps, cutoff, in_temp, in_temp_init, bimodal, &
    disordered_system, thermostat, integrator, in_dt, init_steps, steps, measure_steps, traj_steps, boxes,  outfile1, outfile2, &
    outfile3, outfile4)

    ! open files
    open(12,file=outfile1)  ! observables evolution
    open(13,file=outfile2)  ! mean square displacement
    open(36,file=outfile3)  ! system evolution

    ! determine the integrator
    if (integrator == "Verlet") then
        integrator_num = 1
        print*, "Verlet integrator"
    else ! Euler
        integrator_num = 0
        print*, "Euler integrator"
    endif

    ! determine vectors dimensions
    allocate(newpos(nparts,3), oldpos(nparts,3), vel(nparts,3), force(nparts,3), gofr(boxes))

    ! calculate conversion fators
    dnparts = dble(nparts)
    lconv = LJ_sig
    econv = LJ_eps/dnparts
    rhoconv = 1d24*mass/(navo*LJ_sig**3)
    tempconv = LJ_eps/(navo*kbol)
    pconv = 1d33*LJ_eps/(navo*LJ_sig**3d0)
    timeconv = 0.1d0*(mass*LJ_sig**2d0/LJ_eps)**0.5d0

    ! calculate parameters in reduced units
    rho = in_rho/rhoconv
    temp = in_temp/tempconv
    temp_init = in_temp_init/tempconv
    dt = in_dt/timeconv


    print*, "The simulation started."
    ! -------------------------------------------------------- prepare the system -------------------------------------------------!
    print*, "The initial configuration is being generated."
    ! initialize the geometry
    if(geometry == "SC ") then
        call sc_lattice(nparts,rho,newpos,box_l,box_a,box_m)
        print*, "SC geometry"
    else
        !call fcc_lattice(nparts,rho,newpos,box_l,box_a,box_m)
    endif

    if (thermostat == "Yes") then 
        therm_on = 1
        print*, "Thermostat on, Temperature (K) = ", in_temp
    else
        therm_on = 0
    endif
    
    ! initialize velocities
    if (bimodal == "Yes") then
        call bimodal_vel(nparts,temp_init,vel)
        print*, "Initial velocities: Bimodal distribution at Temperature = ", in_temp_init
    else
        vel = 0d0
        print*, "Initial velocities equal to zero"
    endif

    ! initialize the histogram (pair distribution function)
    call pair_distribution_function(1,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)

    ! mess up the system
    if (disordered_system == "Yes") then
        print*, "Disordered system at Temperature = ", in_temp_init 
        call LJ_potential(nparts,newpos,cutoff,box_l,1, pot, force)
        nu = 5.d0*dt
        sigma = dsqrt(temp_init)
        do iter  = 1, init_steps
            if (integrator_num == 1) then 
                call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,1)
            else
                call euler_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,1)
            endif
        enddo
    endif
    print*, "The initial configuration was generated, the dynamics starts."
    !------------------------------------------------------------- the system is ready --------------------------------------------!
    ! Initial state
    write(12,*) '# N = ',nparts,', rho (g/cm^3) = ',in_rho,', L (A) = ',box_L*lconv,', a (A) = ',box_a*lconv,', M = ',box_m
    write(12,*) '# time (ps), kin/N (kJ/mol), pot/N (kJ/mol), tot/N (kJ/mol), temp (K), pres (Pa)'
    time = 0.d0
    call LJ_potential(nparts,newpos,cutoff,box_l,1, pot, force)
    call Kinetic_Energy(nparts, vel, kin)
    tot = kin + pot
    tempi = kin*2.d0/(3.d0*dnparts-3.d0)
    call pressure(nparts,rho,tempi,box_l,cutoff,newpos,pres)
    write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
    ! time evolution (dynamics)
    oldpos = newpos
    dr2 = 0.d0
    nu = 5.d0*dt
    sigma = dsqrt(temp)
    do iter  = 1, steps
        ! next point
        time = time + dt
        if (integrator_num == 1) then
            call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,therm_on)
        else
            call euler_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,therm_on)
        endif
        ! mean square displacement
        dr2 = dr2 + sum((newpos - oldpos)**2)
        if (mod(iter,100) == 0) then
            write(13,*) time*timeconv, dr2*lconv**2/dble(nparts)/100.d0
            oldpos = newpos
        end if
        if (mod(iter,measure_steps) == 0) then

            ! sample distances (pair distribution function)
            call pair_distribution_function(2,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)

            ! calculate energy, temperature and pressure
            call Kinetic_Energy(nparts, vel, kin)
            tot = kin + pot
            tempi = kin*2.d0/(3.d0*dnparts-3.d0)
            call pressure(nparts,rho,tempi,box_l,cutoff,newpos,pres)
            write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
        end if

        ! trajectories
        if (mod(iter,traj_steps) == 0) then

            ! Write trajectory in xyz format to be read by VMD
            write(36,*) nparts
            write(36,*)
            do l=1,nparts
                write(36,*) 'A', newpos(l,1)*lconv, newpos(l,2)*lconv, newpos(l,3)*lconv
            end do

        end if
    end do

    ! normalize histogram (pair distribution function)
    call pair_distribution_function(3,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)
    ! write results in real units
    dr = dr*lconv
    call pair_distribution_function(4,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)

    print*, "The dynamics ended, the simulation ended."

    ! close files

    close(12)
    close(13)
    close(36)

end program Main


!**********************************************************************************************************************************!
! read input file
subroutine Read_parameters(params_file, nparts, geometry, density, mass, sigma, epsilon, cutoff, temp_room, temp_init, bimodal, &
    disorder_system, thermostat, integrator, time_step, init_steps, steps, measure_steps, traj_steps, boxes, observables_file, &
    MSD_file, positions_file, gofr_file)

    implicit none
    ! input
    character(len=50), intent(in) :: params_file
    ! output
    integer, intent(out) :: nparts
    character(len=3), intent(out) :: geometry ! SC/FCC
    double precision, intent(out) :: density ! g/cm^3
    double precision, intent(out) :: mass ! g/mol
    double precision, intent(out) :: sigma ! A
    double precision, intent(out) :: epsilon ! KJ/mol
    double precision, intent(out) :: cutoff ! relative to cell length
    double precision, intent(out) :: temp_room ! K
    double precision, intent(out) :: temp_init ! K
    character(len=3), intent(out) :: bimodal ! Si/No
    character(len=3), intent(out) :: disorder_system ! Yes/No
    character(len=3), intent(out) :: thermostat ! Yes/No
    character(len=6), intent(out) :: integrator ! Euler/Verlet
    double precision, intent(out) :: time_step ! ps
    integer, intent(out) :: init_steps ! initial steps in case of disordered system
    integer, intent(out) :: steps ! steps of simulation
    integer, intent(out) :: measure_steps ! steps between two measures
    integer, intent(out) :: traj_steps ! steps between two trajectories
    integer, intent(out) :: boxes ! number of boxes histogram (g(r))
    character(len=80), intent(out) :: observables_file ! energy, temperature, etc evolution
    character(len=80), intent(out) :: MSD_file ! mean square displacement
    character(len=80), intent(out) :: positions_file ! system evolution 
    character(len=80), intent(out) :: gofr_file ! gofr.f90 output file
    ! others
    integer :: my_unit

    open(newunit=my_unit, file=params_file)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*) 
    read(my_unit,*) nparts
    read(my_unit,*)
    read(my_unit,*) geometry
    read(my_unit,*)
    read(my_unit,*) density
    read(my_unit,*) 
    read(my_unit,*) mass
    read(my_unit,*) 
    read(my_unit,*) sigma
    read(my_unit,*) 
    read(my_unit,*) epsilon
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*) 
    read(my_unit,*) cutoff
    read(my_unit,*)
    read(my_unit,*) temp_room
    read(my_unit,*)
    read(my_unit,*) temp_init
    read(my_unit,*)
    read(my_unit,*) bimodal
    read(my_unit,*) 
    read(my_unit,*) disorder_system
    read(my_unit,*) 
    read(my_unit,*) thermostat
    read(my_unit,*) 
    read(my_unit,*) integrator
    read(my_unit,*) 
    read(my_unit,*) time_step
    read(my_unit,*)
    read(my_unit,*) init_steps
    read(my_unit,*) 
    read(my_unit,*) steps
    read(my_unit,*) 
    read(my_unit,*) measure_steps
    read(my_unit,*)
    read(my_unit,*) traj_steps
    read(my_unit,*)
    read(my_unit,*) boxes
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*) 
    read(my_unit,*) observables_file
    read(my_unit,*) 
    read(my_unit,*) MSD_file
    read(my_unit,*) 
    read(my_unit,*) positions_file
    read(my_unit,*)
    read(my_unit,*) gofr_file
    close(my_unit)

    return

end subroutine Read_parameters

!**********************************************************************************************************************************!
