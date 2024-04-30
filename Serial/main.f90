program main_simulation
    ! outline of main program
    use initial_positions_module
    use pbc_mod
    use Forces_and_Energies
    use integrators
    use readers_mod
    use gr_module
    implicit none
    real*8, allocatable, dimension(:,:) :: positions, velocities 
    integer :: N,n_steps,n_save_pos
    real*8 :: L,temperature,dt,cutoff,epsilon,sigma,nu
    character (len=500) :: simulation_name
    ! reads input file
    print*, 'START OF SIMULATION.'
    call read_parameters("parameters.nml", dt, N, n_steps, n_save_pos, L,&
    simulation_name, temperature, epsilon, cutoff, nu)
    temperature = temperature/epsilon
    sigma = dsqrt(temperature)
    print*, 'Parameters of the simulation loaded.'
    ! allocates memory
    allocate(positions(N,3),velocities(N,3))
    print*, 'Needed arrays allocated.'
    ! initialize system
    call initial_positions(N, L, positions)
    call initial_velocities(N, temperature, velocities)
    print*, 'System initialized.'
    ! simulation loop
    call main_loop(n_steps,n_save_pos, dt, L, sigma, nu, cutoff, positions, velocities)
    ! deallocates memory

    deallocate(positions,velocities)
    write(*,'(A)') "END OF SIMULATION."
end program main_simulation