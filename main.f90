program main_simulation
    ! outline of main program
    use initial_positions_module
    use pbc_mod
    use Forces_and_Energies
    use integrators
    use readers_mod
    implicit none
    real*8, allocatable, dimension(:,:) :: positions, velocities 
    integer :: N,n_steps
    real*8 :: L,temperature,dt,cutoff
    character (len=500) :: simulation_name
    ! reads input file
    call read_parameters("parameters.nml", dt, N, n_steps, L,&
    simulation_name, temperature, cutoff)
    ! allocates memory
    allocate(positions(N,3),velocities(N,3))
    ! initialize system
    call initial_positions(N, L, positions)
    call initial_velocities(N,temperature,velocities)
    ! simulation loop
    call main_loop(n_steps, dt, L, cutoff, positions, velocities)
    ! deallocates memory
    deallocate(positions,velocities)
    write(*,'(A)') "END OF SIMULATION."
end program main_simulation