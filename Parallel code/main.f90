program main_simulation
    ! outline of main program
    use initial_positions_module
    use pbc_mod
    use Forces_and_Energies
    use integrators
    use readers_mod
    implicit none
    include 'mpif.h'
    real*8, allocatable, dimension(:,:) :: positions, velocities 
    integer :: N,n_steps,n_save_pos,nproc,ierror,iproc
    real*8 :: L,temperature,dt,cutoff,epsilon,sigma,nu
    character (len=500) :: simulation_name

    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
    if (iproc==0) then
        print*, 'START OF SIMULATION.'
        ! reads input file
        call read_parameters("parameters.nml", dt, N, n_steps, n_save_pos, L,&
        simulation_name, temperature, epsilon, cutoff, nu, nproc)
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
    endif
    ! simulation loop
    call main_loop(iproc,n_steps,n_save_pos, dt, L, sigma, nu, nproc, cutoff, positions, velocities)
    ! deallocates memory
    if (iproc==0) then
        deallocate(positions,velocities)
        write(*,'(A)') "END OF SIMULATION."
    endif
    call mpi_finalize(ierror)
end program main_simulation