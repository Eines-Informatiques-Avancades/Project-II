program main_simulation
    ! outline of main program
    use mpi
    use initial_positions_module
    use pbc_mod
    use verlet
    use Forces_and_Energies
    use integrators
    use readers_mod
    !include 'mpif.h'
    implicit none
    real*8, allocatable, dimension(:,:) :: positions, velocities 
    integer :: N,n_steps,n_save_pos,nproc,ierror,iproc
    real*8 :: L,temperature,dt,cutoff,vcutoff,epsilon,sigma,nu
    character (len=500) :: simulation_name

    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
    if (iproc==0) then
        print*, 'START OF SIMULATION.'
        ! reads input file
        call read_parameters("parameters.nml", dt, N, n_steps, n_save_pos, L,&
        simulation_name, temperature, epsilon, cutoff, vcutoff, nu, nproc)

        temperature = temperature/epsilon
        sigma = dsqrt(temperature)
        print*, 'Parameters of the simulation loaded.'
    endif

    call MPI_Bcast(N,           1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(n_steps,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(n_save_pos,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nproc,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    ! REAL64 Broadcasting
    call MPI_Bcast(epsilon,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(cutoff,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nu,          1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(dt,          1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(L,           1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(temperature, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(sigma,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    ! CHARACTER broadcasting
    call MPI_Bcast(simulation_name, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
    ! allocates memory
    allocate(positions(N,3),velocities(N,3))
    print*, 'Needed arrays allocated.'
    ! initialize system
    if (iproc==0) then
        call initial_positions(N, L, positions)
        call initial_velocities(N, temperature, velocities)
        print*, 'System initialized.'
    endif
    call MPI_Bcast(positions,   N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(velocities,  N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)


    ! simulation loop
    call main_loop(MPI_COMM_WORLD,iproc,n_steps,n_save_pos, dt, L, sigma, nu, nproc, cutoff, vcutoff, positions, velocities)
    ! deallocates memory
    if (iproc==0) then
        deallocate(positions,velocities)
        write(*,'(A)') "END OF SIMULATION."
    endif
    call mpi_finalize(ierror)
end program main_simulation