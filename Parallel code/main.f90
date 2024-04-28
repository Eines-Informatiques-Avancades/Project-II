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
    real*8, allocatable, dimension(:,:) :: velocities 
    integer :: N,n_steps,n_save_pos,nproc,ierror,iproc,i
    real*8 :: L,temperature,dt,cutoff,vcutoff,epsilon,sigma,nu
    character (len=500) :: simulation_name
    real*8, allocatable, dimension(:,:) :: local_positions, positions
    integer :: nsub = 16
    integer :: x(40), gather_counts(4), gather_displs(4)
    integer :: j, numprocs

    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierror)
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
    call MPI_Bcast(vcutoff,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nu,          1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(dt,          1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(L,           1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(temperature, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(sigma,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    ! CHARACTER broadcasting
    call MPI_Bcast(simulation_name, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
    ! allocates memory
    allocate(velocities(N,3))
    print*, 'Needed arrays allocated.'
    ! allocate(local_positions(N/nproc,3))

    ! call assign_subsystem(nproc, N, subsystems) !(i,1)=imin, (i,2)=imax
    N=64
    allocate(local_positions(nsub,3))

    call initial_positions(N, L, local_positions, iproc)

    allocate(positions(N,3))
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    gather_displs = (/0, 16, 32, 48/)
    gather_counts = 16

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    do j=1,3
        call MPI_ALLGATHERV(local_positions(:, j), nsub, MPI_DOUBLE_PRECISION, &
            positions(:,j), gather_counts, gather_displs, &
            MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    if (iproc == 0) then
        open(4, file='initial_positions.dat')
        do i = 1, N
            write(4,'(3(f5.2,x))') positions(j,:)
        end do
        close(4)
    end if

    deallocate(local_positions)

    if (iproc==0) then
        ! call initial_positions(N, L, positions)
        call initial_velocities(N, temperature, velocities)
        print*, 'System initialized.'
    endif
    call MPI_Bcast(velocities,  N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    ! simulation loop
    call main_loop(MPI_COMM_WORLD,iproc,n_steps,n_save_pos, dt, L, sigma, nu, nproc, cutoff, vcutoff, positions, velocities)
    ! deallocates memory
    if (iproc==0) then
        deallocate(velocities)
        deallocate(positions)
        write(*,'(A)') "END OF SIMULATION."
    endif
    call mpi_finalize(ierror)
    ! deallocate(subsystems)
end program main_simulation