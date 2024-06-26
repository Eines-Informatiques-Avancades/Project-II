program main_simulation
    ! outline of main program
    use mpi
    use initial_positions_module
    use pbc_mod
    use verlet
    use Forces_and_Energies
    use integrators
    use readers_mod
    use block_average_module
    use gr_module
    !include 'mpif.h'
    !!! --- Author: Paula Sierra and Emma Valdés --- !!!
    !!! --- Contributors: Quim Badosa and Guillem Arasa --- !!!
    implicit none
    real*8, allocatable, dimension(:,:) :: positions,local_positions, all_positions, velocities 
    integer, allocatable, dimension(:) :: gather_counts, gather_displs
    real*8 :: L,temperature,dt,cutoff,vcutoff,epsilon,sigma,nu,wtime
    integer :: N,n_steps,n_save_pos,nproc,ierror,iproc,i,j,nsub
    character (len=500) :: simulation_name

    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
    call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)

    if (iproc==0) then
        wtime = mpi_wtime()
    endif

    if (iproc==0) then
        print*, 'START OF SIMULATION.'
        ! reads input file
        call read_parameters("parameters.nml", dt, N, n_steps, n_save_pos, L,&
        simulation_name, temperature, epsilon, cutoff, vcutoff, nu)

        temperature = temperature/epsilon
        sigma = dsqrt(temperature)
        print*, 'Parameters of the simulation loaded.'
    endif

    call MPI_Bcast(N,           1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(n_steps,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(n_save_pos,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

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
    Nsub = N/nproc
    allocate(local_positions(Nsub,3))
    allocate(positions(N,3))
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    if (N /= nproc**3) then
        print*, 'Initial positions in serial.'
        if (iproc == 0) then
           call initial_positions_serial(N, L, positions)
        endif
        call MPI_Bcast(positions, N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    else
        print*, 'Initial positions in parallel.'
        call initial_positions_parallel(N, L, local_positions, iproc)
        allocate(gather_counts(nproc), gather_displs(nproc))
        gather_counts = N/nproc
        gather_displs(1) = 0
        do i=2,nproc 
            gather_displs(i) = gather_displs(i-1)+gather_counts(i-1)
        enddo
    
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        do j=1,3
        call MPI_ALLGATHERV(local_positions(:, j), gather_counts(iproc+1), MPI_DOUBLE_PRECISION, &
                            positions(:,j), gather_counts, gather_displs, &
                            MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        deallocate(gather_counts,gather_displs)
    endif

    if (iproc==0) then
        call initial_velocities(N, temperature, velocities)
        print*, 'System initialized.'
    endif
    call MPI_Bcast(velocities,  N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    ! simulation loop
    call main_loop(MPI_COMM_WORLD,ierror,iproc,n_steps,n_save_pos, dt, L, sigma, nu, nproc, cutoff, vcutoff, positions, velocities)

    if (iproc==0) then
        wtime = mpi_wtime() - wtime
        write(*,'(A,F10.2,A)') "Elapsed time: ", wtime, " seconds."
    endif

    ! deallocates memory
    deallocate(velocities,local_positions,positions)

    if (iproc==0) then
        ! Note: Adjust block_size
        call compute_and_save_block_averages('tempinst.dat', 1000, 2, 'temp_blockavg.dat')
        call compute_and_save_block_averages('pressure.dat', 1000, 2, 'pressure_blockavg.dat')
        call compute_and_save_block_averages('g_r.dat', 1000, 2, 'gr_blockavg.dat')
        call compute_and_save_block_averages('energies.dat', 1000, 2, 'energykin_blockavg.dat')
        call compute_and_save_block_averages('energies.dat', 1000, 3, 'energypot_blockavg.dat')
        call compute_and_save_block_averages('energies.dat', 1000, 4, 'energytot_blockavg.dat')
    endif

    write(*,'(A)') "END OF SIMULATION."

    call mpi_finalize(ierror)
end program main_simulation
