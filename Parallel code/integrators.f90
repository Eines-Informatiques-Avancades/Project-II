module integrators
    use mpi
    use pbc_mod
    use verlet
    use Forces_and_Energies
    use initial_positions_module
    implicit none

    public :: vv_integrator, boxmuller, therm_Andersen
contains
    subroutine vv_integrator(imin,imax,positions, velocities, forces, vlist,nnlist, cutoff, L, dt, max_dist)
        !
        !  Subroutine to update the positions of all particles using the Velocity Verlet
        ! algorithm = performs a single velocity verlet step

        ! Args:
        !    positions  ([3,N]) : positions of all N particles, in reduced units.
        !    velocities ([3,N]) : velocities of all N partciles, in reduced units.
        !    cutoff          () : cutoff value of the interaction.
        !    L               () : length of the sides of the box.
        !    dt              () : value of the integration timestep.
        
        ! Returns:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        real*8, intent(inout) :: max_dist
        integer, intent(in) :: imin, imax, vlist(:), nnlist(:)
        real*8, dimension(:,:), intent(inout) :: positions, velocities, forces
        real*8, intent(in)                                 :: cutoff, L, dt
        integer :: N

        N = size(positions,dim=1)
        ! forces will require the verlet lists
        call VDW_forces(positions, vlist, nnlist, imin, imax, L, cutoff, max_dist, forces)
        !forces(imin:imax,:)=0.d0
        positions(imin:imax,:) = positions(imin:imax,:) + (dt*velocities(imin:imax,:)) + (0.5d0*dt*dt*forces(imin:imax,:))
        call PBC(imin, imax, positions, L,N)
        !velocities = velocities + (0.5d0*dt*forces)

        !call VDW_forces(positions, L, cutoff, forces)
        !velocities = velocities + 0.5d0*dt*forces
    end subroutine vv_integrator

    subroutine main_loop(comm,iproc,N_steps, N_save_pos, dt, L, sigma, nu, nproc, cutoff, vcutoff, positions, velocities)
    implicit none
    integer, intent(in) :: comm, N_steps, N_save_pos,iproc,nproc
    real*8, intent(in) :: dt, cutoff, vcutoff, L, sigma, nu
    real*8, allocatable, dimension(:,:), intent(inout) :: positions, velocities 
    real*8, allocatable :: forces(:,:),max_dist(:)

    real*8 :: KineticEn, PotentialEn, TotalEn, Tinst, press, vcf2
    integer :: N, i, j,unit_dyn=10,unit_ene=11,unit_tem=12,unit_pre=13,imin,imax,subsystems(nproc,2),Nsub,ierror
    integer, allocatable, dimension(:) :: gather_counts, gather_displs, nnlist, vlist
    real*8 :: time
    logical :: update_vlist = .FALSE.
    vcf2 = vcutoff*vcutoff
    N = size(positions, dim=1)
    ! allocation, open files
    ! write initial positions and velocities at time=0
    allocate(forces(N,3))

    if (iproc==0) then
        open(unit_dyn,file = 'dynamics.dat',status="REPLACE")
        open(unit_ene,file = 'energies.dat',status="REPLACE")
        open(unit_tem,file = 'tempinst.dat',status="REPLACE")
        open(unit_pre,file = 'pressure.dat',status="REPLACE")
        do i=1,N
            write(unit_dyn,'(3(f8.3,x))') positions(i,:)
        enddo
        write(unit_dyn,'(A)') " "
    endif
    ! ----------------- Parallel approach -----------------------
    ! Enter nproc as a parameter
    call assign_subsystem(nproc,N,subsystems)
    ! iproc goes from 0 to nproc-1, indices go from 1 to nproc
    imin=subsystems(iproc+1,1) 
    imax=subsystems(iproc+1,2)
    Nsub = imax-imin+1
    ! allocate parallel related structures needed for MPI_ALLGATHERV
    ! and also for the calculation of forces
    allocate(gather_counts(nproc),gather_displs(nproc), nnlist(Nsub), vlist(Nsub*Nsub),max_dist(nproc))
    max_dist = 0.d0
    call MPI_BARRIER(comm,ierror)
    ! Fer Verlet lists
    call verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
    do i=1,N_steps

        time = i*dt
        call MPI_BARRIER(comm,ierror)

        call MPI_ALLGATHER(Nsub, 1, MPI_INTEGER, gather_counts, 1, MPI_INTEGER, comm, ierror)
        ! Calculate displacements for gather operation 
        ! (tells program where to start writing the positions from each worker)
        ! first processor writes first particle, etc
        gather_displs(1) = 0
        do j = 2, nproc
            gather_displs(j) = gather_displs(j - 1) + gather_counts(j - 1)
        end do
        call vv_integrator(imin,imax,positions,velocities,forces,vlist,nnlist,cutoff,L,dt,max_dist(iproc))
        call MPI_BARRIER(comm,ierror)

        ! Perform MPI_ALLGATHERV to gather positions from all processes
        do j=1,3
            call MPI_ALLGATHERV(positions(imin:imax, j), Nsub, MPI_DOUBLE_PRECISION, &
            positions(:,j), gather_counts, gather_displs, MPI_DOUBLE_PRECISION, &
            comm, ierror)
        enddo

        call therm_Andersen(imin,imax,velocities,nu,sigma,Nsub)
        ! -----------------------------------------------------------
        ! compute kinetic, potential and total energies
        if (iproc==0) then
            call kineticE(velocities,KineticEn)
            call potentialE(positions,cutoff,PotentialEn, boxsize=L)
            TotalEn=KineticEn+PotentialEn
            ! compute Instantaneous temperature
            call Tempinst(KineticEn,N,Tinst)
        endif
        ! communicate Tinst to the other workers so they can compute their partial pressure
        call MPI_Bcast(Tinst,     1, MPI_DOUBLE_PRECISION, 0, comm, ierror)
        ! compute pressure
        ! call Pressure(imin,imax,nproc,iproc,positions,L,cutoff,Tinst,press)
        ! write variables to output - positions, energies
        if (iproc==0) then
            print*, ''
            if (MOD(i,N_save_pos).EQ.0) then
                do j=1,N 
                    write(unit_dyn,'(3(e12.3,x))') positions(j,:)
                enddo 
                write(unit_dyn,*) " "
                write(unit_ene,'(4(e12.3,x))') time, KineticEn, PotentialEn, TotalEn
                write(unit_tem,'(2(e12.3,x))') time, Tinst
                write(unit_pre,'(2(e12.3,x))') time, press
            endif
        endif
        ! check if Verlet lists need to be updated
        if (max_dist(iproc) > vcf2) then 
            update_vlist = .TRUE.
            write(*,'(i3,x,A)') iproc,'speaking: Verlet lists need to be updated.'
            write(*,*) max_dist(iproc)
        endif
        call MPI_BARRIER(comm, ierror)
        if (update_vlist) then 
            max_dist(iproc) = 0.d0
            call verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
        endif
    enddo
    if (iproc==0) then
        deallocate(forces)
        close(unit_dyn)
        close(unit_ene)
        close(unit_tem)
        close(unit_pre)
    endif
    end subroutine main_loop

    subroutine boxmuller(sigma, x1, x2, xout1, xout2)
    implicit none
    real*8 :: pi, sigma, x1, x2, xout1, xout2
    pi = 4d0*datan(1d0)
   
    xout1=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*pi*x2)
    xout2=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dsin(2d0*pi*x2)
   
    end subroutine boxmuller

    subroutine therm_Andersen(imin,imax,velocities,nu,sigma,N)
    implicit none
    integer, intent(in) :: imin,imax
    real*8, intent(in) :: nu, sigma
    real*8, dimension(:,:), intent(inout) :: velocities 
    real*8 :: x1, x2, xout1, xout2
    integer :: i, j, seed, N
   
    do i=imin,imax
        call random_number(x1)
        if (x1 < nu) then
            do j=1,3
                call random_number(x1)
                call random_number(x2)
                call boxmuller(sigma, x1, x2, xout1, xout2)
                velocities(i,j)=xout1
            enddo
        endif
    enddo
    end subroutine therm_Andersen
end module integrators

