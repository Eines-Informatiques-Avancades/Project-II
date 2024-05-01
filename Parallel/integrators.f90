module integrators
    use mpi
    use pbc_mod
    use verlet
    use Forces_and_Energies
    use initial_positions_module
    use gr_module
    implicit none

    public :: vv_integrator1,vv_integrator2, boxmuller, therm_Andersen
contains
    subroutine vv_integrator1(imin,imax,positions, velocities, vlist,nnlist, cutoff, L, dt)
        !!! --- Authors: Emma Valdés and Paula Sierra --- !!!
        !
        !  Subroutine to update the positions of each worker's assigned particles using the first step 
        !  of the Velocity Verlet integrator.

        ! Args:
        !   imin                : Assigned particle with the lowest index.
        !   imax                : Assigned particle with the highest index.
        !   positions  ([3,N])  : positions of all N particles, in reduced units.
        !   velocities ([3,N])  : velocities of all N particles, in reduced units.
        !   forces ([3,N])      : forces perceived by each particle, in reduced units.
        !   vlist ([Nsub**2])   : Worker's Verlet list. Nsub is the number of particles assigned to the worker.
        !   nnlist ([Nsub])     : List with the number of neighbors of each particle of the worker's assigned subsystem.           
        !   cutoff              : cutoff value of the VdW interaction.
        !   L                   : length of the sides of the box.
        !   dt                  : value of the integration timestep.
        
        ! Returns:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        
        integer, intent(in) :: imin, imax, vlist(:), nnlist(:)
        real*8, dimension(:,:), intent(inout) :: positions, velocities
        real*8, intent(in)                                 :: cutoff, L, dt
        real*8, dimension(size(positions,dim=1),3) :: forces
        integer :: N

        N = size(positions,dim=1)
        ! forces will require the verlet lists
        call VDW_forces(positions, vlist, nnlist, imin, imax, cutoff, forces, L)
        positions(imin:imax,:) = positions(imin:imax,:) + (dt*velocities(imin:imax,:)) + (0.5d0*dt*dt*forces(imin:imax,:))
        call PBC(imin, imax, positions, L,N)
        velocities(imin:imax,:) = velocities(imin:imax,:) + (0.5d0*dt*forces(imin:imax,:))
    end subroutine vv_integrator1

    subroutine vv_integrator2(imin,imax,positions, velocities, vlist,nnlist, cutoff, L, dt)
        !!! --- Authors: Emma Valdés and Paula Sierra --- !!!
        !  Subroutine to update the positions of each worker's assigned particles using the second step 
        !  of the Velocity Verlet integrator.

        implicit none
        integer, intent(in) :: imin, imax, vlist(:), nnlist(:)
        real*8, dimension(:,:), intent(inout) :: positions, velocities
        real*8, intent(in)                                 :: cutoff, L, dt
        real*8, dimension(size(positions, dim=1),3) :: forces

        call VDW_forces(positions, vlist, nnlist, imin, imax, cutoff, forces, L)
        velocities(imin:imax,:) = velocities(imin:imax,:) + 0.5d0*dt*forces(imin:imax,:)
    end subroutine vv_integrator2

    subroutine main_loop(comm,ierror,iproc,N_steps, N_save_pos, dt, L, sigma, nu, nproc, cutoff, vcutoff, positions, velocities)
        !!! --- Authors: Emma Valdés and Paula Sierra --- !!!
        !!! --- Contributors: Guillem Arasa and Quim Badosa --- !!!
    implicit none
    integer, intent(inout) :: ierror
    integer, intent(in) ::  comm,N_steps, N_save_pos,iproc,nproc
    real*8, intent(in) :: dt, cutoff, vcutoff, L, sigma, nu
    real*8, allocatable, dimension(:,:), intent(inout) :: positions, velocities 
    real*8, allocatable :: forces(:,:)
    real*8 :: KineticEn, PotentialEn, TotalEn, Tinst, press, vcf2
    double precision, dimension(:), allocatable :: g_r

    integer :: unit_dyn=10,unit_ene=11,unit_tem=12,unit_pre=13,unit_g_r=14
    integer :: imin,imax,subsystems(nproc,2),Nsub,N, i, j, n_update=10
    integer, allocatable, dimension(:) :: gather_counts, gather_displs, nnlist, vlist
    real*8 :: time, max_dist,local_kineticEn,local_potentialEn, global_max_dist
    real*8 :: volume, Virialterm, global_Virialterm
    
    volume=L**3.0 
    
    vcf2 = vcutoff*vcutoff
    kineticEn=0.d0
    potentialEn=0.d0
    N = size(positions, dim=1)
    ! allocation, open files
    ! write initial positions and velocities at time=0
    !allocate(forces(N,3))

    if (iproc==0) then
        open(unit_dyn,file = 'dynamics.XYZ',status="REPLACE")
        open(unit_ene,file = 'energies.dat',status="REPLACE")
        open(unit_tem,file = 'tempinst.dat',status="REPLACE")
        open(unit_pre,file = 'pressure.dat',status="REPLACE")
        open(unit_g_r,file='g_r.dat',status='REPLACE')

        ! open(unit_g_r,file='g_r.dat',status='REPLACE')
        write(unit_dyn,*) N
        write(unit_dyn,*) 0,"TIMESTEP:", 0
        do i=1,N
            write(unit_dyn,*) "Kr",positions(i,:)
        enddo
    endif

    ! Divide system in subsystems
    call assign_subsystem(nproc,N,subsystems)
    imin = subsystems(iproc+1,1)
    imax = subsystems(iproc+1,2)
    Nsub = imax-imin+1
    allocate(gather_counts(nproc),gather_displs(nproc), nnlist(Nsub), vlist(N*Nsub),g_r(1))

    ! Quantities needed to share information between processes
    gather_counts = Nsub
    call MPI_BARRIER(comm,ierror)
    
    ! Build Verlet lists
    call verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
    
    ! allocate(g_r(1))
    do i=1,N_steps
        time = i*dt
        call MPI_ALLGATHER(Nsub, 1, MPI_INTEGER, gather_counts, 1, MPI_INTEGER, comm, ierror)
        ! Calculate displacements for gather operation 
        ! (tells program where to start writing the positions from each worker)
        ! first processor writes first particle, etc
        gather_displs(1) = 0
        do j = 2, nproc
            gather_displs(j) = gather_displs(j - 1) + gather_counts(j - 1)
        end do
        call MPI_BARRIER(comm,ierror)

        ! First step of Verlet integration
        call vv_integrator1(imin,imax,positions,velocities,vlist,nnlist,cutoff,L,dt)
        ! Perform MPI_ALLGATHERV to update positions from all processes
        call MPI_BARRIER(comm,ierror)
        do j=1,3
            call MPI_ALLGATHERV(positions(imin:imax, j), Nsub, MPI_DOUBLE_PRECISION, &
            positions(:,j), gather_counts, gather_displs, MPI_DOUBLE_PRECISION, &
            comm, ierror)
        enddo
        
        ! Second step of Verlet integration, recalculates actual forces
        call vv_integrator2(imin,imax,positions,velocities,vlist,nnlist,cutoff,L,dt)
        ! Andersen thermostat modifies some velocities
        call therm_Andersen(imin,imax,velocities,nu,sigma,Nsub)
      
        ! compute kinetic, potential and total energies
        call kineticE(imin, imax, velocities,local_kineticEn)
        call potentialE(positions, vlist, nnlist, imin, imax,cutoff,local_PotentialEn, boxsize=L)
        call MPI_BARRIER(comm,ierror)
        ! gets each local kinetic energy, sums it all into kineticEn, which is a variable that processor 0 keeps
        call MPI_REDUCE(local_kineticEn, kineticEn, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, &
                        ierror)
        call MPI_REDUCE(local_potentialEn, potentialEn, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, &
        ierror)
        
        if (iproc==0) then
             call calculate_g_r(positions, L, N, 0.1d0, g_r, 2.d0, 200)
        endif


        ! write variables to output - positions, energies
        if (iproc==0) then
        ! Pressure units are J/m³
          ! we multiply ideal gas term *Kb=1.38*10**(-23)
            press= (dble(N)*Tinst)/volume + (1.d0/(3.d0*volume))*global_Virialterm
            !Pressure in reduced units

            if (MOD(i,N_save_pos).EQ.0) then
                write(unit_dyn,*) N
                write(unit_dyn,*) time,"TIMESTEP:", i/100
                do j=1,N 
                    write(unit_dyn,*) "Kr",positions(j,:)
                enddo 
                write(unit_ene,'(4(e12.3,x))') time, KineticEn, PotentialEn, TotalEn
                write(unit_tem,'(2(e12.3,x))') time, Tinst
                write(unit_pre,'(2(e12.3,x))') time, press
                write(unit_g_r,'(2(e12.3,x))') time, g_r(1)
            endif
            ! check if Verlet lists need to be updated
            if (MOD(i, n_update).EQ.0) then
                call MPI_BARRIER(comm, ierror)
                call MPI_Reduce(max_dist,global_max_dist,1,&
                    MPI_DOUBLE_PRECISION,MPI_MAX, 0, comm, ierror)
                call MPI_Bcast(global_max_dist, 1, MPI_DOUBLE_PRECISION, 0, comm, ierror)
    
                if (global_max_dist  > vcf2) then
                    call verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
                endif
            endif
        endif
    enddo
    
    if (iproc==0) then
        close(unit_dyn)
        close(unit_ene)
        close(unit_tem)
        close(unit_pre)
        close(unit_g_r)
    endif

    deallocate(nnlist, vlist, gather_counts, gather_displs,g_r)
    end subroutine main_loop

    subroutine boxmuller(sigma, x1, x2, xout1, xout2)
        !!! --- Authors: Paula Sierra --- !!!
        !!! --- Contributor: Rocío Aragoneses --- !!!

        implicit none
        real*8 :: pi, sigma, x1, x2, xout1, xout2
        pi = 4d0*datan(1d0)

        xout1=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*pi*x2)
        xout2=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dsin(2d0*pi*x2)
    end subroutine boxmuller

    subroutine therm_Andersen(imin,imax,velocities,nu,sigma,N)
        !!! --- Author: Paula Sierra --- !!!
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

