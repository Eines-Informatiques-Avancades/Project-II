module integrators
    use mpi
    use pbc_mod
    use verlet
    use Forces_and_Energies
    use initial_positions_module
    implicit none

    public :: vv_integrator1,vv_integrator2, boxmuller, therm_Andersen
contains
    subroutine vv_integrator1(imin,imax,positions, velocities, forces, vlist,nnlist, cutoff, L, dt)

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
        !   max_dist            :
        
        ! Returns:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none

        real*8 :: max_dist
        
        integer, intent(in) :: imin, imax, vlist(:), nnlist(:)
        real*8, dimension(:,:), intent(inout) :: positions, velocities, forces
        real*8, intent(in)                                 :: cutoff, L, dt
        integer :: N

        N = size(positions,dim=1)
        ! forces will require the verlet lists
        call VDW_forces(positions, vlist, nnlist, imin, imax, L, cutoff, max_dist, forces)
        positions(imin:imax,:) = positions(imin:imax,:) + (dt*velocities(imin:imax,:)) + (0.5d0*dt*dt*forces(imin:imax,:))
        call PBC(imin, imax, positions, L,N)
        velocities(imin:imax,:) = velocities(imin:imax,:) + (0.5d0*dt*forces(imin:imax,:))
    end subroutine vv_integrator1

    subroutine vv_integrator2(imin,imax,positions, velocities, forces, vlist,nnlist, cutoff, L, dt, max_dist)
        !
        !  Subroutine to update the positions of each worker's assigned particles using the second step 
        !  of the Velocity Verlet integrator.

        implicit none
        real*8, intent(out) :: max_dist
        integer, intent(in) :: imin, imax, vlist(:), nnlist(:)
        real*8, dimension(:,:), intent(inout) :: positions, velocities, forces
        real*8, intent(in)                                 :: cutoff, L, dt


        call VDW_forces(positions, vlist, nnlist, imin, imax, L, cutoff, max_dist, forces)
        velocities(imin:imax,:) = velocities(imin:imax,:) + 0.5d0*dt*forces(imin:imax,:)
    end subroutine vv_integrator2

    subroutine main_loop(comm,iproc,N_steps, N_save_pos, dt, L, sigma, nu, nproc, cutoff, vcutoff, positions, velocities)
    implicit none
    integer, intent(in) :: comm, N_steps, N_save_pos,iproc,nproc
    real*8, intent(in) :: dt, cutoff, vcutoff, L, sigma, nu
    real*8, allocatable, dimension(:,:), intent(inout) :: positions, velocities 
    real*8, allocatable :: forces(:,:)

    real*8 :: KineticEn, PotentialEn, TotalEn, Tinst, press, vcf2
    integer :: N, i, j,unit_dyn=10,unit_ene=11,unit_tem=12,unit_pre=13,imin,imax,subsystems(nproc,2),Nsub,ierror
    integer, allocatable, dimension(:) :: gather_counts, gather_displs, nnlist, vlist

    real*8 :: time, max_dist,local_kineticEn
    logical :: update_vlist = .FALSE.
    !----------------- NEW-------------------------------------
    real*8 :: temp,volume, Virialterm, global_Virialterm, ierr
    ! Temp i volume esta definit ja? Per si de cas ho posso, com a prova!
    temp=298.0
    volume=L**3.0 
    
    vcf2 = vcutoff*vcutoff
    kineticEn=0.d0

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

    ! Enter nproc as a parameter
    call assign_subsystem(nproc,N,subsystems)
    ! iproc goes from 0 to nproc-1, indices go from 1 to nproc
    imin=subsystems(iproc+1,1) 
    imax=subsystems(iproc+1,2)
    Nsub = imax-imin+1
    ! allocate parallel related structures needed for MPI_ALLGATHERV
    ! and also for the calculation of forces
    allocate(gather_counts(nproc),gather_displs(nproc), nnlist(Nsub), vlist(Nsub*Nsub))

    call MPI_BARRIER(comm,ierror)
    ! Build Verlet lists
    call verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)

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
        call vv_integrator1(imin,imax,positions,velocities,forces,vlist,nnlist,cutoff,L,dt)
        ! Perform MPI_ALLGATHERV to update positions from all processes
        
        call MPI_BARRIER(comm,ierror)
        do j=1,3
            call MPI_ALLGATHERV(positions(imin:imax, j), Nsub, MPI_DOUBLE_PRECISION, &
            positions(:,j), gather_counts, gather_displs, MPI_DOUBLE_PRECISION, &
            comm, ierror)
        enddo

        ! Second step of Verlet integration, recalculates actual forces
        call vv_integrator2(imin,imax,positions,velocities,forces,vlist,nnlist,cutoff,L,dt,max_dist)
        call therm_Andersen(imin,imax,velocities,nu,sigma,Nsub)
        ! -----------------------------------------------------------
        ! compute kinetic, potential and total energies
        call kineticE(imin, imax,velocities,local_kineticEn)
        call MPI_BARRIER(comm,ierror)
        ! gets each local kinetic energy, sums it all into kineticEn, which is a variable that processor 0 keeps
        call MPI_REDUCE(local_kineticEn, kineticEn, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, &
                        ierror)
        if (iproc==0) then
            call potentialE(positions,cutoff,PotentialEn, boxsize=L)
            TotalEn=KineticEn+PotentialEn
            ! compute Instantaneous temperature
            call Tempinst(KineticEn,N,Tinst)
        endif
        ! communicate Tinst to the other workers so they can compute their partial pressure
        call MPI_Bcast(Tinst,     1, MPI_DOUBLE_PRECISION, 0, comm, ierror)
        
        !!!!! --------------------- NEW ---------------------------------------------
        !!  -------------------------------------------------------------------------
	! compute pressure
        
        call MPI_ALLGATHER(Nsub, 1, MPI_INTEGER, gather_counts, 1, MPI_INTEGER, comm, ierror)
        ! Calculate displacements for gather operation
        ! (tells program where to start writing the positions from each worker)
        ! first processor writes first particle, etc
        gather_displs(1) = 0
        do j = 2, nproc
            gather_displs(j) = gather_displs(j - 1) + gather_counts(j - 1)
        end do
        call MPI_BARRIER(comm,ierror)
        call Pressure(vlist,nnlist,imin,imax,positions,L,cutoff,temp,max_dist,Virialterm)
        print*, Virialterm
        
        call MPI_Reduce(Virialterm, global_Virialterm, 1,& 
            MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierror)
      
        
        ! Pressure units are J/m³
          ! we multiply ideal gas term *Kb=1.38*10**(-23)
         press= (dble(N)*temp)/volume + (1.d0/(3.d0*volume))*Virialterm
          !Pressure in reduced units
           
        print*, press
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
        if (max_dist > vcf2) then 
            update_vlist = .TRUE.
            write(*,'(i3,x,A)') iproc,'speaking: Verlet lists need to be updated.'
            write(*,*) max_dist
        endif
        call MPI_BARRIER(comm, ierror)
        if (update_vlist) then 
            max_dist = 0.d0
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
    deallocate(nnlist, vlist)
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

