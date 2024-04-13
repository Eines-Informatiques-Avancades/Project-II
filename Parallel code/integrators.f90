module integrators
    use pbc_mod
    use Forces_and_Energies
    use initial_positions_module

    implicit none
    public :: vv_integrator, boxmuller, therm_Andersen

contains
    subroutine vv_integrator(positions, velocities, forces, cutoff, L, dt)
        !
        !  Subroutine to update the positions of all particles using the Velocity Verlet
        ! algorithm = performs a single velocity verlet step

        ! Args:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.
        !    cutoff          (REAL64) : cutoff value of the interaction.
        !    L               (REAL64) : length of the sides of the box.
        !    dt              (REAL64) : value of the integration timestep.
        
        ! Returns:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        real*8, allocatable, dimension(:,:), intent(inout) :: positions, velocities, forces
        real*8, intent(in)                                 :: cutoff, L, dt
        integer :: N
        N = size(positions,dim=1)
        call VDW_forces(positions, L, cutoff, forces)
        positions = positions + (dt*velocities) + (0.5d0*dt*dt*forces)
        call PBC(positions, L,N)
        velocities = velocities + (0.5d0*dt*forces)

        call VDW_forces(positions, L, cutoff, forces)
        velocities = velocities + 0.5d0*dt*forces
    end subroutine vv_integrator

    subroutine main_loop(N_steps, N_save_pos, dt, L, sigma, nu, cutoff, positions, velocities)
    implicit none
    integer, intent(in) :: N_steps, N_save_pos
    real*8, intent(in) :: dt, cutoff,L, sigma, nu
    real*8, allocatable, dimension(:,:), intent(inout) :: positions, velocities 
    real*8, allocatable, dimension(:,:) :: forces
    real*8 :: KineticEn, PotentialEn, TotalEn, Tinst, press
    integer :: N, i, j,unit_dyn=10,unit_ene=11,unit_tem=12,unit_pre=13
    real*8 :: time
    N = size(positions, dim=1)
    allocate(forces(N,3))
    ! open files
    open(unit_dyn,file = 'dynamics.dat',status="REPLACE")
    open(unit_ene,file = 'energies.dat',status="REPLACE")
    open(unit_tem,file = 'tempinst.dat',status="REPLACE")
    open(unit_pre,file = 'pressure.dat',status="REPLACE")
    do i=1,N
        write(unit_dyn,'(3(f8.3,x))') positions(i,:)
    enddo
    write(unit_dyn,'(A)') " "
    ! write initial positions and velocities at time=0
    do i=1,N_steps 
        time = i*dt
        ! Enter nproc as a parameter
        call assign_subsystem(nproc,N,subsystems)
        ! ----------------- Parallel approach -----------------------
        call mpi_init(ierror)
        call mpi_comm_rank(MPI_COMM_WORD,iproc,ierror)

        imin=subsystems(iproc,1)
        imax=subsystems(iproc,2)

        ! Fer Verlet lists
        
        call vv_integrator(positions,velocities,forces,cutoff,L,dt)
        call therm_Andersen(velocities,nu,sigma,N)

        call mpi_finalize(ierror)
        ! -----------------------------------------------------------
        ! compute kinetic, potential and total energies
        call kineticE(velocities,KineticEn)
        call potentialE(positions,cutoff,PotentialEn, boxsize=L)
        TotalEn=KineticEn+PotentialEn

        ! compute Instantaneous temperature
        call Tempinst(KineticEn,N,Tinst)
        ! compute pressure
        call Pressure (positions,L,cutoff,Tinst,press)
        ! write variables to output - positions, energies

        if (MOD(i,N_save_pos).EQ.0) then
            do j=1,N 
                write(unit_dyn,'(3(e12.3,x))') positions(j,:)
            enddo 
            write(unit_dyn,'(A)') " "
            write(unit_ene,'(4(e12.3,x))') time, KineticEn, PotentialEn, TotalEn
            write(unit_tem,'(2(e12.3,x))') time, Tinst
            write(unit_pre,'(2(e12.3,x))') time, press
        endif
    enddo
    deallocate(forces)
    close(unit_dyn)
    close(unit_ene)
    close(unit_tem)
    close(unit_pre)
    end subroutine main_loop

    subroutine boxmuller(sigma, x1, x2, xout1, xout2)
    implicit none
    real*8 :: pi, sigma, x1, x2, xout1, xout2
    pi = 4d0*datan(1d0)
   
    xout1=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*pi*x2)
    xout2=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dsin(2d0*pi*x2)
   
    end subroutine boxmuller

    subroutine therm_Andersen(velocities,nu,sigma,N)
    implicit none
    integer :: i, j, seed, N
    real*8, intent(in) :: nu, sigma
    real*8, allocatable, dimension(:,:), intent(inout) :: velocities 
    real*8 :: x1, x2, xout1, xout2
   
    do i=1,N
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

    subroutine verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
        implicit none
        integer, intent(in) :: imin, imax, N
        integer, intent(out) :: nnlist(imin:imax), vlist(:) !(# neighbors x part), (i_neighbor)
        real*8, dimension(:,:), intent(in) :: positions
        real*8, intent(in) :: vcutoff
        integer :: i,j,k
        real*8 :: rij, vcutoff2

        vcutoff2=vcutoff*vcutoff

        do i=imin,imax
            nnlist(i)=0
        enddo

        k=1
        do i=imin,imax
            do j=1,N
                rij=(positions(i,1)-positions(j,1))**2+(positions(i,2)-positions(j,2))**2+(positions(i,3)-positions(j,3))**2
                if (rij<vcutoff2) then
                    nnlist(i)=nnlist(i)+1
                    vlist(k)=j
                    k=k+1
                endif
            enddo
        enddo

    end subroutine verletlist
end module integrators

