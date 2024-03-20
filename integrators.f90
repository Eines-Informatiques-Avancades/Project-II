module integrators
    use pbc_mod
    use Forces_and_Energies

    implicit none
    public :: vv_integrator

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
        real*8, allocatable, dimension(:,:), intent(inout)      :: positions, velocities, forces
        real*8, intent(in)                         :: cutoff, L, dt
        integer :: N
        N = size(positions,dim=1)

        call VDW_forces(positions, L, cutoff, forces)
        positions = positions + (dt*velocities) + (0.5d0*dt*dt*forces)
        call PBC(positions, L,N)

        velocities = velocities + (0.5d0*dt*forces)

        call VDW_forces(positions, L, cutoff, forces)
        velocities = velocities + 0.5d0*dt*forces
    end subroutine vv_integrator

    subroutine main_loop(N_steps, N_save_pos, dt, L, cutoff, positions, velocities)
    implicit none
    integer, intent(in) :: N_steps, N_save_pos
    real*8, intent(in) :: dt, cutoff,L
    real*8, allocatable, dimension(:,:), intent(inout) :: positions, velocities 
    real*8, allocatable, dimension(:,:) :: forces
    integer :: N, i, unit_dyn,j
    real*8 :: time
    N = size(positions, dim=1)
    allocate(forces(N,3))
    ! open files
    open(newunit=unit_dyn,file = 'dynamics.dat',status="REPLACE")
    do i=1,N
        write(unit_dyn,'(3(f8.3,x))') positions(i,:)
    enddo
    write(unit_dyn,'(A)') " "
    ! write initial positions and velocities at time=0
    do i=1,N_steps 
        time = i*dt
        call vv_integrator(positions,velocities,forces,cutoff,L,dt)
        ! write variables to output - positions, energies
        if (MOD(i,N_save_pos).EQ.0) then
            do j=1,N 
                write(unit_dyn,'(3(e12.3,x))') positions(j,:)
            enddo 
            write(unit_dyn,'(A)') " "
        endif
    enddo
    deallocate(forces)
    close(unit_dyn)
    end subroutine main_loop
end module integrators

