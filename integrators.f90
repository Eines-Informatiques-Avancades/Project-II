module integrators
    use module_1
    use module_2
    ! etc etc
    ! obviament de moment no compila
    implicit none
    public :: vv_integrator

contains
    subroutine vv_integrator(positions, velocities, cutoff, L, dt)
        !
        !  Subroutine to update the positions of all particles using the Velocity Verlet
        ! algorithm. 

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
        real*8, dimension(:,:), intent(inout)      :: positions, velocities
        real*8, intent(in)                         :: cutoff, L, dt
        ! local variables
        real*8, dimension(3, size(positions(1,:))) :: forces

        call calc_vdw_force(positions, cutoff, L, forces)
        
        positions = positions + (dt*velocities) + (0.5d0*dt*dt*forces)
        call PBC(positions, L)

        velocities = velocities + (0.5d0*dt*forces)

        call calc_vdw_force(positions, cutoff, L, forces)
        velocities = velocities + 0.5d0*dt*forces

    end subroutine vv_integrator

end module integrators

