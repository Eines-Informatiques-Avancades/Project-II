module Integrator_module
    
    ! parameters
    use pbc_module
    use forces_module
    implicit none

    ! functions and subroutines
    contains
   
    subroutine therm_Andersen(nparts,velocities,nu,sigma)

    !------------------------------------------------------------------------------------------------------------------------------!
    ! Author: Àlex Teruel
    ! About:
    ! The therm_Andersen subroutine applies the Andersen thermostat on a npart system
    ! Input variables:
    !   nparts: total number of particles in the system
    !   nu: probability of accepting a change of velocity
    !   sigma: distribution deviation
    ! Modified input and output variables:
    !   velocities(nparts,3): matrix containing the velocities in each direction for each particle 
    !------------------------------------------------------------------------------------------------------------------------------!

       implicit none
       ! Input and output variables
       integer :: nparts
       double precision :: nu, sigma
       double precision :: velocities(nparts,3)
       ! Internal variables (subroutine)
       double precision :: nu_n, x1, x2, x3, x4
       integer :: n

       do n = 1, nparts
           call random_number(nu_n)
           ! If the random number < nu, new velocities are calculated
           if (nu_n < nu) then
               ! The gaussians are calculated using the Box-Muller method
               call random_number(x1); call random_number(x2); call random_number(x3); call random_number(x4)
               velocities(n,1) = sigma*dsqrt(-2.d0*(dlog(1.d0-x1)))*dcos(2.d0*pi*x2)
               velocities(n,2) = sigma*dsqrt(-2.d0*(dlog(1.d0-x1)))*dsin(2.d0*pi*x2)
               velocities(n,3) = sigma*dsqrt(-2.d0*(dlog(1.d0-x3)))*dcos(2.d0*pi*x4)
           end if
       end do

    end subroutine therm_Andersen

    subroutine velocity_verlet_andersen(nparts,box_L,cutoff,nu,sigma,dt,positions,velocities,pot_E,forces,therm_on)

     !-----------------------------------------------------------------------------------------------------------------------------!
     ! Author: Àlex Teruel
     ! About:
     ! The velocity_verlet_andersen subroutine applies the velocity Verlet algorithm to integrate an instant of time through an
     ! interacting particle under a Lennard-Jones potential and an Andersen thermostat.
     ! Input variables:
     !   nparts: total number of particles in the system
     !   box_L: length of each side of the box
     !   cutoff: interaction radius
     !   nu: probability of accepting a change of velocity
     !   sigma: distribution deviation
     !   dt: time difference
     !   therm_on: if it is 1, we apply the thermostat, if not, no.
     ! Modified input and output variables:
     !   positions(npart,3): matrix containing the positions of the npart system
     !   velocities(nparts,3): matrix containing the velocities in each direction for each particle 
     !   pot_E: potential energy of the system
     !   forces(nparts,3): matrix to obtain the force applied on each particle of the system
     !-----------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        ! Input and output variables
        integer :: nparts, therm_on
        double precision :: box_L, cutoff, nu, sigma, dt
        double precision :: positions(nparts,3), velocities(nparts,3)
        double precision :: pot_E, forces(nparts,3)
        ! Internal variables (subroutine)
        double precision :: newforces(nparts,3)
    
        ! New positions are calculated
        positions = positions + velocities*dt + 0.5d0*forces*dt**2
        ! Periodic Boundary Conditions are applied
        call pbc_nparts(nparts,positions,box_L)
        ! New forces and potential energy are calculated
        call LJ_potential(nparts,positions,cutoff,box_L,1,pot_E,newforces)
        ! New velocities are calculated
        velocities = velocities + 0.5d0*dt*(forces + newforces)
        ! Force values are reassigned
        forces = newforces
        ! Thermostat is applied (sometimes)
        if(therm_on == 1) call therm_Andersen(nparts,velocities,nu,sigma)
    
        return
    
    end subroutine velocity_verlet_andersen

    subroutine euler_andersen(nparts,box_L,cutoff,nu,sigma,dt,positions,velocities,pot_E,forces,therm_on)

     !-----------------------------------------------------------------------------------------------------------------------------!
     ! Author: Àlex Teruel
     ! About:
     ! The euler_andersen subroutine applies the Euler algorithm to integrate an instant of time through an
     ! interacting particle under a Lennard-Jones potential and an Andersen thermostat.
     ! Input variables:
     !   nparts: total number of particles in the system
     !   box_L: length of each side of the box
     !   cutoff: interaction radius
     !   nu: probability of accepting a change of velocity
     !   sigma: distribution deviation
     !   dt: time difference
     !   therm_on: if it is 1, we apply the thermostat, if not, no.
     ! Modified input and output variables:
     !   positions(npart,3): matrix containing the positions of the npart system
     !   velocities(nparts,3): matrix containing the velocities in each direction for each particle 
     !   pot_E: potential energy of the system
     !   forces(nparts,3): matrix to obtain the force applied on each particle of the system
     !-----------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        ! Input and output variables
        integer :: nparts, therm_on
        double precision :: positions(nparts,3), velocities(nparts,3), forces(nparts,3), dt, box_L, cutoff, pot_E, nu, sigma
    
        ! new positions
        positions = positions + velocities*dt + 0.5d0*forces*dt**2
        ! new velocities
        velocities = velocities + forces*dt
        ! periodic boundary conditions
        call pbc_nparts(nparts,positions,box_L)
        ! new forces and potential energy
        call LJ_potential(nparts,positions,cutoff,box_L,1,pot_E,forces)
        ! thermostat (sometimes)
        if (therm_on == 1) call therm_Andersen(nparts,velocities,nu,sigma)

        return
    
    end subroutine euler_andersen

end module Integrator_module