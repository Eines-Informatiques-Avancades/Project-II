!---------------------------------------------------------Initial state module----------------------------------------------------------!
module Initial_state_module

    implicit none
    
    contains

    subroutine sc_lattice(npart,density,positions,length,space,celdim)

        !------------------------------------------------------------------------------------------------------------------------------!
        ! Author: Oliver Loveday
        ! About:
        ! The sc_lattice subroutine generates a crystal system (disposition-like) with a simple cubic unit cell, centered at the 
        ! origin of coordinates.
        ! Input variables:
        !   npart: total number of particles in the system
        !   density: density
        ! Output variables:
        !   positions(npart,3): matrix containing the positions of all the particles (npart)
        !   length: length of each side of the box
        !   space: distance between two nodes in the same dimension
        !   celdim: number of nodes in each dimension
        !------------------------------------------------------------------------------------------------------------------------------!
        
        implicit none

        integer, intent(in) :: npart
        double precision, intent(in) :: density
        integer, intent(out) :: celdim
        double precision, intent(out) :: positions(npart,3), length
        double precision :: space
        integer :: i, j, k, numpart
        
        ! Parameters are determined
        length = (dble(npart)/density)**(1.d0/3.d0)
        celdim = int(npart**(1.d0/3.d0)+0.5d0)
        space = length/dble(celdim)

        ! Particles are introduced in each position
        numpart = 0
        do i = 1, celdim
            do j = 1, celdim
                do k = 1, celdim
                    numpart = numpart + 1
                    ! Counter is updated each time a new particle is introduced
                    positions(numpart,1) = -0.5d0*length + dble(i-1)*space
                    positions(numpart,2) = -0.5d0*length + dble(j-1)*space
                    positions(numpart,3) = -0.5d0*length + dble(k-1)*space
    
                end do
            end do
        end do
    
        return
    
    end subroutine sc_lattice

    subroutine bimodal_vel(nparts,temp,vel)

        !------------------------------------------------------------------------------------------------------------------------------!
        ! Author: Oliver Loveday
        ! Information
        ! The bimodal_vel initializes the velocities of a particle system following the bimodal distribution, which is compatible
        ! for a given temperature
        ! Input variables:
        !   nparts: total number of particles in the system
        !   temp: temperature of the system
        ! Input and output variables (modified):
        !   vel(nparts,3): matrix which contains the velocities of the particles for each direction (x, y and z)
        !------------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        ! Input and output variables
        integer :: nparts
        double precision :: temp, vel(nparts,3)
        ! Other variables
        integer :: ipos, ineg, i, j
        double precision :: vpos, vneg, xi
    
        vpos = dsqrt(temp)
        vneg = -vpos
        ipos = 0
        ineg = 0
        do i = 1, nparts
            call random_number(xi)
            if ((xi < 0.5d0).and.(ipos < nparts/2)) then
                ipos = ipos + 1
                do j = 1, 3
                    vel(i,j) = vpos
                end do
            else if (ineg < nparts/2) then
                ineg = ineg + 1
                do j = 1, 3
                    vel(i,j) = vneg
                end do
            else
                ipos = ipos + 1
                do j = 1, 3
                    vel(i,j) = vpos
                end do
            end if
        end do
    
        return
    
    end subroutine bimodal_vel

end module Initial_state_module
