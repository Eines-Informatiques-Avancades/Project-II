module pbc_module

    !parameters
    implicit none

    !functions and subroutines
    contains

    subroutine pbc(position,length)

        !------------------------------------------------------------------------------------------------------------------------------!
        ! Author: Adrià Calzada
        ! About:
        ! The pbc subroutine applies the periodic boundary conditions on one particle
        ! Input variables:
        !   length: length of each side of the box
        ! Modified input and output variables:
        !   position(3): vector containing the position of one particle
        !------------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        double precision, intent(in) :: length
        double precision :: position(3)
        integer :: i
       
        do i = 1, 3
            if(position(i) >  0.5d0*length) position(i) = position(i) - length
            if(position(i) < -0.5d0*length) position(i) = position(i) + length
        end do
    
        return
       
    end subroutine pbc
    
    subroutine pbc_nparts(nparts,positions,length)
    
        !------------------------------------------------------------------------------------------------------------------------------!
        ! Author: Adrià Calzada
        ! About:
        ! The pbc_nparts subroutine applies the periodic boundary conditions (pbc) on all the particles of the system.
        ! Input variables:
        !   nparts: total number of particles in the system
        !   length: length of each side of the box
        ! Modified input and output variables:
        !   position(nparts,3): matrix containing the positions of the npart system
        !------------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        integer, intent(in) :: nparts
        double precision, intent(in) :: length
        double precision :: positions(nparts,3)
        integer :: i, n
    
        do n = 1, nparts
            do i = 1, 3
                if(positions(n,i) >  0.5d0*length) positions(n,i) = positions(n,i) - length
                if(positions(n,i) < -0.5d0*length) positions(n,i) = positions(n,i) + length
            end do
        end do
    
        return
    
    end subroutine pbc_nparts
end module pbc_module
