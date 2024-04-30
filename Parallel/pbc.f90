module pbc_mod
    implicit none
    public :: PBC
    public :: minimum_image
contains

    subroutine PBC(imin, imax, position, L, N)
        !
        ! This subroutine applies the periodic boundary conditions to the positions of the particles,
        ! so that the particles are always inside the box. The box is centered at the origin.
        !
        ! Args:
        !    position:        (N,3): array with the positions of the particles
        !    L:               double precision with the size of the box
        !
        ! Returns:
        !    position:        (N,3): array with the positions of the particles with the PBC applied
        integer, intent(in) :: imin, imax
        double precision, dimension(:,:), intent(inout) :: position
        double precision, intent(in) :: L
        integer :: i, j, M, N
        double precision :: center

        center = L/ 2.0d0

        N = size(position(:,1))
        M = size(position(1,:))

        do i = imin, imax
            do j = 1, M
                if (position(i,j) > center + L/2.0d0) then
                    position(i,j) = position(i,j) - L
                else if (position(i,j) < center - L/2.0d0) then
                    position(i,j) = position(i,j) + L
                end if
            end do
        end do
        
    end subroutine PBC

    subroutine minimum_image(dx, L)
        !
        ! This subroutine calculates the minimum image convention for the distance between two particles.
        ! Args:
        !   dx:             double precision, intent(inout) with the distance between two particles
        !   L:              double precision with the size of the box
        !
        ! Returns:
        !   dx:             double precision, intent(inout) with the minimum image convention applied

        implicit none
        double precision, intent(inout) :: dx
        double precision, intent(in) :: L

        if (dx > L/2.0d0) then
            dx = dx - L
        else if (dx < -L/2.0d0) then
            dx = dx + L
        end if

    end subroutine minimum_image

end module pbc_mod
