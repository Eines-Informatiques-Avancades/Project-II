module pbc_mod
    implicit none
contains

    subroutine PBC(position, L, N)
        double precision, dimension(:,:), intent(inout) :: position
        double precision, intent(in) :: L
        integer :: i, j, M, N


        N = size(position(:,1))
        M = size(position(1,:))

        do i = 1, N
            do j = 1, M
                if (position(i,j) > L) then
                    position(i,j) = position(i,j) - L
                else if (position(i,j) < 0) then
                    position(i,j) = position(i,j) + L
                end if
            end do
        end do
    end subroutine PBC

    ! a les periodic boundary conditions afegir la convencio de minima imatge (L/2)

end module pbc_mod