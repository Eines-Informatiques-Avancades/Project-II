module gr_module
    implicit none

contains 
    interface
        subroutine minimum_image(r, L)
            double precision, intent(inout) :: r
            double precision, intent(in) :: L
        end subroutine minimum_image
    end interface

    subroutine calculate_g_r(positions, L, n_particles, dr, g_r, max_r, n_bins)
        implicit none
        double precision, allocatable, dimension(:,:) :: positions
        double precision, intent(in) :: L, dr, max_r
        integer, intent(in) :: n_particles, n_bins
        double precision, allocatable, dimension(:) :: g_r
        integer :: i, j, bin
        double precision :: r, delta_r
        double precision, parameter :: pi = 3.1415926535897932384626433832795


        ! Initialize g(r) array
        if (.not. allocated(g_r)) then
            allocate(g_r(n_bins))
        end if
        g_r = 0.0

        do i = 1, n_particles-1
            do j = i+1, n_particles
                r = sqrt(sum((positions(i,:) - positions(j,:))**2))
                
                ! Apply minimum image convention
                call minimum_image(r, L)
                
                if (r <= max_r) then
                    bin = int(r / dr) + 1
                    if (bin >= 1 .and. bin <= n_bins) then
                        g_r(bin) = g_r(bin) + 2
                    end if
                end if
            end do
        end do
        
        ! Normalize g(r)
        do bin = 1, n_bins
            delta_r = (bin - 0.5) * dr
            g_r(bin) = g_r(bin) / (n_particles * (4.0/3.0 * pi * (delta_r**3 - (delta_r-dr)**3)))
        end do
    end subroutine calculate_g_r
end module gr_module