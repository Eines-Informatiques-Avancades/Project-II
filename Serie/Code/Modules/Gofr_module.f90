module pair_distribution_function_module

    use pbc_module
    implicit none

    contains

! pair distribution function
subroutine pair_distribution_function(class,nparts,measures,rho,positions,nboxes,gofr,length,cutoff,dr,results)

    ! information:
    !   case 1: initialize gofr
    !   case 2: sample distances
    !   case 3: normalize gofr
    !   case 4: write results

    implicit none
    ! in/out variables
    character(len=50), intent(in)   :: results
    integer, intent(in)             :: class, nparts, nboxes
    double precision, intent(in)    :: rho, positions(nparts,3), length, cutoff
    integer, intent(inout)          :: measures
    double precision, intent(inout) :: gofr(nboxes), dr
    !---------------------------------------------------------------------------!
    ! internal variables
    integer :: i, j, xi
    double precision :: x1(3), x2(3), xr, vol
    double precision, parameter :: pi = 4.d0*datan(1.d0)

    select case (class)

        case (1)

        measures = 0
        gofr = 0.d0
        dr = cutoff/dble(nboxes)

        case (2)

        measures = measures + 1
        do i = 1, nparts-1
            do j = i+1, nparts
                ! calculate the distance between particles i and j with pbc
                x1 = positions(i,:)
                x2 = positions(j,:)
                xr = distance(x1,x2,length)
                ! if the distance is less than L/2 add 1 in the corresponding box
                if (xr < cutoff) then
                    xi = int(xr/dr) + 1
                    gofr(xi) = gofr(xi) + 2
                endif
            enddo
        enddo

        case (3)

        do i = 1, nboxes
            vol = 4.d0/3.d0*pi*dr**3*dble(i**3-(i-1)**3)
            gofr(i) = gofr(i)/(dble(measures*nparts)*vol*rho)
        end do

        case (4)

        open(99,file=results)
        do i = 1, nboxes
            xr = dr*i - 0.5d0*dr
            write(99,*) xr, gofr(i)
        end do
        close(99)

        case default

        print*, "Pair distribution function: Wrong case. Case must be 1, 2, 3, or 4."

    end select

    return

end subroutine pair_distribution_function

! distance between 2 particles with pbc
double precision function distance(x1,x2,l)

implicit none
double precision, intent(in) :: x1(3), x2(3), l
double precision, external   :: pbc_func
double precision             :: x12(3)

x12 = x1-x2
call pbc(x12,l)
distance = dsqrt(sum(x12**2))

return

end function distance

end module pair_distribution_function_module