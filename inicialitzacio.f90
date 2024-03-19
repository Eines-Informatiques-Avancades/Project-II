! First module of the molecular dynamics simulation program.
! This program implements the initial positions of the particles,
! the initial velocities.

! Subrutina inicializar posiciones en una red cubica
module initial_positions_module
    implicit none
contains

    subroutine initial_positions(N, L, position)
        implicit none
        integer, intent(in) :: N
        real*8, intent(in) :: L
        real(8), dimension(N,3), intent(out) :: position
        real(8) :: a
        integer :: i, j, k, M

        M = NINT( N**(1.0d0/3.0d0) )
        a = L / dble(M)

        do i = 1, M
            do j = 1, M
                do k = 1, M
                    position((i-1)*M**2 + (j-1)*M + k, 1) = a * (i-0.5d0)
                    position((i-1)*M**2 + (j-1)*M + k, 2) = a * (j-0.5d0)
                    position((i-1)*M**2 + (j-1)*M + k, 3) = a * (k-0.5d0)
                end do
            end do
        end do
    end subroutine initial_positions

    ! Subrutina para leer los parametros de entrada: N, L, T
    subroutine input_parameters(N, L, T)
        implicit none
        integer, intent(out) :: N
        real(8), intent(out) :: T,L
        integer :: i
        character(len=20) :: filename
        character(len=20) :: line
        real :: value

        filename = 'input.dat'
        open(unit=1, file=filename, status='old', action='read')

        do i= 1, 3
            read(1, *) line, value
            line = adjustl(line)
            select case(line)
                case('N')
                    N = int(value)
                case('L')
                    L = dble(value)
                case('T')
                    T = value
                case default
                    print *, 'Error: unknown parameter'
            end select
        end do
        close(1)
    end subroutine input_parameters

    ! Initial velocities
    subroutine initial_velocities(N, T, velocity)
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: T
        real(8), dimension(N,3), intent(out) :: velocity
        real(8) :: v, vcm(3)
        integer :: i, j
        real(8), dimension(3) :: gaussian

        call random_number(gaussian)
        do i = 1, 3
            vcm(i) = 0.0d0
            do j = 1, N
                vcm(i) = vcm(i) + gaussian(i)
            end do
            vcm(i) = vcm(i) / dble(N)
        end do

        do i = 1, N
            call random_number(gaussian)
            do j = 1, 3
                velocity(i,j) = gaussian(j) - vcm(j)
            end do
            v = 0.0d0
            do j = 1, 3
                v = v + velocity(i,j)**2
            end do
            v = sqrt(v)
            velocity(i,:) = velocity(i,:) / v * sqrt(T)
        end do

        open(3, file='initial_velocities.dat')
        do i = 1, N
            write(3,*) velocity(i,1), velocity(i,2), velocity(i,3)
        end do
        close(3)
    end subroutine initial_velocities

    subroutine main
        implicit none
        integer :: N
        real(8) :: T,L
        integer :: i
        real(8), dimension(:,:), allocatable :: position

        call input_parameters(N, L, T)
        print *, 'N = ', N
        print *, 'L = ', L
        print *, 'T = ', T
        allocate(position(N,3))
        call initial_positions(N, L, position)
        open(2, file='initial_positions.xyz')
        do i = 1, N
            write(2,*) position(i,1), position(i,2), position(i,3)
        end do
        close(2)
        deallocate(position)
        allocate(position(N,3))
        call initial_velocities(N, T, position)
        deallocate(position)
    end subroutine main

end module initial_positions_module