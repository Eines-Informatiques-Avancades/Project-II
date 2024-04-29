! First module of the molecular dynamics simulation program.
! This program implements the initial positions of the particles,
! the initial velocities.

module initial_positions_module
    implicit none
    public :: initial_positions
    public :: input_parameters
    public :: initial_velocities
    public :: main
contains

    subroutine initial_positions(N, L, position, iproc)
        ! 
        ! This subroutine initializes the positions of the particles
        ! in a cubic lattice.

        ! Args:
        !   N           (integer, intent(in)):  Number of particles
        !   L           (integer, intent(in))  Length of the sides of the box
        !   position    (real(8), dimension(N,3), intent(out)):  Array with the positions of the particles

        ! Returns:
        !   None

        implicit none
        integer, intent(in) :: N, iproc
        real*8, intent(in) :: L
        real(8), dimension(:,:), intent(out) :: position
        real(8) :: a
        integer :: j, k, M

        position(:,:) = 0.d0

        M = NINT( N**(1.0d0/3.0d0) )
        a = L / dble(M)

        do j = 1, M
            do k = 1, M
                position((j-1)*M + k, 1) = a * (iproc + 0.5d0)
                position((j-1)*M + k, 2) = a * (j-0.5d0)
                position((j-1)*M + k, 3) = a * (k-0.5d0)
            end do
        end do


    end subroutine initial_positions

    subroutine input_parameters(N, L, T)
        !
        ! This subroutine reads the input parameters from a file.
        ! 
        ! Args:
        !   None   
        ! 
        ! Returns:
        !   N   (integer, intent(out)):  Number of particles
        !   L   (integer, intent(out)):  Length of the sides of the box
        !   T   (real(8), intent(out)):  Temperature of the system

        implicit none
        integer, intent(out) :: N
        real(8), intent(out) :: T,L
        integer :: i
        character(len=20) :: filename
        character(len=20) :: line
        real :: value

        filename = 'parameters.nml'
        open(unit=1, file=filename, status='old', action='read')

        do i= 1, 12
            read(1, *) line, value
            line = adjustl(line)
            select case(line)
                case('n_particles =')
                    N = int(value)
                case('L =')
                    L = int(value)
                case('temperature =')
                    T = value
                case default
                    print *, 'Error: unknown parameter'
            end select
        end do
        
        close(1)
    end subroutine input_parameters

    subroutine initial_velocities(N, T, velocity)
        !
        ! This subroutine initializes the velocities of the particles
        ! using a Gaussian distribution.
        !
        ! Args:
        !   N           (integer, intent(in)):  Number of particles
        !   T           (real(8), intent(in)):  Temperature of the system
        !
        ! Returns:
        !   velocity    (real(8), dimension(N,3), intent(out)):  Array with the velocities of the particles

        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: T
        real(8), dimension(N,3), intent(out) :: velocity
        integer :: i, j
        real(8), dimension(3) :: rand
        do i=1, N
            do j=1, 3
                call random_number(rand)
                if (rand(j) > 0.5d0) then
                    velocity(i,j) = dsqrt(T)
                else
                    velocity(i,j) = -dsqrt(T)
                end if
            end do
        end do

        open(3, file='initial_velocities.dat')
        do i = 1, N
            write(3,*) velocity(i,1), velocity(i,2), velocity(i,3)
        end do
        close(3)
    end subroutine initial_velocities

    subroutine assign_subsystem(nproc,N,subsystems)
        integer, intent(in) :: nproc, N
        integer, intent(out) :: subsystems(nproc,2)
        integer :: i, Nsubsystem, Nrest

        Nsubsystem=N/nproc
        Nrest=mod(N,nproc)

        do i=1,nproc
            subsystems(i,1)=1+(i-1)*Nsubsystem !imin
            subsystems(i,2)=subsystems(i,1)+(Nsubsystem-1) !imax
        enddo

        subsystems(nproc,2)=subsystems(nproc,2)+Nrest

    end subroutine

    subroutine main
        !
        ! This subroutine is the main routine of the initializations.
        !
        ! Args:
        !   None
        !
        ! Returns:
        !   None

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
        call initial_positions(N, L, position, 1)
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