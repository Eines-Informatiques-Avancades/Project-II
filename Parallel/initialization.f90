! First module of the molecular dynamics simulation program.
! This program implements the initial positions of the particles,
! the initial velocities.

module initial_positions_module
    implicit none
    public :: initial_positions_parallel
    public :: initial_positions_serial
    public :: input_parameters
    public :: initial_velocities
    public :: main
contains

    subroutine initial_positions_parallel(N, L, position, iproc)
        !!! --- Author: Guillem Arasa --- !!!
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


    end subroutine initial_positions_parallel

    subroutine initial_positions_serial(N, L, position)
        !!! --- Author: Guillem Arasa --- !!!
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

    end subroutine initial_positions_serial

    subroutine input_parameters(N, L, T)
        !!! --- Author: Guillem Arasa --- !!!
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
        !!! --- Author: Guillem Arasa --- !!!
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
        !!! --- Authors: Emma Valdés and Guillem Arasa --- !!!
        ! 
        ! Subroutine assigns min and max particle labels assigned to each numbered processor.
        ! Processor with rank 0 will always have the least possible amount of particles.
        ! 
        ! Arguments
        ! nproc     (integer, intent(in)): number of processors used in the simulation run.
        ! N         (integer, intent(in)): number of particles in the simulated system.
        
        ! Returns
        ! subsystems([nproc,2]): array with min and max labels (columns) for each processor (rows)
        integer, intent(in) :: nproc, N
        integer, intent(out) :: subsystems(nproc,2)
        integer :: i,j, Nsub(nproc), Nrest

        Nsub = N/Nproc ! minimum number of particles in subsystem
        Nrest=mod(N,nproc) ! remainder

        ! If N is not a multiple of nproc, the remaining particles are equally distributed among the 
        ! processors starting with the highest rank, until there are no particles left
        if (nrest/=0) then
            do i = nproc-nrest+1, nproc
                Nsub(i) = Nsub(i)+1
            enddo
        endif

        subsystems(nproc,2) = N                 ! label of first particle imin
        subsystems(nproc,1) = N - Nsub(nproc)+1 ! label of last particle imax
        do i=1,nproc-1 
            j = nproc-i
            subsystems(j,2) = subsystems(j+1,1)-1           ! imax depends on imin of processor with rank above
            subsystems(j,1) = subsystems(j,2)-Nsub(j) +1    ! imin depends on imax of same processor and Nsub of same processor
        enddo 


    end subroutine

    subroutine main
        !!! --- Author: Guillem Arasa --- !!!
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
        call initial_positions_serial(N, L, position)
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