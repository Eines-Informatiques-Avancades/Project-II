module forces_module

use Constants_module
use pbc_module
!--------------------------!
use mpi_vars
!--------------------------!
implicit none

contains

subroutine LJ_potential(npart,positions,cutoff,length,pbc_on, Upot, force)
    implicit none
    integer, intent(in) :: npart, pbc_on
    double precision, intent(in) :: positions(npart,3), cutoff, length
    double precision :: Upot, force(npart,3)
    integer :: i, j, k
    double precision :: dr(3), dr2, dr6, dr8, dr12, dr14
    
   
   
    forces(first_part:last_part,:)=0.d0   
    do i=first_part,last_part
    do j = i+1, i-1
                dr(1) = (/positions(i,1) - positions(j,1),&
                & (positions(i,2) - positions(j,2)),&
                & positions(i,3) - positions(j,3)/)
            end do
            ! PBC are applied
            if(pbc_on == 1) call pbc(dr,length)
            ! Distance calculation
            dr2 = dr(1)**2 + dr(2)**2 + dr(3)**2
            if(dr2 <= cutoff**2) then
                force(i,:) = force(i,:) + (48.d0/dr2**7 - 24.d0/dr2**4)*dxyz
            end if
        end do
    end do
    !$omp ddend do
    !$omp end parallel
    return 
end subroutine LJ_potential

subroutine Kinetic_Energy(nparts,velocity,kin_E)

    !------------------------------------------------------------------------------------------------------------------------------!
    ! Author: Daniel Conde
    ! About:
    ! The Kinetic_Energy subroutine calculates the kinetic energy for a particle system.
    ! Input variables:
    !   nparts: total number of particles in the system
    !   velocity(nparts,3): matrix containing the velocities in each direction for each particle 
    ! Output variables:
    !   kin_E: kinetic energy of the system
    !------------------------------------------------------------------------------------------------------------------------------!

    implicit none
    ! Input and output variables
    integer, intent(in) :: nparts
    double precision, intent(in) :: velocity(nparts,3)
    double precision :: kin_E
    ! Internal variables (subroutine)
    integer :: n, i
    double precision :: vel2
   
    kin_E = 0.d0
    do n = 1, nparts
        vel2 = 0.d0
        do i = 1, 3
            vel2 = vel2 + velocity(n,i)**2
        end do
        kin_E = kin_E + 0.5d0*vel2
    end do
   
    return
   
end subroutine Kinetic_Energy





subroutine pressure(n,rho,temp,l,cutoff,pos,pres)
    implicit none
    ! Variables
    integer :: n
    double precision :: rho, temp, l, cutoff, pos(n,3)
    double precision :: pres,UPot
    integer :: i, j
    double precision :: dr(3), dr2, force(3)
    Upot=0
    pres = 0.d0
    do i=first_part,last_part
        do j=1,i-1
                dr(1) = (/positions(i,1) - positions(j,1),&
                & (positions(i,2) - positions(j,2)),&
                & positions(i,3) - positions(j,3)/)
                call pbc(dr,l)
            dr2 = sum(dr**2)
            if (dr2 < cutoff**2) then
                pres = pres+ (48.d0/dr2**6 - 24.d0/dr2**3)
                Upot=Upot+4.d0/dr2**2-4.d0/dr2**3- 4.d0*( 1.d0/cutoff**12 - 1.d0/cutoff**6)
            end if
        end do
    end do
    call MPI_REDUCE(pres,pres,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(Upot,Upot,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
    !Eliminamos suma doble:
    if (workerid==master) then
            pres=pres/2.d0
            Upot=Upot/(2.d0*Nparts)
            pres = rho*temp + pres/(3.d0*l**3)
    end if

    return

end subroutine pressure


end module forces_module

