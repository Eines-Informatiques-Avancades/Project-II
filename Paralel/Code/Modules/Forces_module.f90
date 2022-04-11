module forces_module

use Constants_module
use pbc_module
!--------------------------!
use paralel
!--------------------------!

implicit none
contains
SUBROUTINE L_J(dr,force,Upot,cutoff)
    IMPLICIT NONE
    REAL*8 dr,cutoff,force,Upot
    force=0d0 !Establecemos la fuerza inicial a 0           
    Upot=0d0  !Establecemos la energia potencial a 0
    IF (dr<cutoff) THEN !Si es menor que el cut-off calculamos la interaccion
      force=(48d0/dr**14d0)-(24d0/dr**8d0)
      Upot=4d0*((1d0/dr**12d0)-(1d0/dr**6d0))
    END IF
    RETURN
END SUBROUTINE L_J




SUBROUTINE LJ_potential(n_particles,positions,cutoff,length,pbc_on, Upot, F) 
! Compute the forces, the potential energy of the
! system and the pressure taking into account PBC
    IMPLICIT NONE
    integer, intent(in):: n_particles,pbc_on
    double precision, intent(in) :: positions(n_particles,3), length
    INTEGER :: i,j,par
    REAL*8 :: cutoff,pot,Upot,pressure
    REAL*8 :: dx,dy,dz,d,ff
    double precision :: dr(3)
    REAL*8, DIMENSION(:,:) :: F
    F=0d0
    Upot=0d0
    pressure=0.0
    !Symmetric Matrix Energy
    DO par=index_matrix2(workerid+1,1), index_matrix2(workerid+1,2)
      i=pairindex(par,1) !Particle i on the interaction
      j=pairindex(par,2) !Particle j on the interaction
      dr(1) = (positions(i,1) - positions(j,1))
      dr(2) = (positions(i,2) - positions(j,2))
      dr(3) =  (positions(i,3) - positions(j,3))
      print*,'dr'
      if(pbc_on == 1) call pbc(dr,length)
      ! Distance calculation
      d = (dr(1)**2 + dr(2)**2 + dr(3)**2)**0.5
      CALL L_J(d,ff,pot,cutoff)
            F(i,1)=F(i,1)+ff*dr(1)
            F(i,2)=F(i,2)+ff*dr(2)
            F(i,3)=F(i,3)+ff*dr(3)
            F(j,1)=F(j,1)-ff*dr(1)
            F(j,2)=F(j,2)-ff*dr(2)
            F(j,3)=F(j,3)-ff*dr(3)
            Upot=Upot+pot
            pressure=pressure+(ff*dr(1)**2d0+ff*dr(2)**2d0+ff*dr(3)**2d0)
    END DO
    !Compute the total interaction force for each particle
    call MPI_ALLREDUCE( F, F,n_particles*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
    ! Compute the sum for all the workers
    call MPI_REDUCE(Upot, Upot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(pressure,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)

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
end module forces_module
