module forces_module

use Constants_module
use pbc_module
implicit none

contains

subroutine LJ_potential(npart,positions,cutoff,length,pbc_on, Upot, force)

    !------------------------------------------------------------------------------------------------------------------------------!
    ! Author: Daniel Conde
    ! About:
    ! The LJ_potential subroutine calculates the potential energy and the forces between particles, which are subjected to a Lennard-
    ! Jones potential.
    ! Input variables:
    !   npart: total number of particles in the system
    !   positions(npart,3): matrix containing the positions of all the particles (npart)
    !   cutoff: interaction radius
    !   length: length of each side of the box
    !   pbc_on: if 1, periodic boundary conditions are applied
    ! Output variables:
    !   Upot: potential energy of the system
    !   force(npart,3): matrix to obtain the force applied on each particle
    !------------------------------------------------------------------------------------------------------------------------------!

    implicit none
    integer, intent(in) :: npart, pbc_on
    double precision, intent(in) :: positions(npart,3), cutoff, length
    double precision :: Upot, force(npart,3)
    integer :: i, j, k
    double precision :: dr(3), dr2, dr6, dr8, dr12, dr14
   
    ! Energy is initialized
    Upot = 0.0
    ! Force is initialized
    force = 0.d0
    ! Loop that loops through all interactions
     !$omp parallel private(dr,dr6,dr8, dr12, dr14,dr2) 
!$omp do schedule(dynamic,4)  reduction(+:Upot) reduction(+:force)
!algo mal en la parelizacion? No va más rápido pero parece que da los mismos resultados, probar con sistema más grande y más tiempo
!La actual version del programa no funciona, el error parece estar en el main (pues una anterior iba)
    do i = 1, npart-1
        do j = i+1, npart
        ! Difference between two particles calculation
            do k = 1, 3
                dr(k) = positions(i,k) - positions(j,k)
            end do
            ! PBC are applied
            if(pbc_on == 1) call pbc(dr,length)
            ! Distance calculation
            dr2 = dr(1)**2 + dr(2)**2 + dr(3)**2
            ! Cutoff correction
            if(dr2 <= cutoff**2) then
                dr6 = dr2**3
                dr8 = dr2**4
                dr12 = dr2**6
                dr14 = dr2**7
                ! Potential energy
                Upot = Upot+4.d0*(1.d0/dr12 - 1.d0/dr6) - 4.d0*( 1.d0/cutoff**12 - 1.d0/cutoff**6)
                ! Force particle i
                force(i,1) = force(i,1) + (48.d0/dr14 - 24.d0/dr8)*dr(1)
                force(i,2) = force(i,2) + (48.d0/dr14 - 24.d0/dr8)*dr(2)
                force(i,3) = force(i,3) + (48.d0/dr14 - 24.d0/dr8)*dr(3)
                ! Force particle j
                force(j,1) = force(j,1) - (48.d0/dr14 - 24.d0/dr8)*dr(1)
                force(j,2) = force(j,2) - (48.d0/dr14 - 24.d0/dr8)*dr(2)
                force(j,3) = force(j,3) - (48.d0/dr14 - 24.d0/dr8)*dr(3)
            end if
        end do
    end do
    !$omp end do
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

    !------------------------------------------------------------------------------------------------------------------------------!
    ! Author: Daniel Conde
    ! About:
    ! The pressure subroutine calculates the pressure of a npart system.
    ! Input variables:
    !   n: number of total particles
    !   rho: density of the system
    !   temp: temperature of the system
    !   l: length of the box
    !   cutoff: interaction radius
    !   pos: position of each particle
    ! Output variables:
    !   pres: system pressure
    !------------------------------------------------------------------------------------------------------------------------------!

    implicit none
    ! Variables
    integer :: n
    double precision :: rho, temp, l, cutoff, pos(n,3)
    double precision :: pres
    integer :: i, j
    double precision :: dr(3), dr2, force(3)

    pres = 0.d0
    do i=1,n-1
        do j=i+1,n
            dr = pos(i,:)-pos(j,:)
            call pbc(dr,l)
            dr2 = sum(dr**2)
            if (dr2 < cutoff**2) then
                force = (48.d0/dr2**7 - 24.d0/dr2**4)*dr
                pres = pres + sum(dr*force)
            end if
        end do
    end do
    pres=pres+16/3.*pi*rho**2*(2/3.*(1./cutoff)**9)-(1/cutoff)**3/(1/(3*l**3))! Added for cut-off corrections
    pres = rho*temp + pres/(3.d0*l**3)

    return

end subroutine pressure


end module forces_module
