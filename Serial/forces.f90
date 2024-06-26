module Forces_and_Energies
use pbc_mod 
contains
    subroutine Pressure (positions,boxsize,cutoff,temp,press)
        !!! --- Author: Quim Badosa --- !!!

        ! This function calculates the pressure done by a number of particles (given by the number of position elements) in a
        ! box of a certain size at a certain temperature. The cutoff is used to set an interacting range of the particles.

        ! Args:
        ! positions (double precision, dimension(:,:)): positions of the particles in the box
        ! vlist (integer, dimension(:)): list of particles that are inside the cutoff distance of a particle
        ! nnlist (integer, dimension(:)): number of neighbors of each particle
        ! imin (integer): index of the first particle in the positions array
        ! imax (integer): index of the last particle in the positions array
        ! cutoff (double precision): cutoff distance for the interaction between particles
        ! boxsize (double precision): size of the box

        ! Returns:
        ! Virialterm (double precision): virial contribution term of the pressure
        ! max_dist (double precision): maximum distance between particles
       
        implicit none

        ! IdealGas contribution term : N*Kb*T/V
        ! Virial contribution term : (1/3V)*sum(r_ij*f_ij)  

        Double precision,allocatable, dimension(:,:), intent(in) :: positions
        Double precision, intent(in) :: cutoff,temp,boxsize
        Double precision, intent(out) :: press
        Double precision, dimension(3,1) :: r_ij
        Double precision, dimension(3) :: f_ij
        Double precision :: d_ij, volume, Virialterm, cf2
        Integer ::  npart,i,j

        ! Initial values:
        Press=0.d0
        Virialterm=0.d0
        !Volume and number of particles:
        volume= dble(boxsize*boxsize*boxsize)
        npart= int(size(positions,dim=1))
        cf2 = cutoff*cutoff
        !Force between particles:
        do i=1,npart-1
            do j=i+1,npart
                !Distance between particles:
                r_ij(1,1)=positions(i,1)-positions(j,1)
                r_ij(2,1)=positions(i,2)-positions(j,2)
                r_ij(3,1)=positions(i,3)-positions(j,3)

        call minimum_image(r_ij(1,1), boxsize)
        call minimum_image(r_ij(2,1), boxsize)
        call minimum_image(r_ij(3,1), boxsize)

        !Computes distance squared
        d_ij=(r_ij(1,1)**2.d0)+(r_ij(2,1)**2.d0)+(r_ij(3,1)**2.d0)
        !Now we compare this distance with the cutoff

                if (d_ij< cf2) then
                    f_ij(1) = ((48.d0 / (d_ij**7)) -( 24.d0 / (d_ij**4))) * r_ij(1,1)
                    f_ij(2) = ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(2,1)
                    f_ij(3) = ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(3,1)
                    !redefine virial Term value:
                    Virialterm = Virialterm + dot_product(r_ij(:,1),f_ij)
                end if
            end do
        end do
    
    press= (dble(npart)*temp)/volume + (1.d0/(3.d0*volume))*Virialterm

    end subroutine Pressure

    subroutine VDW_forces (positions,boxsize,cutoff, VDW_force)
        !!! --- Author: Quim Badosa --- !!!
        ! Subroutine that calculates the interaction between particles using the Lennard Jones potential
        ! The interacting range is limited by a cutoff distance.          
        implicit none
        !ARGUMENTS:
        !Positions : Positions of the particles DIM= (d,npart)  (d usually=3)
        !Boxsize: We will supose cubic system, double precision
        !Cutoff: A range of interaction real double precision 
        !VDW_force: Total interaction force                       
        Double precision,allocatable, dimension(:,:),intent(inout) :: positions
        double precision, allocatable, dimension(:,:), intent(inout) :: vdw_force
        Double precision, intent(in) :: cutoff,boxsize
        !VARIABLES:
        !r_ij : relative position vector between pair of particles, double precision dim= (3,1)
        !d_ij : distance module between the particles, double precision 
        !f_ij: Force between the particles
        !npart: particle number. Integer
        Double precision, dimension(3,1) :: r_ij=0.d0
        Double precision :: d_ij,cf2
        Integer ::  npart,i,j
        npart= int(size(positions,dim=1))
        cf2 = cutoff*cutoff
        vdw_force = 0.d0
        do i=1,npart-1
            do j=i+1,npart
                        
                !Distance between particles:
                r_ij(1,1)=positions(i,1)-positions(j,1)
                r_ij(2,1)=positions(i,2)-positions(j,2)
                r_ij(3,1)=positions(i,3)-positions(j,3)
    
                call minimum_image(r_ij(1,1), boxsize)
                call minimum_image(r_ij(2,1), boxsize)
                call minimum_image(r_ij(3,1), boxsize)
                !Module of r_ij
                !Compute distance squared
                d_ij=(r_ij(1,1)*r_ij(1,1))+(r_ij(2,1)*r_ij(2,1))+(r_ij(3,1)*r_ij(3,1))
                !Now we compare this distance with the cutoff
                if (d_ij< cf2) then
                    !Force made by j to i
                    vdw_force(i,1) = vdw_force(i,1) + ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(1,1)
                    vdw_force(i,2) = vdw_force(i,2) + ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(2,1)
                    vdw_force(i,3) = vdw_force(i,3) + ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(3,1)
        
                    !Force made by i to j
                    vdw_force(j,1) = vdw_force(j,1) - ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(1,1)
                    vdw_force(j,2) = vdw_force(j,2) - ((48.d0 / (d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(2,1)
                    vdw_force(j,3) = vdw_force(j,3) - ((48.d0 /( d_ij**7)) - (24.d0 / (d_ij**4))) * r_ij(3,1)
    
                end if
            end do
        end do
    end subroutine VDW_forces


    subroutine potentialE (positions,cutoff,PotentialEn, boxsize)
        !!! --- Author: Quim Badosa --- !!!
        !Subroutine to calculate the potential energy between particles using Lennard-Jones potential
        implicit none
        !ARGUMENTS:
        !Positions : Positions of the particles DIM= (d,npart)  (d usually=3)
        !Boxsize: We will supose cubic system, double precision
        !Cutoff: A range of interaction real double precision
        !PotentialE: Total potential energy
        Double precision,allocatable, dimension(:,:), intent(in) :: positions
        Double precision, intent(in) :: cutoff, boxsize
        Double precision, intent(out) :: PotentialEn
        !VARIABLES:
        !r_ij : relative position vector between pair of particles, double precision dim= (3,1)
        !d_ij : distance module between the particles, double precision
        !f_ij: Force between the particles
        !npart: particle number. Integer
        Double precision, dimension(3,1) :: r_ij
        Double precision:: e_ij,d_ij2, cf2
        Integer ::  npart,i,j

        potentialEn=0.d0
        npart= int(size(positions,dim=1))
        cf2 = cutoff*cutoff
        !call PBC(positions, boxsize, npart)
        do i=1,npart-1
            do j=i+1,npart
                !Distance between particles:
                r_ij(1,1)=positions(i,1)-positions(j,1)
                r_ij(2,1)=positions(i,2)-positions(j,2)
                r_ij(3,1)=positions(i,3)-positions(j,3)

                call minimum_image(r_ij(1,1), boxsize)
                call minimum_image(r_ij(2,1), boxsize)
                call minimum_image(r_ij(3,1), boxsize)
                !THIS WAY WE MAKE SURE THAT r_ij is inside our box
                !Module of r_ij
                d_ij2=(r_ij(1,1)**2)+(r_ij(2,1)**2)+(r_ij(3,1)**2)
                !Now we compare this distance with the cutoff
                if (d_ij2 < cf2) then
                    !d6 = d_ij2*d_ij2*d_ij2 
                    !d12 = d6*d6 
                    !e_ij = 4.d0*((1.d0/d12) - (1.d0/d6)) - 4.d0*((1.d0/(cf2**6)) - (1.d0/(cf2**3)))
                    e_ij = - 4.d0*(1.d0/cutoff**12.d0 - 1.d0 /cutoff**6.d0)+ 4*(1.d0/d_ij2**6.d0-1.d0/d_ij2**3.d0)
                    potentialEn = potentialEn + e_ij
                end if
            end do
        end do
    end subroutine potentialE

    subroutine  kineticE (velocities,KineticEn)
        !!! --- Author: Quim Badosa --- !!!
        !!! --- Contributor: Paula Sierra --- !!!

        !Subroutine to calculate the kinetic energy of the particles
        implicit none
        !ARGUMENTS:
        !Velocities : Velocities of the particles DIM= (npart,d)  (d usually=3)
        !KineticEn: Total kinetic energy
        Double precision,allocatable, dimension(:,:),intent(in) :: velocities
        Double precision, intent(out) :: KineticEn

        !VARIABLES:
        !npart: particle number. Integer
        ! Double precision, dimension(3,1) :: v_ij
        Double precision:: e_ij
        Double precision :: vel_ij
        Integer ::  npart,i
        npart= int(size(velocities,dim=1))
        KineticEn=0.d0
        do i=1,npart
            vel_ij=dsqrt((velocities(i,1)**2.d0)+(velocities(i,2)**2.d0)+(velocities(i,3)**2.d0))
            e_ij= (1.d0/2.d0)*vel_ij*vel_ij
            KineticEn =KineticEn + e_ij
        enddo
    end subroutine kineticE

    subroutine Tempinst (KE,npart,Tinst)
        !!! --- Author: Quim Badosa --- !!!
        ! SUbroutine to calculate the instant Temperature
        ! Kinetic energy has to be called earlier, so Kinetic energy and number of particles are the parameters of this
        ! function

        implicit none
        !ARGUMENTS
        !Kinetic energy. double precision
        !Number of particles. Integer
        Double precision, intent(in) :: KE
        Integer, intent(in) :: npart
        !InstantTemperature, Tinst     
        Double precision, intent(out) :: Tinst
        Tinst =(2.d0/(3.d0*real(npart)-3.d0))*ke

    end subroutine Tempinst

end module Forces_and_Energies

