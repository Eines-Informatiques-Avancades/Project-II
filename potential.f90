module Forces_and_Energies
use pbc_mod 
contains
        subroutine Pressure (positions,boxsize,cutoff,temp,press)
                !This function calculates the pressure done by a number of particles (given by the number of position elements) in a
                !box of a certain size at a certain temperature. The cutoff is used to set an interacting range of the particles.
             ! Pressure has to terms: Ideal gas contribution Term (thermal motion) and the virial contribution term (interaction
             ! between particles). The second one is based on  the interaction between pairs of particles
             ! (Lennard-Jones potential)

           ! IdealGas contribution term : N*Kb*T/V
           ! Virial contribution term : (1/3V)*sum(r_ij*f_ij)  
       
                implicit none

       !ARGUMENTS:
       !Positions : Positions of the particles DIM= (d,npart)  (d usually=3)
       !Boxsize: We will supose cubic system, double precision
       !Cutoff: A range of interaction real double precision 
       !Temperature: in Kelvin rea, double precision
       ! Press: Result pressure, double precision

                Double precision,allocatable, dimension(:,:) :: positions
                Double precision :: cutoff,temp,Press,boxsize
                
        !VARIABLES:
        !r_ij : relative position vector between pair of particles, double precision dim= (3,1)
        !d_ij : distance module between the particles, double precision 
        !f_ij: Force between the particles
        !Virialterm : contribution for each pair, double precision
        !Volume : boxsize**3.d0 (double precision)
        !npart : Particles number, integer
                
                Double precision, dimension(3,1) :: r_ij
                Double precision, dimension(3) :: f_ij
                Double precision :: d_ij, volume, Virialterm
                Integer ::  npart,i,j


                       

        ! Initial values:
        Press=0.0
        Virialterm=0.0
        !Volume and number of particles:
        volume= boxsize**3.d0
        npart= int(size(positions,dim=2))

        call PBC(positions, boxsize, npart)
        !Force between particles:
        do i=1,npart-1
                do j=i+1,npart
                        !Distance between particles:
                        r_ij(1,1)=positions(1,i)-positions(1,j)
                        r_ij(2,1)=positions(2,i)-positions(2,j)
                        r_ij(3,1)=positions(3,i)-positions(3,j)

                        call minimum_image(r_ij(1,1), boxsize)
                        call minimum_image(r_ij(2,1), boxsize)
                        call minimum_image(r_ij(3,1), boxsize)

                        !Module of r_ij

                        d_ij=((r_ij(1,1)**2.d0)+(r_ij(2,1)**2.d0)+(r_ij(3,1)**2.d0))**(1.d0/2.d0)
                        !Now we compare this distance with the cutoff

                        if (d_ij< cutoff) then
                                 f_ij(1) = (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(1,1)
                                 f_ij(2) = (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(2,1)
                                 f_ij(3) = (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(3,1)

                         !redefine virial Term value:
                                Virialterm = Virialterm + dot_product(r_ij(:,1),f_ij)
                         end if
                   end do
          end do
          ! Paula: I added a ".d0" to the 10 because, if not added, is considered as integer and gives 0.
          press= (real(npart)*(1.38*10.d0**(-23))*temp)/volume + (1.d0/(3.d0*volume))*Virialterm
          
          ! Pressure units are J/m³
          ! we multiply ideal gas term *Kb=1.38*10**(-23)
          
                  
      end subroutine Pressure
      

        


      !!!!!---------------------------------------------------------------------------------------------------------------------


      subroutine VDW_forces (positions,boxsize,cutoff, VDW_force)
              ! Subroutine that calculates the interaction between particles using the Lennard Jones potential
              ! The interacting range is limited by a cutoff distance.

                            
              implicit none

                !ARGUMENTS:
                   !Positions : Positions of the particles DIM= (d,npart)  (d usually=3)
                   !Boxsize: We will supose cubic system, double precision
                   !Cutoff: A range of interaction real double precision 
                   !VDW_force: Total interaction force                       

                Double precision,allocatable, dimension(:,:) :: positions
                Double precision :: cutoff,boxsize
                Double precision,allocatable,dimension(:,:) :: vdw_force

                 !VARIABLES:
                    !r_ij : relative position vector between pair of particles, double precision dim= (3,1)
                    !d_ij : distance module between the particles, double precision 
                    !f_ij: Force between the particles
                    !npart: particle number. Integer

                Double precision, dimension(3,1) :: r_ij
                Double precision, dimension(3,1) :: f_ij
                Double precision :: d_ij
                Integer ::  npart,i,j

                call PBC(positions, boxsize, npart)
                vdw_force=0.0
                npart= int(size(positions,dim=2))

                do i=1,npart-1
                do j=i+1,npart
                        !Distance between particles:
                        r_ij(1,1)=positions(1,i)-positions(1,j)
                        r_ij(2,1)=positions(2,i)-positions(2,j)
                        r_ij(1,1)=positions(3,i)-positions(3,j)

                        !Module of r_ij

                        d_ij=((r_ij(1,1)**2.d0)+(r_ij(2,1)**2.d0)+(r_ij(3,1)**2.d0))**(1.d0/2.d0)
                        !Now we compare this distance with the cutoff

                        if (d_ij< cutoff) then
                                 !Force made by j to i
                                 vdw_force(1,i) = vdw_force(1,i) + (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(1,1)
                                 vdw_force(2,i) = vdw_force(2,i) + (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(2,1)
                                 vdw_force(3,i) = vdw_force(3,i) + (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(3,1)

                                 !Force made by i to j
                                 vdw_force(1,j) = vdw_force(1,j) - (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(1,1)
                                 vdw_force(2,j) = vdw_force(2,j) - (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(2,1)
                                 vdw_force(3,j) = vdw_force(3,j) - (48.d0 / d_ij**14.d0 - 24.d0 / d_ij**8.d0) * r_ij(3,1)

                         end if
                   end do
          end do

          end subroutine VDW_forces

             


          !!!!!!!!!!!!-------------------------------------------------------------------------------



          !!!!!!!!!!!         ENERGIES   !!!!!!!!!!!!!!

          !!!!!!!!!!!-------------------------------------------------------------------------------

          subroutine potentialE (positions,cutoff,boxsize,PotentialEn)

                  !Subroutine to calculate the potential energy between particles using Lennard-Jones potential

                implicit none
               

                !ARGUMENTS:
                   !Positions : Positions of the particles DIM= (d,npart)  (d usually=3)
                   !Boxsize: We will supose cubic system, double precision
                   !Cutoff: A range of interaction real double precision
                   !PotentialE: Total potential energy

                Double precision,allocatable, dimension(:,:) :: positions
                Double precision :: cutoff,boxsize
                Double precision :: PotentialEn

                 !VARIABLES:
                    !r_ij : relative position vector between pair of particles, double precision dim= (3,1)
                    !d_ij : distance module between the particles, double precision
                    !f_ij: Force between the particles
                    !npart: particle number. Integer

                Double precision, dimension(3,1) :: r_ij
                Double precision:: e_ij
                Double precision :: d_ij
                Integer ::  npart,i,j

                potentialEn=0.0
                npart= int(size(positions,dim=2))

                call PBC(positions, boxsize, npart)
                do i=1,npart-1
                do j=i+1,npart
                        !Distance between particles:
                        r_ij(1,1)=positions(1,i)-positions(1,j)
                        r_ij(2,1)=positions(2,i)-positions(2,j)
                        r_ij(1,1)=positions(3,i)-positions(3,j)

                       
                        !THIS WAY WE MAKE SURE THAT r_ij is inside our box

                        !Module of r_ij

                        d_ij=((r_ij(1,1)**2.d0)+(r_ij(2,1)**2.d0)+(r_ij(3,1)**2.d0))**(1.d0/2.d0)
                        !Now we compare this distance with the cutoff

                        if (d_ij< cutoff) then
                                 e_ij = 4*(1.d0 / d_ij**12.d0 - 1.d0 / d_ij**6.d0) 
                                

                        
                                potentialEn = potentialEn + e_ij
                         end if
                   end do
          end do
        end subroutine potentialE
        


        !!!!!!!!-------------------------------------------------------------------------------------------



        subroutine  kineticE (velocities,KineticEn)

                !Subroutine to calculate the kinetic energy of the particles

                implicit none

                !ARGUMENTS:
                   !Velocities : Velocities of the particles DIM= (d,npart)  (d usually=3)
                   !KineticE: Total kinetic energy

                Double precision,allocatable, dimension(:,:) :: velocities
                Double precision :: KineticEn

                 !VARIABLES:
                    !r_ij : relative position vector between pair of particles, double precision dim= (3,1)
                    !d_ij : distance module between the particles, double precision
                    !f_ij: Force between the particles
                    !npart: particle number. Integer

                Double precision, dimension(3,1) :: v_ij
                Double precision:: e_ij
                Double precision :: vel_ij
                Integer ::  npart,i
                
                npart= int(size(velocities,dim=2))
                
                KineticEn=0.0

                do i=1,npart

                        vel_ij=((v_ij(1,1)**2.d0)+(v_ij(2,1)**2.d0)+(v_ij(3,1)**2.d0))**(1.d0/2.d0)
                        
                        e_ij= (1.d0/2.d0)*vel_ij**2.d0

                        KineticEn =KineticEn + e_ij
                enddo
        end subroutine kineticE


end module Forces_and_Energies

