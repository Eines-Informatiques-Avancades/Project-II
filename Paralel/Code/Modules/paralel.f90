module paralel
 include 'mpif.h'
! Open MPI variables and distribution of particles between workers 
integer comm, workerid, numproc, ierror, partner,request!Variable del mpi
integer reslen
integer message
integer stat(MPI_STATUS_SIZE)
integer, dimension(:,:), allocatable :: index_matrix,index_matrix2,pairindex!, double_matrix
integer, dimension(:), allocatable :: desplac,num_send
  

integer,parameter    :: master=0 !PROCESADOR JEFE
integer              :: mida, N
integer              :: first_part, last_part
! required variables for MPI functions
integer, allocatable :: sizes_part(:), displs_part(:)


integer  :: Npart, chunksize, remainder, counter, i,ii,j,b,b_res !VARIABLES EMPLEADAS EN EL PROGRAMA, NPART SE LEE DEL MAIN
contains
! Definition of particles distribution
subroutine paralel_particle_distribution(Npart)
!-------------------------------------------------------------!
!-------------------------------------------------------------!
!-----FUNCIONES QUE SE NECESITAN PARA PARALELIZAR BASICAS-----!
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, workerid, ierror)!Worker id: cada procesador
call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)!Numproc numero total de proc
!-------------------------------------------------------------!
!-------------------------------------------------------------!
!-------------------------------------------------------------!


    
chunksize = int(Npart/real(numproc)) ! number of particles for each worker for "clean" distribution. 
remainder = mod(Npart,numproc) ! number of workers which will have (chunksize - 1) particles. Cuantas particulas quedan sin procesador
!-------------------------------------------------------------!
!---------ASIGNAMOS TAMAÑO A LAS MATRICES NECESARIAS----------!
!-------------------------------------------------------------!
allocate(index_matrix(numproc,2),desplac(numproc),num_send(numproc)) !Index matrix contiene en numproc,1 la primera particula y en numproc 2 la ultima
allocate(index_matrix2(numproc,2))
allocate(pairindex((n_particles*(n_particles-1))/2,2))
!-------------------------------------------------------------!
!-------------------------------------------------------------!
!-------------------------------------------------------------!

!-------------------------------------------------------------!
!-----------Numero de particulas para cada trabajador---------!
!-------------------------------------------------------------!
   if (workerid<remainder) then !si el identificador del trabajador ess más pequño que el numero de particulas sin trabajador
      index_matrix(i,1) = workerid*(chunksize+1)+1!entonces la primera particula (por empezar en 0 el id se suma el 1) 
      !chunksize+1 es para empezar en la ultima particula del anterior
      index_matrix(i,2) = index_matrix(i,1) + chunksize !la utima es la primera mas el numero q le corresponda

    else!si el identificador es mayor que el numero de particulas que quedan sin procesador
      index_matrix(i,1)= workerid*chunksize+remainder+1  !ramainder siempre será inferior al numero total de procesadores
      !la primera particula es el trabajador por el tamaño de particulas mass el resto 
      index_matrix(i,2) = index_matrix(i,1) + chunksize -1
    endif
!-------------------------------------------------------------!
!-------------------------------------------------------------!
!-------------------------------------------------------------!
 
!-------------------------------------------------------------!   
!-----------Compute de relative displacement------------------! 
!-------------------------------------------------------------!
counter=0
do i=1,numproc
	if(i.eq.1)then
		desplac(i)=0
	else
		desplac(i)=index_matrix(i-1,2)
	end if
		num_send(i)=index_matrix(i,2)-index_matrix(i,1)+1
end do
!-------------------------------------------------------------!
!-------------------------------------------------------------!
!-------------------------------------------------------------!

!-------------------------------------------------------------!
!--------Definimos vector para los pares de particulas--------!
!-------------------------------------------------------------!

ii = 1
do i = 1,Npart-1
    do j = i+1,Npart
        pairindex(ii,:) = (/i,j/)
        ii = ii + 1
    enddo
enddo
b =  INT(REAL((Npart*(Npart-1))/2)/REAL(numproc))
b_res=mod((Npart*(Npart-1))/2,numproc)
!----Computing the number of pairs for each worker------------!
do i=1,numproc
	if ((i-1)<b_res)then
		index_matrix2(i,1)=(b+1)*(i-1)+1
		index_matrix2(i,2)=index_matrix2(i,1)+b
	else
		index_matrix2(i,1)=(i-1)*b+b_res+1
		index_matrix2(i,2)=index_matrix2(i,1)+b-1
	end if
end do

!-------------------------------------------------------------!
!-------------------------------------------------------------!
!-------------------------------------------------------------!

 end subroutine paralel_particle_distribution
end module paralel
