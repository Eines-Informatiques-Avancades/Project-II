module paralel
 include 'mpif.h'
 
! Open MPI variables and distribution of particles between workers

  integer              :: ierror, workerid, numproc
  integer,parameter    :: master=0
  integer              :: mida, N
  integer              :: first_part, last_part
! required variables for MPI functions
  integer, allocatable :: sizes_part(:), displs_part(:)

  contains
! Definition of particles distribution
subroutine paralel_particle_distribution(Npart)
integer  :: Npart, chunksize, remainder, counter, i
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, workerid, ierror)!Worker id: cada procesador
call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)!Numproc numero total de proc
allocate( sizes_part(numproc))
    allocate(displs_part(numproc))
    chunksize = Npart/real(numproc) ! number of particles for each worker for "clean" distribution
    remainder = mod(Npart,numproc) ! number of workers which will have (chunksize - 1) particles

   if (workerid<remainder) then
      first_part = workerid*(chunksize+1)+1
      last_part = first_part + chunksize

    else
      first_part = workerid*chunksize+remainder+1
      last_part = first_part + chunksize -1
    endif
print*,'CS',chunksize ! sizes_part contains the number of elements that are received/sent from each process.
    counter=0
    do i=1,numproc
     displs_part(i) = counter
      if ((i-1)<remainder) then
       sizes_part(i)= chunksize+1
     else
       sizes_part(i)= chunksize
      end if
     counter = counter + sizes_part(i)
    end do
print*,'Worker ID',workerid,'Primera : ',first_part,'U: ',last_part 
print*,'Procesador',numproc
  end subroutine paralel_particle_distribution
end module paralel
