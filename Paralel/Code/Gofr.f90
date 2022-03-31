program gofr_3D

    implicit none
    ! inputs
    character (len=80) :: infile, datafile, outfile
    integer :: boxes, nconfig, nparts
    double precision :: length, rho
    double precision, allocatable :: positions(:,:)
    ! outputs
    double precision, allocatable :: gofr(:)
    double precision :: r, r_mostprob
    ! others
    character :: dummy
    integer :: i, j
    double precision :: dr
    
    infile = "gofr_params.txt"

    ! read the parameters
    open(10,file=infile)
    read(10,*) nconfig
    read(10,*) nparts
    read(10,*) length
    read(10,*) rho
    read(10,*) boxes
    read(10,*) datafile
    read(10,*) outfile
    close(10)


    ! determine the dimension of the postioins array and histogram
    allocate(positions(nparts,3), gofr(boxes))

    ! set all boxes to zero
    gofr = 0.d0
    
    open(10,file=datafile)
    do i = 1, nconfig
        positions = 0.d0
        read(10,*)
        read(10,*)
        do j = 1, nparts
            ! read the positions of the ith configuration
            read(10,*) dummy, positions(j,1), positions(j,2), positions(j,3)
        enddo
        ! sample the distances between the npart particles
        call sample_distances(positions,gofr,nparts,boxes,length)
    enddo
    close(10)

    ! normalize the histogram
    call normalization(gofr,dr,nconfig,boxes,length,rho)

    ! calculate the most probable distance
    call most_probable(gofr,boxes,dr,r_mostprob)

    ! write the results
    open(10,file=outfile)
    write(10,*) "# maximum at r(A) = ", r_mostprob
    write(10,*) "# r, g(r)"
    do i = 1, boxes
        r = dr*i - dr*0.5d0
        write(10,*) r, gofr(i)
    enddo
    close(10)
    
    print*, "Radial distribution function was calculated."
    
end program gofr_3D

!**********************************************************************************************************************************!
!                                                        Subroutines                                                               !
!**********************************************************************************************************************************!
!
! periodic boundary conditions
subroutine pbc(dr,length)

    implicit none
    double precision :: dr(3)
    double precision :: length
    integer :: k

    do k = 1, 3
        if (dr(k) > +0.5d0*length) dr(k) = dr(k) - length
        if (dr(k) < -0.5d0*length) dr(k) = dr(k) + length
    enddo

    return

end subroutine pbc
!
! distance between 2 particles
double precision function distance(x1,x2,l)

    implicit none
    double precision, intent(in) :: x1(3), x2(3), l
    double precision, external   :: pbc_func
    double precision             :: x12(3)

    x12 = x1-x2
    call pbc(x12,l)
    distance = dsqrt(sum(x12**2))

    return

end function distance
!
! sample distances
subroutine sample_distances(positions,gofr,nparts,boxes,length)

    implicit none
    integer, intent(in) :: nparts, boxes
    double precision, intent(in) :: positions(nparts,3), length
    double precision, intent(inout) :: gofr(boxes)
    double precision, external :: distance
    integer :: i, j, xi
    double precision :: x1(3), x2(3), xr, rmax, dr

    ! distances greater than L/2 are neglected
    rmax = 0.5d0*length
    dr = rmax/dble(boxes)

    do i = 1, nparts-1
        do j = i+1, nparts
            ! calculate the distance between particles i and j with pbc
            x1 = positions(i,:)
            x2 = positions(j,:)
            xr = distance(x1,x2,length)
            ! if the distance is less than L/2 add 1 in the corresponding box
            if (xr < rmax) then
                xi = int(xr/dr) + 1
                gofr(xi) = gofr(xi) + 1
            endif
        enddo
    enddo

    return

end subroutine sample_distances
!
! normalization: sum(r^2*g(r)*dr) = (N-1)/(4pi*rho)
subroutine normalization(gofr,dr,nconfig,boxes,length,rho)

    implicit none
    integer, intent(in) :: nconfig, boxes
    double precision, intent(in) :: length, rho
    double precision, intent(inout) :: gofr(boxes)
    double precision, intent(out) :: dr
    double precision, parameter :: pi = 4.d0*datan(1.d0)
    integer :: i
    double precision :: r

    dr = 0.5d0*length/dble(boxes)

    ! apply the normalization for gofr in 3D
    do i = 1, boxes
        r = dr*i - dr*0.5d0
        gofr(i) = gofr(i)/(4.d0*pi*r**2*dr*rho*dble(nconfig))
    enddo
    ! ad hoc correction: due to the aproximations g(r) doesn't go exactly to 1 for large r (but usually very close). Here we force
    ! the histogram to do it.
    gofr = gofr/gofr(boxes)
    return

end subroutine normalization
!
! most probable distance
subroutine most_probable(gofr,boxes,dr,r)

    implicit none
    integer, intent(in) :: boxes
    double precision, intent(inout) :: gofr(boxes), dr
    double precision, intent(out) :: r
    integer :: i
    double precision :: maxprob, ri

    ! save the position of the box with a greater contribution
    maxprob = 0.d0
    do i = 1, boxes
        ri = dr*i - dr*0.5d0
        if (gofr(i) > maxprob) then
            maxprob = gofr(i)
            r = ri
        endif
    enddo

    return

end subroutine most_probable
!
!**********************************************************************************************************************************!
!                                                                                                                                  !
!**********************************************************************************************************************************!
