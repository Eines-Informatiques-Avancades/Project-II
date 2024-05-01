module verlet
    implicit none
    public :: verletlist

contains

    subroutine verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
        !!! --- Authors: Paula Sierra and Emma Valdés --- !!!
        implicit none
        integer, intent(in) :: imin, imax, N
        integer, intent(out) :: nnlist(:), vlist(:) !(# neighbors x part), (i_neighbor)
        real*8, dimension(:,:), intent(in) :: positions
        real*8, intent(in) :: vcutoff
        integer :: i,j,k
        real*8 :: rij, vcutoff2
        vcutoff2=vcutoff*vcutoff
        k=1
        vlist = 0
        ! i, j are real particle labels
        ! nnlist is indexed 1 - Nsub regardless of the worker rank
        do i=imin,imax
            nnlist(i-imin+1) = 0
            do j=1,N
                if (i.NE.j) then
                    rij=(positions(i,1)-positions(j,1))**2+(positions(i,2)-positions(j,2))**2+(positions(i,3)-positions(j,3))**2
                    if (rij<vcutoff2) then
                        nnlist(i-imin+1)=nnlist(i-imin+1)+1
                        vlist(k)=j
                        k=k+1
                    endif
                endif
            enddo
        enddo

    end subroutine verletlist

end module verlet
