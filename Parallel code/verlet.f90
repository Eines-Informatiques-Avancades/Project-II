module verlet
    implicit none
    public :: verletlist

contains

    subroutine verletlist(imin,imax,N,positions,vcutoff,nnlist,vlist)
        implicit none
        integer, intent(in) :: imin, imax, N
        integer, intent(out) :: nnlist(imin:imax), vlist(:) !(# neighbors x part), (i_neighbor)
        real*8, dimension(:,:), intent(in) :: positions
        real*8, intent(in) :: vcutoff
        integer :: i,j,k
        real*8 :: rij, vcutoff2

        vcutoff2=vcutoff*vcutoff

        do i=imin,imax
            nnlist(i)=0
        enddo

        k=1
        do i=imin,imax
            do j=1,N
                rij=(positions(i,1)-positions(j,1))**2+(positions(i,2)-positions(j,2))**2+(positions(i,3)-positions(j,3))**2
                if (rij<vcutoff2) then
                    nnlist(i)=nnlist(i)+1
                    vlist(k)=j
                    k=k+1
                endif
            enddo
        enddo

    end subroutine verletlist
end module verlet
