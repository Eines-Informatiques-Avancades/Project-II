module block_average
    implicit none

    contains

    subroutine block_average(data, block_size, averages, variances, num_blocks)
        implicit none
        ! Input data and parameters
        double precision, intent(in) :: data(:)        ! Array of data points to be averaged
        integer, intent(in) :: block_size              ! Number of data points per block

        ! Output data
        double precision, allocatable, intent(out) :: averages(:), variances(:)
        integer, intent(out) :: num_blocks

        ! Local variables
        integer :: i, j
        double precision :: sum_block, sum_sq_block, mean, variance

        ! Calculate the number of complete blocks
        num_blocks = size(data) / block_size
        if (num_blocks * block_size < size(data)) then
            num_blocks = num_blocks + 1
        end if

        ! Allocate memory for output arrays
        allocate(averages(num_blocks))
        allocate(variances(num_blocks))

        ! Compute averages and variances for each block
        do i = 1, num_blocks
            sum_block = 0.0
            sum_sq_block = 0.0
            do j = 1, block_size
                if ((i - 1) * block_size + j <= size(data)) then
                    sum_block = sum_block + data((i - 1) * block_size + j)
                    sum_sq_block = sum_sq_block + data((i - 1) * block_size + j)**2
                end if
            end do
            mean = sum_block / min(block_size, size(data) - (i - 1) * block_size)
            variance = (sum_sq_block / min(block_size, size(data) - (i - 1) * block_size)) - mean**2
            averages(i) = mean
            variances(i) = variance
        end do
    end subroutine block_average

end module block_average