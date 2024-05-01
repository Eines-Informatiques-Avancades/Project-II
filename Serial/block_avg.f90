module block_average_module
    implicit none
    public :: block_average, read_data_from_file, write_results, compute_and_save_block_averages
contains

    subroutine block_average(array, num_data, block_size, total_mean, std_dev)
        implicit none
        integer :: i, j, num_blocks
        integer, intent(in) :: num_data, block_size
        real*8, intent(in) :: array(:)
        real*8, intent(out) :: total_mean, std_dev
        real*8 :: block_avg, block_avg2, total_mean2

        num_blocks = num_data / block_size
        total_mean = 0.0
        total_mean2 = 0.0

        do i = 1, num_blocks
            block_avg = 0.0
            block_avg2 = 0.0

            do j = 1, block_size
                block_avg = block_avg + array((i-1) * block_size + j)
                block_avg2 = block_avg2 + array((i-1) * block_size + j)**2.d0
            end do

            block_avg = block_avg / block_size
            block_avg2 = block_avg2 / block_size
            total_mean = total_mean + block_avg
            total_mean2 = total_mean2 + block_avg2
        end do

        total_mean = total_mean / num_blocks
        total_mean2 = total_mean2 / num_blocks
        std_dev = sqrt(total_mean2 - total_mean**2)
        std_dev = std_dev / sqrt(real(num_blocks, kind=8))
    end subroutine block_average

    subroutine read_data_from_file(filename, data, num_data, column)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: column
        real*8, allocatable, intent(out) :: data(:)
        integer, intent(out) :: num_data
        integer :: i, j, ios
        real*8 :: read_data, discard

        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Failed to open file"
            stop
        end if

        num_data = 0
        do
            read(10, *, iostat=ios) (discard, j=1, column-1), read_data
            if (ios /= 0) exit
            num_data = num_data + 1
        end do
        rewind(10)

        allocate(data(num_data))
        do i = 1, num_data
            read(10, *) (discard, j=1, column-1), data(i)
        end do

        close(10)
    end subroutine read_data_from_file

    subroutine write_results(output_file, total_mean, std_dev)
        implicit none
        character(len=*), intent(in) :: output_file
        real*8, intent(in) :: total_mean, std_dev

        open(unit=20, file=output_file, status='replace', action='write')
        write(20, '(F10.5, 2X, F10.5)') total_mean, std_dev
        close(20)
    end subroutine write_results

    subroutine compute_and_save_block_averages(filename, block_size, column, output_file)
        implicit none
        real*8, allocatable :: data(:)
        character(len=*), intent(in) :: filename, output_file
        integer, intent(in) :: block_size, column
        integer :: num_data
        real*8 :: total_mean, std_dev

        call read_data_from_file(filename, data, num_data, column)
        call block_average(data, num_data, block_size, total_mean, std_dev)
        call write_results(output_file, total_mean, std_dev)

        deallocate(data)
    end subroutine compute_and_save_block_averages

end module block_average_module
