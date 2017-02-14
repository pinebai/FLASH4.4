!!****if* source/monitors/Logfile/LogfileMain/log_allocateStrArr
!!
!! NAME
!!
!!  log_allocateStrArr
!!
!! SYNOPSIS
!!
!!  log_allocateStrArr(integer, intent(in)  :: length,
!!                     integer, intent(in)  :: dim)
!!
!! DESCRIPTION
!!  Allocates a 2 dimensional array and initializes each entry to an empty string
!!
!! ARGUMENTS
!!
!!   length : size of first dimension of string to allocate
!!
!!   dim : size of second dimension of log string to allocate
!!
!!
!!
!!***

subroutine log_allocateStrArr(length, dim)

  use Logfile_data, ONLY : log_strArr
  
implicit none
  integer, intent(in) :: length, dim
  integer :: i, j
  
  allocate (log_strArr(length, dim))
  
  do i = 1, length
     do j = 1, dim
        log_strArr(i,j) = ''
     end do
  end do

end subroutine log_allocateStrArr
