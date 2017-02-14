!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_printGlobalStatusVector
!!
!! NAME
!!
!!  pl_printGlobalStatusVector
!!
!! SYNOPSIS
!!
!!  call pl_printGlobalStatusVector (integer, intent (in) :: fileUnit)
!!
!! DESCRIPTION
!!
!!  Prints the global status vector to the file associated with the passed file unit number.
!!
!! ARGUMENTS
!!
!!  fileUnit : the file unit number
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine pl_printGlobalStatusVector (fileUnit)

  use Pipeline_data, ONLY : pl_procStatusGlobal, &
                            pl_size

#include "Pipeline.h"

  implicit none

  integer, intent (in) :: fileUnit

  integer :: rank
!
!
!     ...Print out the vector.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') ' Pipeline Global Status Vector '
  write (fileUnit,'(/)')
  write (fileUnit,'(/,3x,4(a))') ' Rank ',' Comm Status ',' Recv Status ',' Item Status '
  write (fileUnit,'(/)')

  do rank = 0, pl_size-1
     write (fileUnit,'(i6,3(5x,i3,5x))') rank, pl_procStatusGlobal (rank,PL_STATUS_COMM), &
                                               pl_procStatusGlobal (rank,PL_STATUS_RECV), &
                                               pl_procStatusGlobal (rank,PL_STATUS_ITEM)
  end do

  write (fileUnit,'(/)')
!
!
!     ...Ready!
!
!
  return
end subroutine pl_printGlobalStatusVector
