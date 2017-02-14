!!****f* source/Grid/Grid_solidBodyUnitTest
!!
!! NAME
!!
!!  Grid_solidBodyUnitTest
!!
!! SYNOPSIS
!!
!!  Grid_solidBodyUnitTest(integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  Stub
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***


subroutine Grid_solidBodyUnitTest(fileUnit,perfect)
  implicit none
  integer, intent(IN)           :: fileUnit ! Output to file
  logical, intent(INOUT)        :: perfect  ! Flag to indicate errors
end subroutine Grid_solidBodyUnitTest
