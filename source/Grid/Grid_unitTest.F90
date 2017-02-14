!!****f* source/Grid/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  Grid_unitTest(integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  The Grid unit test has several implementations. There is a general test,
!!  which exercises the data access interface of the Grid unit, for example
!!  the get/putData functions, and can be used with all GridMain
!!  implementations. There is also a test specific to
!!  the Uniform Grid which tests internal working of the 
!!  UG such as parallel data exchanges, domain setup, etc. 
!!   
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!! NOTES
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  implicit none

  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  perfect = .false.

  return
 
end subroutine Grid_unitTest
