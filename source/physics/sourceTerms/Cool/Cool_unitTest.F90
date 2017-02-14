!!****f* source/physics/sourceTerms/Cool/Cool_unitTest
!!
!! NAME
!!
!!  Cool_unitTest
!!
!! SYNOPSIS
!!
!!  Cool_unitTest(
!!                integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This function will implement a unit test for the cooling source term
!!   
!!
!! ARGUMENTS
!!  
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!! NOTES
!!
!!***
subroutine Cool_unitTest (  fileUnit, perfect )

  implicit none

  integer, intent(IN) ::  fileUnit
  logical, intent(INOUT) :: perfect

  return

end subroutine Cool_unitTest
