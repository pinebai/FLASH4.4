!!****if* source/flashUtilities/Pipeline/localAPI/pl_printGlobalStatusVector
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

  implicit none

  integer, intent (in) :: fileUnit

  return
end subroutine pl_printGlobalStatusVector
