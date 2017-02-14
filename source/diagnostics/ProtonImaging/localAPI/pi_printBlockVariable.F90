!!****if* source/diagnostics/ProtonImaging/localAPI/pi_printBlockVariable
!!
!! NAME
!!
!!  pi_printBlockVariable
!!
!! SYNOPSIS
!!
!!  call pi_printBlockVariable (integer, intent (in) :: blockID,
!!                              integer, intent (in) :: variable,
!!                              integer, intent (in) :: fileUnit)
!!
!! DESCRIPTION
!!
!!  Prints a certain cell variable content of a block on the current processor to the file
!!  associated with the unit number 'fileUnit'. To help locating the block, its bounding
!!  Box is also printed. The whole block is printed, excluding the guard cells.
!!
!! ARGUMENTS
!!
!!  blockID  : the block number ID on the current processor
!!  variable : the integer handle for the variable to be printed
!!  fileUnit : the unit number for the printout file
!!
!! NOTES
!!
!!  Only certain variables (those relevant for proton imaging) can be printed right now.
!!  Any attempt to print a different variable will result in an informal message, but the
!!  calculation will continue. The variables that can be printed are:
!!
!!               1) magnetic x-component (handle MAGX_VAR)
!!               2) magnetic y-component (handle MAGY_VAR)
!!               3) magnetic z-component (handle MAGZ_VAR)
!!               4) electric x-component (handle ELEX_VAR)
!!               5) electric y-component (handle ELEY_VAR)
!!               6) electric z-component (handle ELEZ_VAR)
!!
!!***

subroutine pi_printBlockVariable (blockID, variable, fileUnit)

  implicit none

  integer, intent (in) :: blockID
  integer, intent (in) :: variable
  integer, intent (in) :: fileUnit

  return
end subroutine pi_printBlockVariable
