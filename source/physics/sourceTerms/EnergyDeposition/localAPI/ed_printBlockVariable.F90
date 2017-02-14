!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_printBlockVariable
!!
!! NAME
!!
!!  ed_printBlockVariable
!!
!! SYNOPSIS
!!
!!  call ed_printBlockVariable (integer, intent (in) :: blockID,
!!                              integer, intent (in) :: variable,
!!                              integer, intent (in) :: fileUnit)
!!
!! DESCRIPTION
!!
!!  Prints a certain cell variable content of a block on the current processor to the file
!!  associated with the unit number 'fileUnit'. To help locating the block, its bounding
!!  Box is also printed. The whole block is printed, including the guard cells.
!!
!! ARGUMENTS
!!
!!  blockID  : the block number ID on the current processor
!!  variable : the integer handle for the variable to be printed
!!  fileUnit : the unit number for the printout file
!!
!! NOTES
!!
!!  Only certain variables (those relevant for the laser deposition) can be printed
!!  right now. Any attempt to print a different variable will result in an informal
!!  message, but the calculation will continue. The variables that can be printed
!!  are:
!!
!!               1) energy deposition    (handle DEPO_VAR)
!!               2) mass density         (handle DENS_VAR)
!!               3) electron temperature (handle TELE_VAR)
!!
!!***

subroutine ed_printBlockVariable (blockID, variable, fileUnit)

  implicit none

  integer, intent (in) :: blockID
  integer, intent (in) :: variable
  integer, intent (in) :: fileUnit

  return
end subroutine ed_printBlockVariable
