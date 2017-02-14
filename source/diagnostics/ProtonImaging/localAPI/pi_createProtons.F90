!!****if* source/diagnostics/ProtonImaging/localAPI/pi_createProtons
!!
!! NAME
!!
!!  pi_createProtons
!!
!! SYNOPSIS
!!
!!  call pi_createProtons (integer, intent (in) :: blockCount, 
!!                         integer, intent (in) :: blockList (:),
!!                         real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates protons and places them in their initial blocks. This routine calls the
!!  appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeSimulation : current simulation time
!!
!!***

subroutine pi_createProtons (blockCount, blockList, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeSimulation

  return
end subroutine pi_createProtons
