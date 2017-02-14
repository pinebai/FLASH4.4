!!****if* source/physics/Gravity/GravityMain/UserDefined/Gravity_potentialListOfBlocks
!!
!! NAME
!! 
!!  Gravity_potentialListOfBlocks 
!!
!! SYNOPSIS
!!
!!  call Gravity_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                     integer(IN) :: blockList(blockCount),
!!                            optional,integer(IN) :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  Compute the zone-averaged gravitational potential on all blocks
!!  specified in the list.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!   potentialIndex : If present, determines which variable in UNK to use
!!                    for storing the updated potential.  If not present,
!!                    GPOT_VAR is assumed.
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  gravitational potential.  Invokes a solver (of the Poisson equation)
!!  if necessary. On return, if potentialIndex is not present,
!!     GPOT_VAR:  contains potential for the current simulation time.
!!     GPOL_VAR (if defined): contains potential at the previous simulation time.
!!  On return, if potentialIndex is present, the UNK variable given by
!!  potentialIndex contains the newly computed potential.
!!
!!  May affect other variables related to particle properties if particles
!!  are included in the simulation.  In particular,
!!     PDEN_VAR (if defined): may get updated to the current density from
!!                particles if particles have mass.
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Gravity implementation.
!!  The following information is subject to change without notice.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!     IMGM_VAR (image mass)
!!     IMGP_VAR (image potential)
!!  For the Multipole implementation:
!!     (none)
!!
!! NOTES
!!
!!  This implementation is a stub to be replaced by the end-user
!!  within their Simulation directory.
!!***

subroutine Gravity_potentialListOfBlocks(blockCount,blockList, potentialIndex)

!=============================================================================
  implicit none


  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN), optional :: potentialIndex

  return
end subroutine Gravity_potentialListOfBlocks
