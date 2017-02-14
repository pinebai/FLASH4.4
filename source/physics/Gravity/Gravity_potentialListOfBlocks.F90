!!****f* source/physics/Gravity/Gravity_potentialListOfBlocks
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
!!      This routine computes the gravitational potential on all
!!      blocks specified in the list, for the gravity implementations
!!      (i.e., various Poisson implementations), which make use of it
!!      in computing the gravitational acceleration.
!!
!!      The resulting potential can be considered as zone-averaged.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.  For some implementation of Gravity,
!!      in particular with Barnes-Hut tee solver, additional combinations
!!      of boundary conditions may be supported.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!   potentialIndex : If present, determines which variable in UNK to use
!!                    for storing the updated potential.  If not present,
!!                    GPOT_VAR is assumed.
!!                    Presence or absense of this optional dummy argument
!!                    also determines whether some side effects are enabled
!!                    or disabled; see discussion of two modes under NOTES
!!                    below.
!!
!! NOTES
!!
!!  Gravity_potentialListOfBlocks can operate in one of two modes:
!!  * automatic mode  - when called without the optional potentialIndex.
!!    Such a call will usually be made once per time step, usually
!!    from the main time advancement loop in Driver_evolveFlash.
!!    Various side effects are enabled in this mode, see SIDE EFFECT below.
!!
!!  * explicit mode  - when called with the optional potentialIndex.
!!    The potential is stored in the variable explicitly given, and
!!    side effects like saving the previous potential in GPOL_VAR
!!    and updating some sink particle state and properties are
!!    suppressed.
!!
!!  This call does nothing for those Gravity implementations provided with
!!  FLASH that just provide a given (time-independent) potential for other
!!  physics modules and don't compute the potential by invoking a solver:
!!  Constant, PointMass, and PlanePar.
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
!!  There are additional side effects if sink particles are used.
!!  These effects happen by calls to Particles_sinkAccelGasOnSinks and
!!  Particles_sinkAccelSinksOnGas, which may update sink particle
!!  properties and additional UNK variables that store
!!  accelerations. The calls are only made in automatic mode.
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
!!***

subroutine Gravity_potentialListOfBlocks(blockCount,blockList, potentialIndex)

!=============================================================================
  implicit none


  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN), optional :: potentialIndex

  return
end subroutine Gravity_potentialListOfBlocks
