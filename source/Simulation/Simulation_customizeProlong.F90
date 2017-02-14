!!****f* source/Simulation/Simulation_customizeProlong
!!
!! NAME
!!
!!  Simulation_customizeProlong
!!
!! SYNOPSIS
!!
!!  call Simulation_customizeProlong(integer(IN)  :: beforeOrAfter)
!!
!! DESCRIPTION
!!
!!  The Simulation_customizeProlong interface provides a way to
!!  customize the prolongation of Grid data that normally happens
!!  after an AMR Grid has changed - in particular, the interpolation
!!  of data into blocks that were newly created by refining existing
!!  blocks.
!!
!!  After the refinement, apply the user-defined routine for a 
!!  special prolongation routine from the coarse parents to 
!!  the newly created child blocks.
!!
!!  The interface is called twice for each time that the global
!!  prolongation operation is applied: once just before prolongation
!!  gets applied, and then again after prolongation is done.
!!  The single argument beforeOrAfter is used to distinguish the calls.
!!
!!  Generally, this routine is only present as a stub and does
!!  nothing.
!!
!! ARGUMENTS
!!
!!   beforeOrAfter : BEFORE when called before prolongation,
!!                   AFTER when called after prolongation.
!!
!!
!! EXAMPLE  
!!
!!  This is how this interface is normally called, in one of the
!!  subroutines that implement Grid_updateRefinement processing:
!!
!!    #include "constants.h" 
!!    #include "Flash.h" 
!!    ...
!!       use Grid_data,ONLY: gr_meshMe            ! my rank
!!    ...
!!       ! call before to modify prolongation method:
!!       call Simulation_customizeProlong(BEFORE)
!!       ! PARAMESH routine that does prolongation:
!!       call amr_prolong (gr_meshMe, 1, NGUARD)
!!       ! call after to restore prolongation method:
!!       call Simulation_customizeProlong(AFTER)
!!    ...
!!
!! NOTES
!!
!!   The constants BEFORE and AFTER are defined in constants.h.
!!
!!   As of FLASH 3.1.1, a non-stub implementation for use by
!!   MHD simulations is provided in the directory tree location
!!   source/Simulation/SimulationMain/magnetoHD . All simulations
!!   placed under the magnetoHD directory will therefore this
!!   implementation by default.  This is usually desired for
!!   simulations using MHD.  For this reason, simulations using
!!   MHD should be placed under the magnetoHD directory.
!!
!!   See documentation headers in the magnetoHD implementation for
!!   details about that implementation.
!!***

subroutine Simulation_customizeProlong(beforeOrAfter)

  implicit none

  integer, intent(IN) :: beforeOrAfter ! BEFORE for before, AFTER for after
     

end subroutine Simulation_customizeProlong
