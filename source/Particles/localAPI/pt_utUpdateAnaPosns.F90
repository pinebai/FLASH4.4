!!****if* source/Particles/localAPI/pt_utUpdateAnaPosns
!!
!! NAME
!!
!!  pt_utUpdateAnaPosns
!!
!! SYNOPSIS
!!
!!  pt_utUpdateAnaPosns(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!
!!  Updates particles' POS{X,Y,Z}_PART_PROP and VEL{X,Y,Z}_PART_PROP
!!  properties.
!!
!!  This version just sets the new coordinates and velocities
!!  based on an analytically known solution.
!!
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!   t     -- time for which solution is sought
!!  
!!***

!===============================================================================

subroutine pt_utUpdateAnaPosns (dtOld,dtNew,t)
    
  
  implicit none

  
  real, INTENT(in)  :: dtOld, dtNew, t
     
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_utUpdateAnaPosns


