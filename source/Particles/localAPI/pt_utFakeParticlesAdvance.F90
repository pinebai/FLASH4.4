!!****if* source/Particles/localAPI/pt_utFakeParticlesAdvance
!!
!! NAME
!!
!!  pt_utFakeParticlesAdvance
!!
!! SYNOPSIS
!!
!!  pt_utFakeParticlesAdvance(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t,
!!                           integer(in) :: ind)
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
!!   ind   -- index into pt_typeInfo
!!  
!! PARAMETERS
!!
!!
!!***

!===============================================================================

subroutine pt_utFakeParticlesAdvance (dtOld,dtNew,t, ind)
    
  implicit none

  real, INTENT(in)  :: dtOld, dtNew, t
  integer, intent(in) :: ind
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_utFakeParticlesAdvance


