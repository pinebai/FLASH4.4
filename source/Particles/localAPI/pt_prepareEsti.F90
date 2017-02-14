!!****if* source/Particles/localAPI/pt_prepareEsti
!!
!! NAME
!!
!!  pt_prepareEsti
!!
!! SYNOPSIS
!!
!!  pt_prepareEsti(real(in) :: dtOld,
!!                    real(in) :: dtNew,
!!                    real(inout):: particles(:,p_count),
!!                    integer(in):: p_count,
!!                    integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the passive particle module.
!!  This portion cleans up and prepares for the next time step.
!!  
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time interval
!!   dtNew -- current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind  -- index of this particle type into pt_typeInfo
!!  
!! PARAMETERS
!!
!!
!!***

!===============================================================================

subroutine pt_prepareEsti (dtOld,dtNew,particles,p_count, ind)
  
  implicit none
  
#include "Flash.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
  
  return
  
end subroutine pt_prepareEsti

