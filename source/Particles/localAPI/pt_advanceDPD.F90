!!****if* source/Particles/localAPI/pt_advanceDPD
!!
!! NAME
!!
!!  pt_advanceDPD
!!
!! SYNOPSIS
!!
!!  pt_advanceDPD(real(IN) :: dtOld,
!!                  real(IN) :: dtNew,
!!                  real(inout):: particles(:,p_count),
!!                  integer(in):: p_count)
!!                  integer(in):: ind)
!!
!! ARGUMENTS
!!  
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind    --- index for type into pt_typeInfo
!!  
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version is the forward Euler advancement for the active 
!!  submodule.
!!
!!***

!===============================================================================

subroutine pt_advanceDPD (dtOld,dtNew,particles,p_count, ind)
    
  implicit none

#include "Flash.h"
    
  integer       :: i
  real, INTENT(in) :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles

  return

end subroutine pt_advanceDPD

!===============================================================================

