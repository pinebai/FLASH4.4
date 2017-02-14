!!****f* source/Particles/Particles_longRangeForce
!!
!! NAME
!!
!!  Particles_longRangeForce()
!!
!! SYNOPSIS
!!
!!  Particles_longRangeForce(real,intent(inout) :: particles(:,:),
!!                           integer,intent(in) :: p_count,
!!                           integer,intent(in) :: mapType)
!!
!! DESCRIPTION
!!
!!  Computes long-range forces on particles, ie. forces which couple all
!!  particles to each other.  This version is for particle-mesh gravitation and
!!  maps the gravitational acceleration on the mesh to the particle positions.
!!  
!! ARGUMENTS
!!
!!   particles :: the list of particles to be operated on
!!   p_count   :: count of the particles in the list
!!   mapType   :: when mapping grid quantities to particle, method to use
!!  
!! PARAMETERS
!!
!!***

!===============================================================================

subroutine Particles_longRangeForce (particles,p_count,mapType)

    
!-------------------------------------------------------------------------------
#include "Flash.h"

  implicit none

!-------------------------------------------------------------------------------
  integer, intent(IN) :: p_count,mapType
  real,dimension(NPART_TYPES,p_count),intent(INOUT) :: particles

!-------------------------------------------------------------------------------

  return

end subroutine Particles_longRangeForce

!===============================================================================
