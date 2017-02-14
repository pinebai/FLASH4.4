!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/fl_fsTFIFlameSpeedBlock
!!
!! NAME
!!
!!  fl_fsTFIFlameSpeedBlock
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIFlameSpeedBlock(real, dimension(:,:,:,:), pointer  :: solndata,
!!                               real, dimension(:,:,:)(inout) :: s,
!!                               real, dimension(:,:,:)(in) :: ds,
!!                               real(in) :: dx,
!!                               integer, dimension(LOW:HIGH,MDIM)(in) :: complimits,
!!                               real, dimension(:,:,:)(out), optional  :: quench_limit)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!!
!! The default behavior is to do nothing
!! the laminar flame speed is supplied and the turbulent
!! flame speed is returned.
!!
!! ARGUMENTS
!!
!!   solndata : solution data
!! 
!!   s : not used 
!!
!!   ds : not used 
!!
!!   dx : not used
!!
!!   complimits : not used
!!
!!   quench_limit : quench limit array
!!
!!
!!
!!***


subroutine fl_fsTFIFlameSpeedBlock(solnData, s, ds, dx, compLimits, &
                                   quench_limit)

#include "Flash.h"
#include "constants.h"

  implicit none

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(:,:,:), intent(inout) :: s
  real, dimension(:,:,:), intent(in) :: ds
  real, intent(in) :: dx
  integer, dimension(LOW:HIGH,MDIM), intent(in) :: compLimits
  real, dimension(:,:,:), intent(out), optional :: quench_limit

  if (present(quench_limit)) quench_limit(:,:,:) = -1.0e0

end subroutine fl_fsTFIFlameSpeedBlock
