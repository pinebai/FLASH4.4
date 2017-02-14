!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/fl_fsTFIFlameSpeedBlock
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
!! This subroutine implements Colin et al. (2000)
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   s : 
!!
!!   ds : 
!!
!!   dx : 
!!
!!   complimits : 
!!
!!   quench_limit : 
!!
!!
!!
!!***


subroutine fl_fsTFIFlameSpeedBlock(solnData, s, ds, dx, compLimits, &
                                   quench_limit)

#include "constants.h"
#include "Flash.h"

  use Flame_interface, ONLY : Flame_getWidth
  use Turb_interface, ONLY : Turb_getFilterScale
  use fl_fsTFIInterface, ONLY : fl_fsTFIEnhance

  implicit none

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(:,:,:), intent(inout) :: s
  real, dimension(:,:,:), intent(in) :: ds
  real, intent(in) :: dx
  integer, dimension(LOW:HIGH,MDIM), intent(in) :: compLimits
  real, dimension(:,:,:), intent(out), optional :: quench_limit

  integer :: istat

  integer :: i,j,k
  real :: E, up, dl0, dl1, de_over_dl1
  real :: de

  ! get filter scale associated with turbulence operator
  call Turb_getFilterScale(dx, de)

  ! get thickened flame width
  call Flame_getWidth(dl1)

  de_over_dl1 = de / dl1

  do k = compLimits(LOW,KAXIS), compLimits(HIGH,KAXIS)
     do j = compLimits(LOW,JAXIS), compLimits(HIGH,JAXIS)
        do i = compLimits(LOW,IAXIS), compLimits(HIGH,IAXIS)

           ! OP2 from Colin et al. 2000
           ! this is the turbulent component of the velocity field
           ! stored in TURB_VAR by Turb_calc
           up = solnData(TURB_VAR,i,j,k)

           ! get enhancement factor
           if (present(quench_limit)) then
              call fl_fsTFIEnhance(E, up, s(i,j,k), de, ds(i,j,k), & 
                                   de_over_dl1, quench_limit(i,j,k))
              quench_limit(i,j,k) = quench_limit(i,j,k) * s(i,j,k)
           else
              call fl_fsTFIEnhance(E, up, s(i,j,k), de, ds(i,j,k), & 
                                   de_over_dl1)
           endif

           ! apply enhancement
           s(i,j,k) = E * s(i,j,k)
           
        enddo
     enddo
  enddo

  return
end subroutine fl_fsTFIFlameSpeedBlock
