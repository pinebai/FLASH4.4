!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_accExtrnPotPoint
!!
!! NAME
!!
!!  grv_accExtrnPotPoint
!!
!!
!! SYNOPSIS
!!
!!  call grv_accExtrnPotPoint
!!           real(in) :: pos(MDIM)
!!           real(in) :: acc(MDIM)
!!  )
!!
!! DESCRIPTION
!!
!!   Interpolates in the external gravitational field and returns
!!   all components of the acceleration at a given position (pos).
!!
!! ARGUMENTS
!!
!!   pos     - Physical coordinates of the position for which the 
!!             acceleration is returned
!!   acc     - Array to receive result
!!
!! RESULT
!!
!!   Vector of the external gravitational acceleration at position pos.
!!
!!***

subroutine grv_accExtrnPotPoint(pos, acc)
  use Driver_interface, ONLY : Driver_abortFlash
  use Gravity_data, ONLY : grv_bhExtrnPotN, &
    grv_bhExtrnPotCoord, grv_bhExtrnPotPot, grv_bhExtrnPotAcc, &
    grv_bhExtrnPotDel, grv_useExternalPotential, grv_bhExtrnPotIType, &
    grv_bhEPTypeR, grv_bhEPTypeX, grv_bhEPTypeY, grv_bhEPTypeZ, &
    grv_bhExtrnPotCenterX, grv_bhExtrnPotCenterY, grv_bhExtrnPotCenterZ

  implicit none
#include "constants.h"
#include "Flash.h"
  real, dimension(MDIM), intent(in)  :: pos
  real, dimension(MDIM), intent(out) :: acc
  integer :: ind
  real :: delr, Fr, p
  real :: xCenter, yCenter, zCenter

  acc(:) = 0.0
  if (.not. grv_useExternalPotential) return

  select case (grv_bhExtrnPotIType)
    case(grv_bhEPTypeR)
      zCenter = 0.
      yCenter = 0.
      if (NDIM == 3) then 
         zCenter = pos(KAXIS) - grv_bhExtrnPotCenterZ
      endif
      if (NDIM >= 2) then
         yCenter = pos(JAXIS) - grv_bhExtrnPotCenterY
      endif
      xCenter = pos(IAXIS) - grv_bhExtrnPotCenterX

      delr = sqrt(xCenter*xCenter + yCenter*yCenter + zCenter*zCenter)
      ind = floor(delr/grv_bhExtrnPotDel) + 1
      p = (delr - grv_bhExtrnPotCoord(ind)) / (grv_bhExtrnPotCoord(ind+1) - grv_bhExtrnPotCoord(ind))
      Fr = (1.-p)*grv_bhExtrnPotAcc(ind) + p*grv_bhExtrnPotAcc(ind+1)
      acc(IAXIS) = Fr * xCenter/delr
      acc(JAXIS) = Fr * yCenter/delr
      acc(KAXIS) = Fr * zCenter/delr
     
    case(grv_bhEPTypeX)
      call Driver_abortFlash("grv_readExtrnPotential: planex symmetry not supported yet")
    case(grv_bhEPTypeY)
      call Driver_abortFlash("grv_readExtrnPotential: planey symmetry not supported yet")
    case(grv_bhEPTypeZ)
      acc(IAXIS) = 0.0
      acc(JAXIS) = 0.0
      zCenter = pos(KAXIS) - grv_bhExtrnPotCenterZ
     
      ind = floor(abs(zCenter)/grv_bhExtrnPotDel) + 1
      if (ind > grv_bhExtrnPotN) &
      & call Driver_abortFlash("grv_accExtrnPotPoint: out of external pot. array")
      p = (abs(zCenter) - grv_bhExtrnPotCoord(ind)) &
      & / (grv_bhExtrnPotCoord(ind+1) - grv_bhExtrnPotCoord(ind))
      if (zCenter < 0.0) then
        acc(KAXIS) =  (1.-p)*grv_bhExtrnPotAcc(ind) + p*grv_bhExtrnPotAcc(ind+1)
      else
        acc(KAXIS) = -(1.-p)*grv_bhExtrnPotAcc(ind) - p*grv_bhExtrnPotAcc(ind+1)
      endif

    case default
      call Driver_abortFlash("grv_readExtrnPotential: unrecognized potential type")
  end select

end subroutine grv_accExtrnPotPoint


