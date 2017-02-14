!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleFinalize
!!
!! NAME
!!
!!  gr_mpoleFinalize
!!
!! SYNOPSIS
!!
!!  gr_mpoleFinalize ()
!!
!! DESCRIPTION
!!
!!  Finalize the multipole Poisson solver.  Deallocate all storage.
!!
!!***

subroutine gr_mpoleFinalize ()

  use gr_mpoleData, ONLY : gr_mpoleMaxL,                  &
                           gr_mpoleNumberInv,             &
                           gr_mpoleZoneRmax,              &
                           gr_mpoleZoneQmax,              &
                           gr_mpoleZoneType,              &
                           gr_mpoleZoneScalar,            &
                           gr_mpoleZoneLogNorm,           &
                           gr_mpoleZoneExponent,          &
                           gr_mpoleZoneScalarInv,         &
                           gr_mpoleZoneLogNormInv,        &
                           gr_mpoleZoneExponentInv,       &
                           gr_mpoleZoneMaxRadiusFraction, &
                           gr_mpoleBlockList

  implicit none

  if (allocated (gr_mpoleZoneRmax)             ) deallocate (gr_mpoleZoneRmax)
  if (allocated (gr_mpoleZoneQmax)             ) deallocate (gr_mpoleZoneQmax)
  if (allocated (gr_mpoleZoneType)             ) deallocate (gr_mpoleZoneType)
  if (allocated (gr_mpoleZoneScalar)           ) deallocate (gr_mpoleZoneScalar)
  if (allocated (gr_mpoleZoneLogNorm)          ) deallocate (gr_mpoleZoneLogNorm)
  if (allocated (gr_mpoleZoneExponent)         ) deallocate (gr_mpoleZoneExponent)
  if (allocated (gr_mpoleZoneScalarInv)        ) deallocate (gr_mpoleZoneScalarInv)
  if (allocated (gr_mpoleZoneLogNormInv)       ) deallocate (gr_mpoleZoneLogNormInv)
  if (allocated (gr_mpoleZoneExponentInv)      ) deallocate (gr_mpoleZoneExponentInv)
  if (allocated (gr_mpoleZoneMaxRadiusFraction)) deallocate (gr_mpoleZoneMaxRadiusFraction)
  if (allocated (gr_mpoleBlockList)            ) deallocate (gr_mpoleBlockList)

  if (gr_mpoleMaxL > 0) then
      if (allocated (gr_mpoleNumberInv)) deallocate (gr_mpoleNumberInv)
  end if

  return
end subroutine gr_mpoleFinalize
