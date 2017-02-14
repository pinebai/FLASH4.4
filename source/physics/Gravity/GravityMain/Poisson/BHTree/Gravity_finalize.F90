!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_finalize
!!
!! NAME
!!
!!  Gravity_finalize
!!  
!! SYNOPSIS
!!
!!  Gravity_finalize()
!!
!! DESCRIPTION
!!  deallocate any memory that might have allocated in Gravity unit
!!  and other housekeeping to prepare for the end of the unit.
!!
!!
!!***
subroutine Gravity_finalize()

  use Gravity_data, ONLY : grv_bhExtrnPotCoord, grv_bhExtrnPotPot, &
    grv_bhExtrnPotAcc, grv_bhTreeEwaldAccV42, grv_bhTreeEwaldPotV42, &
    grv_bhLogfield, grv_bhLogfieldDer, grv_bhTreeEwald
  implicit none


  if (allocated(grv_bhExtrnPotCoord)) deallocate(grv_bhExtrnPotCoord)
  if (allocated(grv_bhExtrnPotPot)) deallocate(grv_bhExtrnPotPot)
  if (allocated(grv_bhExtrnPotAcc)) deallocate(grv_bhExtrnPotAcc)
  if (allocated(grv_bhTreeEwaldAccV42)) deallocate(grv_bhTreeEwaldAccV42)
  if (allocated(grv_bhTreeEwaldPotV42)) deallocate(grv_bhTreeEwaldPotV42)
  if (allocated(grv_bhTreeEwald)) deallocate(grv_bhTreeEwald)
  if (allocated(grv_bhLogfield)) deallocate(grv_bhLogfield)
  if (allocated(grv_bhLogfieldDer)) deallocate(grv_bhLogfieldDer)

  return

end subroutine Gravity_finalize
