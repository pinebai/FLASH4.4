!!****f* source/physics/sourceTerms/Flame/Flame_computeAbarZbar
!!
!! NAME
!!
!!  Flame_computeAbarZbar
!!
!! SYNOPSIS
!!
!!  Flame_computeAbarZbar( real(in),dimension(:,:)  :: solnScalars,
!!                         real(inout),dimension(:) :: abarData, zbarData)
!!
!! DESCRIPTION
!!
!!  A callback, typically called by Eos unit implementations to get
!!  values for EOS_ABAR and EOS_ZBAR input elements of eosData
!!  before the EOS computation proper.
!!
!! ARGUMENTS
!!
!!   solnScalars - 2D array of simulation mass scalars etc used to calculate abar, zbar
!!                 same order as in main grid arrays
!!      abarData - returned array of abar values
!!      zbarData - returned array of zbar values
!!
!! SEE ALSO
!!
!!  See Flame_interface.F90 for possible updates
!!
!!***

! this is a stub for when the Flame Unit is not included
!
subroutine Flame_computeAbarZbar(solnScalars, abarData, zbarData)

  implicit none

  real, dimension(:,:), intent(in)  :: solnScalars
  real, dimension(:), intent(inout) :: abarData, zbarData

end subroutine Flame_computeAbarZbar
