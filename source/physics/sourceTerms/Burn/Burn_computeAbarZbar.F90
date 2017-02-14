!!****f* source/physics/sourceTerms/Burn/Burn_computeAbarZbar
!!
!! NAME
!!
!!  Burn_computeAbarZbar
!!
!! SYNOPSIS
!!
!!  call Burn_computeAbarZbar(real, dimension(:,:)(in) :: solnscalars,
!!                            real, dimension(:)(inout) :: abardata,
!!                            real, dimension(:)(inout) :: zbardata)
!!
!! DESCRIPTION
!!
!! Stub
!!
!! ARGUMENTS
!!
!!   solnscalars : solution scalars
!!
!!   abardata : abar info 
!!
!!   zbardata : zbar info
!!
!!
!!
!!***

subroutine Burn_computeAbarZbar(solnScalars, abarData, zbarData)

  implicit none

  real, dimension(:,:), intent(in)  :: solnScalars
  real, dimension(:), intent(inout) :: abarData, zbarData

end subroutine Burn_computeAbarZbar
