!!****f* source/physics/RadTrans/RadTrans_computeFluxLimiter
!!
!! NAME
!!
!!  RadTrans_computeFluxLimiter
!!
!! SYNOPSIS
!!
!!  call RadTrans_computeFluxLimiter(integer(in) :: ifl,
!!                                   integer(in) :: iflOut,
!!                                   integer(in) :: ieddi3,
!!                                   real(INOUT) :: solnData,
!!                                   integer(IN) :: blockID,
!!                                   integer(IN),OPTIONAL :: gcLayers)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ifl : 
!!
!!   iflOut : 
!!
!!   ieddi3 : 
!!
!!   solnData : 
!!
!!   blockID : ID of block in current processor
!!
!!   gcLayers : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine RadTrans_computeFluxLimiter(ifl, iflOut, ieddi3, solnData, blockID,gcLayers)
  implicit none

  integer, intent(in) :: ifl
  integer, intent(in) :: iflOut
  integer, intent(in) :: ieddi3
  real,    intent(INOUT) :: solnData(:,1:,1:,1:)
  integer, intent(IN) :: blockID
  integer, intent(IN),OPTIONAL :: gcLayers

end subroutine RadTrans_computeFluxLimiter
