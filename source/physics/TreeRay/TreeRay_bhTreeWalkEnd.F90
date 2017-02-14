!!****f* source/physics/TreeRay/TreeRay_bhTreeWalkEnd
!!
!! NAME
!!
!!  TreeRay_bhTreeWalkEnd
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine TreeRay_bhTreeWalkEnd(iterate)
  implicit none
#include "constants.h"
#include "Flash.h"
  logical,intent(OUT) :: iterate

  iterate = .false.
  return
end subroutine TreeRay_bhTreeWalkEnd

