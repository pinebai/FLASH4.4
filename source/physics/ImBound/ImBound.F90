!!****f* source/physics/ImBound/ImBound
!!
!!
!! NAME
!!
!!  ImBound
!!
!!
!! SYNOPSIS
!!
!!  ImBound()
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine ImBound(blockCount, blockList, dt, forcflag)
  implicit none
#include "Flash.h"
  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  integer, INTENT(IN) :: forcflag
  !! -----------------------------------------------------
  return
end subroutine ImBound
