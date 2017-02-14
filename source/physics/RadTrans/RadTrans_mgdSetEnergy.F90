!!****f* source/physics/RadTrans/RadTrans_mgdSetEnergy
!!
!!  NAME 
!!
!!  RadTrans_mgdSetEnergy
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdSetEnergy( integer(IN) :: blockId,
!!                              integer(IN) :: axis(MDIM),
!!                              integer(IN) :: grpNum,
!!                              real(IN)    :: eg )
!!
!!  DESCRIPTION 
!!
!!      Set the specific energy for a particular energy group in a
!!      particular cell. This routine has been created to make it easy
!!      for users to specify the energy in a group. This can be a
!!      little complicated because of mesh replication - but all of
!!      the details are handled internally in RadTrans
!!
!! ARGUMENTS
!!
!!    blockId : The blockId of the cell
!!    axis    : An array storing the i,j,k coordinate of the cell
!!    grpNum  : The energy group number
!!    eg      : The specific internal energy to use [ergs/g]
!! 
!!***
subroutine RadTrans_mgdSetEnergy(blockId, axis, grpNum, eg)
  implicit none

#include "constants.h"

  integer, intent(in) :: blockId
  integer, intent(in) :: axis(MDIM)
  integer, intent(in) :: grpNum
  real,    intent(in) :: eg

  return

end subroutine RadTrans_mgdSetEnergy
