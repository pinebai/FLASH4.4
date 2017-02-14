!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_mgdSetEnergy
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
  use RadTrans_data, ONLY: rt_boltz, rt_speedlt, rt_radconst, rt_acrossMe, &
       rt_meshCopyCount
  use rt_data, ONLY: rt_mgdNumGroups, rt_useMGD
  
  use RadTrans_interface, ONLY: RadTrans_mgdGetBound, RadTrans_planckInt
  use Grid_interface, ONLY: Grid_putPointData, Grid_getPointData
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockId
  integer, intent(in) :: axis(MDIM)
  integer, intent(in) :: grpNum
  real,    intent(in) :: eg

  ! Local Variables:
  integer :: n        ! Loop counter for number of groups
  integer :: g        ! The group number
  integer :: ig       ! The group unk index
  real    :: erad_tot ! The total specific radiation energy in the cell

  ! Loop over all groups that this mesh owns...
  do n = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     ! Get the group number for this group:
     g = NONREP_LOC2GLOB(n, rt_acrossMe, rt_meshCopyCount)
     
     ! Check to see if this mesh owns this group:
     if(g == grpNum) then
        ! Get the UNK index for this group:
        ig = MGDR_NONREP_LOC2UNK(n)
        call Grid_putPointData(blockId, CENTER, ig, EXTERIOR, axis, eg)

        ! Adjust the total cell specific radiation energy for consistency:
        call Grid_getPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, erad_tot)
        erad_tot = erad_tot + eg
        call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, erad_tot)
     end if
  end do

  return

end subroutine RadTrans_mgdSetEnergy
