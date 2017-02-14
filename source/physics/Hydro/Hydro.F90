!!****f* source/physics/Hydro/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!! SYNOPSIS
!!
!!  Hydro( integer(IN) :: blockCount, 
!!         integer(IN) :: blockList(blockCount), 
!!         real(IN)    :: timeEndAdv, 
!!         real(IN)    :: dt, 
!!         real(IN)    :: dtOld, 
!!         integer(IN) :: sweepOrder )
!!
!! DESCRIPTION
!!
!!  Perform a 1, 2, or 3D hydro update.  This version handles
!!  directionally split hydro schemes.  The input
!!  parameter sweepOrder determines the ordering of sweep
!!  directions.  For example, in 3d, SWEEP_XYZ means to perform 
!!  1d sweeps first in the x direction, then the y, then z direction, 
!!  while SWEEP_ZYX means to performs the 1d sweeps in the opposite order.
!!
!!  The blockList and blockCount arguments tell this routine on 
!!  which blocks and on how many to operate.  blockList is an 
!!  integer array of size blockCount that contains the local 
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.  
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList -   array holding local IDs of blocks on which to advance
!!  timeEndAdv -  end time
!!  dt -          timestep
!!  dtOld -       old timestep
!!  sweepOrder -  direction of hydro sweep, can be: SWEEP_XYZ or SWEEP_ZYX
!!                as defined in  constants.h
!!
!!
!!***


subroutine Hydro(  blockCount, blockList, &
                   timeEndAdv, dt, dtOld, &
                   sweepOrder )
implicit none
#include "Flash.h"
  
  integer, INTENT(IN) :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld
  integer, INTENT(IN) :: sweepOrder
  
end subroutine Hydro
