!!****f* source/physics/IncompNS/IncompNS
!!
!!
!! NAME
!!
!!  IncompNS
!!
!!
!! SYNOPSIS
!!
!!  IncompNS(integer(IN) :: blockCount, 
!!      integer(IN) :: blockList(blockCount)
!!      real(IN)    :: timeEndAdv,
!!      real(IN)    :: dt,
!!      real(IN)    :: dtOld,
!!      integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs INS timestep advancement.
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
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - dummy consistent with toplayer stub function
!!  dt         - timestep
!!  dtOld      - dummy consistent with toplayer stub function
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a top-layer stub function
!!
!!***

subroutine IncompNS( blockCount, blockList, &
                     timeEndAdv, dt,  dtOld,&
                     sweepOrder)

  implicit none

#include "Flash.h"

  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(IN) :: sweepOrder
  integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld


end subroutine IncompNS
