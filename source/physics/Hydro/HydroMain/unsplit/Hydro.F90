!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(integer(IN) :: blockCount, 
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs physics update in a directionally unsplit fashion.
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
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***

Subroutine Hydro( blockCount, blockList, &
                  timeEndAdv, dt,  dtOld,&
                  sweepOrder)

  use Hydro_data,       ONLY : hy_useHydro, hy_gpotAlreadyUpToDate
  use hy_uhd_interface, ONLY : hy_uhd_unsplit
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(IN) :: blockCount,sweepOrder
  integer, INTENT(IN) :: blockList(blockCount)
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  hy_gpotAlreadyUpToDate = .FALSE. ! reset this flag, may be set .TRUE. below if warranted.

  if (.not. hy_useHydro) return

  call Timers_start("hydro_unsplit")

  call hy_uhd_unsplit(blockCount, blockList, dt, dtOld)

  call Timers_stop("hydro_unsplit")


End Subroutine Hydro
