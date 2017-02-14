!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/Hydro
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
!!  Top level function that switches directions of mhd sweeps
!!  that are specific to an 8-wave splitting algorithm.
!!
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - sweep order
!!
!!
!! NOTE
!!
!!  Most of the arguments are dummy in this routine and they are present
!!  to be consistent with the very top layer stub function.
!!
!!***

subroutine Hydro(  blockCount, blockList, &
                  timeEndAdv, dt, dtOld, sweepOrder )


  use Hydro_data, ONLY : hy_killDivb
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use hy_8wv_interface, ONLY : hy_8wv_divb, hy_8wv_sweep

  implicit none

#include "constants.h"
#include "Flash.h"


  integer, INTENT(in) ::  blockCount, sweepOrder
  integer, INTENT(in), dimension(blockCount) :: blockList
  real,    INTENT(in) :: timeEndAdv, dtOld
  real,    INTENT(inout) :: dt

  call Timers_start("MHD")
  call Timers_start("MHD_sweep")


#if NDIM >= 2

  if( hy_killdivb ) then
     call Timers_start("divb")
    call hy_8wv_divb(blockCount,blockList,dt)
    call Timers_stop("divb")
  end if

#endif


  select case (sweepOrder)
  case (SWEEP_XYZ)
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_X)
#if NDIM >= 2
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Y)
#if NDIM == 3
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Z)
#endif
#endif


  case (SWEEP_ZYX)
#if NDIM >= 2
#if NDIM == 3
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Z)
#endif
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Y)
#endif
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_X)


  case (SWEEP_XZY)
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_X)
#if NDIM >= 2
#if NDIM == 3
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Z)
#endif
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Y)
#endif


  case (SWEEP_YZX)
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Y)
#if NDIM >= 2
#if NDIM == 3
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Z)
#endif
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_X)
#endif


  case (SWEEP_YXZ)
#if NDIM >= 2
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Y)
#endif
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_X)
#if NDIM == 3
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Z)
#endif


  case (SWEEP_ZXY)
#if NDIM == 3
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Z)
#endif
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_X)
#if NDIM >= 2
     call hy_8wv_sweep( blockCount, blockList, dt, SWEEP_Y)
#endif

  end select


  call Timers_stop("MHD_sweep")

#if NDIM >= 2

  if( hy_killdivb ) then
    call Timers_start("divb")
    call hy_8wv_divb(blockCount,blockList,dt)
    call Timers_stop("divb")
  end if

#endif

  call Timers_stop ("MHD")

end subroutine Hydro
