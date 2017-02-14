!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS
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
!!  numProcs   - total number of processors
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - dummy consistent with toplayer stub function
!!  dt         - timestep
!!  dtOld      - dummy consistent with toplayer stub function
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***

subroutine IncompNS( blockCount, blockList, &
                     timeEndAdv, dt,  dtOld,&
                     sweepOrder)

  use ins_interface, ONLY : ins_ab2rk3
  use Driver_interface, ONLY : Driver_getNStep
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use IncompNS_data, ONLY : ins_useIncompNS, ins_nstep

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(IN) :: sweepOrder
  integer, INTENT(INOUT) :: blockList(MAXBLOCKS) 
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld
  integer :: x, i

  if (.NOT. ins_useIncompNS) RETURN

  call Driver_getNStep(ins_nstep)

  call Timers_start("ins_ab2rk3")

  call ins_ab2rk3(  blockCount, blockList, timeEndAdv, dt)

  call Timers_stop("ins_ab2rk3")


end subroutine IncompNS
