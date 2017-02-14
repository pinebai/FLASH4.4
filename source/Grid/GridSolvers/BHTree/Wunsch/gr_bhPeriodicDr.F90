!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPeriodicDr
!!
!! NAME
!!
!!  gr_bhPeriodicDr
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhPeriodicDr()
!!
!! DESCRIPTION
!!   Determines levels up to which individual block trees need to be sent to
!!   different CPUs. For a given CPU, copies all trees up to appropriate levels
!!   to a single message and sends it to the CPU.
!!   
!!
!! ARGUMENTS
!!
!!
!!***



subroutine gr_bhPeriodicDr(dr)

  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhLxHalf, gr_bhLyHalf, gr_bhLzHalf, &
    gr_bhBndType
  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"
  real, dimension(MDIM+2), intent(INOUT) :: dr

  if (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) then
    if (dr(IAXIS) < -gr_bhLxHalf) then
      dr(IAXIS) = dr(IAXIS) + gr_bhLx
    else if (dr(IAXIS) > gr_bhLxHalf) then
      dr(IAXIS) = dr(IAXIS) - gr_bhLx
    endif
    !dr(IAXIS) = min(abs(dr(IAXIS)), abs(dr(IAXIS)+gr_bhLx), abs(dr(IAXIS)-gr_bhLx))
  endif
  if (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) then
    if (dr(JAXIS) < -gr_bhLyHalf) then
      dr(JAXIS) = dr(JAXIS) + gr_bhLy
    else if (dr(JAXIS) > gr_bhLyHalf) then
      dr(JAXIS) = dr(JAXIS) - gr_bhLy
    endif
    !dr(JAXIS) = min(abs(dr(JAXIS)), abs(dr(JAXIS)+gr_bhLy), abs(dr(JAXIS)-gr_bhLy))
  endif
  if (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC) then
    if (dr(KAXIS) < -gr_bhLzHalf) then
      dr(KAXIS) = dr(KAXIS) + gr_bhLz
    else if (dr(KAXIS) > gr_bhLzHalf) then
      dr(KAXIS) = dr(KAXIS) - gr_bhLz
    endif
    !dr(KAXIS) = min(abs(dr(KAXIS)), abs(dr(KAXIS)+gr_bhLz), abs(dr(KAXIS)-gr_bhLz))
  endif


end subroutine gr_bhPeriodicDr
