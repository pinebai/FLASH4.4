!!****if* source/Simulation/SimulationMain/Orbit/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer,intent(IN)  :: blockID)
!!
!! DESCRIPTION Initializes the fluid variables in the zones of a specified
!!             block.  This routine is called by the initialization code
!!             after the grid is set up.
!!             This version sets up an approximate gridded particle density
!!             field for the two-particle orbit problem.  The approximate
!!             density guides the initial refinement of the mesh before
!!             the particles themselves are initialized.
!!
!!
!! ARGUMENTS
!!
!!  blockID -       The number of the block to initialize
!!
!!
!!***


subroutine Simulation_initBlock(blockID)

!===================================================================

  use Simulation_data, ONLY: sim_xPosns, sim_yPosns, sim_zPosns, &
       sim_nPtot, sim_EPSX, sim_EPSY, sim_EPSZ
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
    Grid_getBlkCenterCoords, Grid_getBlkIndexLimits, Grid_putPointData

  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockID
  

  real,dimension(MDIM) :: blockSize, blockCenter,pos
  real                 :: bxl, byl, bzl, bxu, byu, bzu
  integer              :: p, i, j, k
  integer,dimension(MDIM)   :: ipos,localSize,guard
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  !===============================================================================



  ! Get block dimensions
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)

  !get the index limits. blkLimits holds interior cell indices 
  !blkLimitsGC holds cell indices including guardcells
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  !holds interior block size in each dimension
  !(if using paramesh, values are just NXB, NYB and NZB)
  localSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1

  !holds nguard cells in each dimension
  !(if using paramesh just, NGUARD for all 3 entries)
  guard = blkLimits(LOW,:)-blkLimitsGC(LOW,:)


  bxl    = blockCenter(IAXIS) - 0.5*blockSize(IAXIS)
  bxu    = blockCenter(IAXIS) + 0.5*blockSize(IAXIS)
  if (NDIM >= 2) then
    byl = blockCenter(JAXIS) - 0.5*blockSize(JAXIS)
    byu = blockCenter(JAXIS) + 0.5*blockSize(JAXIS)
  else
    byl = 0.
    byu = 0.
  endif
  if (NDIM == 3) then
    bzl = blockCenter(KAXIS) - 0.5*blockSize(KAXIS)
    bzu = blockCenter(KAXIS) + 0.5*blockSize(KAXIS)
  else
    bzl = 0.
    bzu = 0.
  endif



! Set approximate gridded particle density

  do p = 1, sim_nPtot

    if ((sim_xPosns(p) >= bxl-sim_EPSX) .and. (sim_xPosns(p) <= bxu+sim_EPSX) .and. &
        (sim_yPosns(p) >= byl-sim_EPSY) .and. (sim_yPosns(p) <= byu+sim_EPSY) .and. &
        (sim_zPosns(p) >= bzl-sim_EPSZ) .and. (sim_zPosns(p) <= bzu+sim_EPSZ)) then

       pos(IAXIS) = guard(IAXIS) + &
                    (sim_xPosns(p)-bxl)/blockSize(IAXIS)*localSize(IAXIS)
       pos(JAXIS) = guard(JAXIS) + &
                    (sim_yPosns(p)-byl)/blockSize(JAXIS)*localSize(JAXIS)
       pos(KAXIS) = guard(KAXIS) + &
                    (sim_zPosns(p)-bzl)/blockSize(KAXIS)*localSize(KAXIS)


       ipos = ceiling(pos)

       call Grid_putPointData(blockID, CENTER, PDEN_VAR, EXTERIOR, ipos, 1.)
       
    endif
    
 enddo
 
 !===============================================================================
 
 return
 
end subroutine Simulation_initBlock
