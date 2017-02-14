!!****f* source/physics/Hydro/Hydro_shockStrength
!!
!!
!! NAME
!!
!!  Hydro_shockStrength
!!
!!
!! SYNOPSIS
!!  call Hydro_shockStrength(real(IN),pointer :: solnData(:,:,:,:), 
!!                     real(INOUT)          :: shock(:,:,:), 
!!                     integer(IN)          :: blkLimits(2,MDIM),
!!                     integer(IN)          :: blkLimitsGC(2,MDIM),
!!                     integer(IN)          :: guardCells(MDIM),
!!                     real(IN)             :: primaryCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),
!!                     real(IN)             :: secondCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),
!!                     real(IN)             :: thirdCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),
!!                     real(IN)             :: threshold,
!!                     integer(IN)          :: mode)
!!
!! DESCRIPTION
!!
!!  Hydro_shockStrength computes a measure of shock strength for each cell
!!  in a block. 
!!
!!  Hydro_shockStrength includes the same multidimensional shock detection algorithm as
!!  Hydro_detectShhock, from which it is derived. 
!!
!!
!! ARGUMENTS
!!
!!  solnData --     pointer to a block of solution data to inspect
!!  shock --        an array returning a measure of shock strenght for cells where a shock (of sufficient
!!                  strength, see threshold) is detected.
!!                  A combination of current shock strenght with the previous value may be returned
!!                  if mode > 1, see mode below.
!!  blkLimits  --   index limits of block, interior cells only; see Grid_getBlkIndexLimits 
!!  blkLimitsGC --  index limits of block, including guard cells.   
!!  guardCells  --  number of layers of guard cells which output is requested in addition
!!                  to interior cells
!!  primaryCoord -- x coordinate of solnData, i.e., coordinate in IAXIS direction, including guardcells
!!  secondCoord --  y coordinate of solnData, i.e., coordinate in JAXIS direction, including guardcells
!!  thirdCoord --   z coordinate of solnData, i.e., coordinate in KAXIS direction, including guardcells
!!  threshold -     threshold value for shock strength
!!  mode      -     1 to forget previous values in shock array, 2 to add to prev value,
!!                  3 for max of previous and current value
!!  
!!
!!  NOTE
!!
!!  The guard cells need to be filled before calling this routine.
!!  The pressure needs to be updated before calling this routine
!!  
!!  SEE ALSO
!!
!!  Grid_getBlkIndexLimits
!!
!!***


subroutine Hydro_shockStrength(solnData, shock, blkLimits, blkLimitsGC, &
                             guardCells, &
                             primaryCoord,secondCoord,thirdCoord, &
                             threshold, mode)
  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, intent(IN) :: guardCells(MDIM)
  real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, intent(inout),dimension(GRID_ILO_GC:GRID_IHI_GC,&
                              GRID_JLO_GC:GRID_JHI_GC,&
                              GRID_KLO_GC:GRID_KHI_GC):: shock
  real,intent(IN),dimension(GRID_ILO_GC:GRID_IHI_GC) :: primaryCoord
  real,intent(IN),dimension(GRID_JLO_GC:GRID_JHI_GC) :: secondCoord
  real,intent(IN),dimension(GRID_KLO_GC:GRID_KHI_GC) :: thirdCoord
#else
  real,intent(inout),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
  real,intent(IN),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) :: primaryCoord
  real,intent(IN),dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) :: secondCoord
  real,intent(IN),dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: thirdCoord
#endif
  real, intent(IN) :: threshold
  integer, intent(IN) :: mode
 
 
end subroutine Hydro_shockStrength




