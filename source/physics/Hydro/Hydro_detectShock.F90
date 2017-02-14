!!****f* source/physics/Hydro/Hydro_detectShock
!!
!!
!! NAME
!!
!!  Hydro_detectShock
!!
!!
!! SYNOPSIS
!!
!!  call Hydro_detectShock(real(IN),pointer :: solnData(:,:,:,:),
!!                         real(OUT)        :: shock(:,:,:),
!!                         integer(IN)      :: blkLimits(2,MDIM),
!!                         integer(IN)      :: blkLimitsGC(2,MDIM),
!!                         integer(IN)      :: guardCells(MDIM),
!!                         real(IN)         :: primaryCoord(:),
!!                         real(IN)         :: secondCoord(:),
!!                         real(IN)         :: thirdCoord(:))
!!
!!
!! DESCRIPTION
!!
!!
!!  Hydro_detectShock is a multidimensional shock detection algorithm.  This
!!  is currently used by the burner to cut off burning in shocks if desired.
!!  This is helpful for detonations in some circumstances.
!!
!!
!! ARGUMENTS
!!
!!  solnData --     pointer to a block of solution data to inspect
!!  shock --        an array indicating if there is a shock in a zone (=1) =0 otherwise
!!  blkLimits  --   index limits of block, interior cells only; see Grid_getBlkIndexLimits 
!!  blkLimitsGC --  index limits of block, including guard cells.   
!!  guardCells  --  number of layers of guard cells which output is requested in addition
!!                  to interior cells
!!  primaryCoord -- x coordinate of solnData, i.e., coordinate in IAXIS direction, including guardcells
!!  secondCoord --  y coordinate of solnData, i.e., coordinate in JAXIS direction, including guardcells
!!  thirdCoord --   z coordinate of solnData, i.e., coordinate in KAXIS direction, including guardcells
!!
!!  NOTE
!!
!!  The guard cells need to be filled before calling this routine.
!!  The pressure needs to be updated before calling this routine
!!
!!  {For future reference : The stub of this routine used to include 
!!   a final argument called    integer, intent(IN)  :: igc)
!!   if =1, include the first row of guardcells in the detection
!!   However, the implementation does not! May need to be addressed in
!!   future}
!!  
!!  SEE ALSO
!!
!!  Grid_getBlkIndexLimits
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif

subroutine Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, &
                             guardCells, &
                             primaryCoord,secondCoord,thirdCoord)


  implicit none
#include "constants.h"


  integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, intent(IN) :: guardCells(MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData 


#ifdef FIXEDBLOCKSIZE
  real, intent(out),dimension(GRID_ILO_GC:GRID_IHI_GC,&
                              GRID_JLO_GC:GRID_JHI_GC,&
                              GRID_KLO_GC:GRID_KHI_GC):: shock
  real,intent(IN),dimension(GRID_ILO_GC:GRID_IHI_GC) :: primaryCoord
  real,intent(IN),dimension(GRID_JLO_GC:GRID_JHI_GC) :: secondCoord
  real,intent(IN),dimension(GRID_KLO_GC:GRID_KHI_GC) :: thirdCoord
#else
  real,intent(out), dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
  real,intent(IN),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) :: primaryCoord
  real,intent(IN),dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) :: secondCoord
  real,intent(IN),dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: thirdCoord

#endif

  shock = 0.0

end subroutine Hydro_detectShock




