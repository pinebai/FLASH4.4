!!****if* source/physics/Hydro/HydroMain/split/RHD/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!! SYNOPSIS
!!
!!  Hydro_computeDt(integer(IN):: blockID, 
!!                  real(IN) :: x(:), 
!!                  real(IN) :: dx(:), 
!!                  real(IN) :: uxgrid(:),
!!                  real(IN) :: y(:), 
!!                  real(IN) :: dy(:), 
!!                  real(IN) :: uygrid(:), 
!!                  real(IN) :: z(:), 
!!                  real(IN) :: dz(:), 
!!                  real(IN) :: uzgrid(:), 
!!                  integer(IN) :: blkLimits(2,MDIM)
!!                  integer(IN) :: blkLimitsGC(2,MDIM)
!!                  real,pointer ::  solnData(:,:,:,:),   
!!                  real,(INOUT) ::   dtCheck, 
!!                  integer(INOUT) :: dtMinLoc(:),
!!                  real(INOUT), optional :: extraInfo)
!!
!! DESCRIPTION
!!
!!  Minimum dt is actually computed in hy_rhd_sweep, where all required
!!  information is evaluated anyway. This routine exchanges results
!!  with Hydro through private data members. Most of the arguements
!!  are redundant for this implementation, they are there only
!!  to match the interface. The relevant ones are described below.
!!
!!
!! ARGUMENTS
!!
!!  blockID -       local block ID
!!
!!  __________begin redundant arguments _______________________
!!  x, y, z -       coordinates 
!!  dx, dy, dz -    deltas in each {x, y z} directions
!!  uxgrid, uygrid, uzgrid - velocity of grid expansion in {x, y z} directions
!!  blkLimits -    the indices for the interior endpoints of the block
!!  blkLimitsGC - the indices for endpoints including the guardcells
!!  solnData -      the physical, solution data from grid
!!  __________end redundant arguments __________________________
!!
!!  dtCheck -      variable to hold timestep constraint
!!  dtMinLoc(5) -  array to hold location of cell responsible for minimum dt:
!!                 dtMinLoc(1) = i index
!!                 dtMinLoc(2) = j index
!!                 dtMinLoc(3) = k index
!!                 dtMinLoc(4) = blockID
!!                 dtMinLoc(5) = hy_meshMe
!! extraInfo    -  Driver_computeDt can provide extra info to the caller
!!                 using this argument.
!!
!!
!!****************************************************************************

subroutine Hydro_computeDt(blockID,  &
                           x, dx, uxgrid, &
                           y, dy, uygrid, &
                           z, dz, uzgrid, &
                           blkLimits,     &
                           blkLimitsGC,   &
                           solnData,      &
                           dtCheck, dtMinLoc,&
                           extraInfo )

  use Hydro_data, only : hy_dtmin, hy_dtminloc
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockID 
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dtCheck
  integer,INTENT(INOUT)    :: dtMinLoc(5)
  real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif
  real,OPTIONAL,intent(INOUT) :: extraInfo


  integer, dimension(MAXBLOCKS) :: blockList
  integer :: localNumBlocks, i

#ifndef FLASH_GRID_UG
     call Grid_getListOfBlocks(PARENT_BLK,blockList,localNumBlocks)
  do i=1,localNumBlocks
     if( (hy_dtminloc(4)== blockID) .and. (blockID == blockList(i)) ) then
        hy_dtmin = 0.5*hy_dtmin
     end if
  enddo
#endif

  if( hy_dtmin < dtcheck ) then
    dtcheck  = hy_dtmin
    dtMinloc = hy_dtminloc
  end if

  if(dtCheck <= 0.0) call Driver_abortFlash("[Hydro]: Computed dt is not > 0! Aborting!")

  !! Set default value to be null in not needed.
  if (present(extraInfo)) extraInfo = 0.

end subroutine Hydro_computeDt
