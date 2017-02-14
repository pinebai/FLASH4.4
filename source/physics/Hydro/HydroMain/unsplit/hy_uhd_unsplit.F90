!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_unsplit
!!
!! NAME
!!
!!  hy_uhd_unsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_unsplit( integer (IN) :: blockCount,
!!                  integer (IN) :: blockList(blockCount),
!!                  real    (IN) :: dt,
!!                  real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!
!!  Performs MHD/Hydro update in a directionally unsplit fashion over a set
!!  of blocks. Divergence cleaning control for the magnetic fields is
!!  handled by using the staggered mesh algorithm with a call to
!!  hy_uhd_staggeredDivb. Before calling this routine, electric fields are
!!  to be calculated with a call to hy_uhd_getElectricFields.
!!  blockList is an integer array of size blockCount that contains 
!!  the blocks over which to update.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - applies an eos to the guard cells; 
!!   - computes fluxes using a call to hy_uhd_getFaceFlux
!!   - if we're not doing flux correction (as controlled by the flux_correct
!!     runtime parameter), then we update all the cell values from the fluxes 
!!     (with a call to hy_uhd_unsplitUpdate), otherwise, we update just cells 
!!     not on the boundaries, and save fluxes for cells on the boundary;
!!   - and finally, we apply an eos to the block.
!! 
!!  After the main block loop, if doing flux correction, we have
!!  the Grid correct boundary fluxes for all blocks where approriate,
!!  and do another loop over blocks, updating the cell values for
!!  cells on the block boundaries using the corrected fluxes, and
!!  apply an eos on the block. 
!!  The same is true for the electric fields correction where
!!  there are block boundaries that are sharing different levels of
!!  refinements.
!!
!!  This implementation is a stub implementation at this level.
!!
!! REFERENCES
!!
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (needed for temporal extrapolations of gravity)
!!
!!***

subroutine hy_uhd_unsplit ( blockCount, blockList, dt, dtOld )

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "UHD.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  real,    INTENT(IN) :: dt, dtOld
  !! -----------------------------------------------------


End Subroutine hy_uhd_unsplit
