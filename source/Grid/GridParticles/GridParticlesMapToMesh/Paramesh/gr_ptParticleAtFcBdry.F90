!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptParticleAtFcBdry
!!
!! NAME
!!  gr_ptParticleAtFcBdry
!!
!! SYNOPSIS
!!
!!  gr_ptParticleAtFcBdry(integer,intent(IN) :: partID, &
!!                        logical, intent(OUT) :: fcBdry)
!!
!! DESCRIPTION
!!
!! Routine which determines whether a particle smears across a fine-coarse boundary.  
!! It is not currently used, but the logical return argument could be used to 
!! select a particle mapping scheme that behaves differently when the particle 
!! smears across a fine-coarse boundary.
!! 
!! ARGUMENTS
!!               partID:  Index of particle ID in particles data structure.
!!               fcBdry:  Logical value based on whether a portion of the particle 
!!                        cloud exists on a block at a different refinement.
!! PARAMETERS
!! 
!!***

subroutine gr_ptParticleAtFcBdry(partID, fcBdry)

#include "Flash.h"
#include "constants.h"
#include "gr_ptMapToMesh.h"

  use tree, ONLY : lnblocks, nodetype
  use Particles_data, ONLY : particles
  use gr_ptMapData, ONLY : gr_ptDomain, gr_ptSmearLen
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: partID
  logical, intent(OUT) :: fcBdry

  real, dimension(1:NDIM) :: maxSpread
  real, dimension(LOW:HIGH,1:NDIM) :: blockBoundary
  integer, dimension(1:MDIM) :: guardCellRegion
  integer :: lb, eachAxis, blkNo, blockID, posInd
  integer :: regionID, i, j, k, allCenters
  logical :: blockIDMatch


  blockID = int(particles(BLK_PART_PROP, partID))

  !----------------------------------------------------------------------------
  ! First determine the index, blkNo, into the structure where we can find 
  ! information about blockID.
  !----------------------------------------------------------------------------
  lb=0
  blkNo=0
  blockIDMatch = .false.

  do while (lb < lnblocks)
     lb = lb + 1
     if(nodetype(lb)==LEAF) then
        blkNo = blkNo + 1

        if((gr_ptDomain(blkNo) % blockID) == blockID) then
           blockIDMatch = .true.
           exit
        end if

     end if
  end do

  if (blockIDMatch.eqv..false.) then
     call Driver_abortFlash("gr_ptParticleAtFcBdry error. block ID not found!!!!")
  end if


  !----------------------------------------------------------------------------
  ! Determine which guard cell region that the particle can smear into.
  !----------------------------------------------------------------------------
  maxSpread(1:NDIM) = &
       real(gr_ptSmearLen) * (gr_ptDomain(blkNo) % cellSpacing(1:NDIM))
  blockBoundary(LOW:HIGH,1:NDIM) = gr_ptDomain(blkNo) % bndBlk(LOW:HIGH,1:NDIM)

  posInd = POSX_PART_PROP
  guardCellRegion(1:MDIM) = 1 !For the case when NDIM < MDIM
  do eachAxis = 1, NDIM

     if( (particles(posInd, partID) + maxSpread(eachAxis)) > & 
          blockBoundary(HIGH,eachAxis) ) then

        guardCellRegion(eachAxis) = RIGHT_EDGE

     else if( (particles(posInd, partID) - maxSpread(eachAxis)) < & 
          blockBoundary(LOW,eachAxis) ) then

        guardCellRegion(eachAxis) = LEFT_EDGE 

     else
     
        guardCellRegion(eachAxis) = CENTER

     end if
     posInd = posInd + 1

  end do


  !----------------------------------------------------------------------------
  ! Locate whether this guard cell region has a neighbor at a different 
  ! refinement level.
  !----------------------------------------------------------------------------
  regionID = 0
  allCenters = 2**NDIM
  do i = 1, guardCellRegion(IAXIS)
     do j = 1, guardCellRegion(JAXIS)
        do k = 1, guardCellRegion(KAXIS)
           if((i*j*k)/=allCenters) then 
              regionID = regionID + 1
           end if
        end do
     end do
  end do      


  !If we have at least one neighbor then check the refinement of any one of these neighbors.
  !Each neighbor that corresponds to a certain guard cell region has the same refinement level.
  fcBdry = .false.
  if ( gr_ptDomain(blkNo) % haloRegion(regionID) % numNegh > 0 ) then
     if ( gr_ptDomain(blkNo) % haloRegion(regionID) % neighbor(1) % negh(REFLEVELDIF) /= 0) then
        fcBdry = .true.
     end if
  end if


  return

end subroutine gr_ptParticleAtFcBdry
