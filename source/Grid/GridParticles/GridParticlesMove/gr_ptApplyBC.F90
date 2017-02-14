!!****if* source/Grid/GridParticles/GridParticlesMove/gr_ptApplyBC
!!
!! NAME
!!
!!  gr_ptApplyBC
!!
!! SYNOPSIS
!!
!!  gr_ptApplyBC(real(INOUT)     :: dataElement(propCount),
!!               integer(IN)     :: propCount,
!!               integer(IN)     :: blockID,
!!               integer(INOUT)  :: lostElements,
!!               logical(INOUT)  :: moved,
!!               integer(INOUT)  :: negh(:))
!!
!! DESCRIPTION
!!   This routine figures out if a particle has hit a global domain boundary.  If
!!   so it applies the correct boundary conditions.  
!!
!!   On return, the particle may be changed in one of the following ways:
!!    o  The particle is marked for dropping by setting its associated blk to NONEXISTENT.
!!       When this happens, lostElements is also decremented by 1.
!!    o  The POS{X,Y,Z}and possibly other properties like
!!       POSPRED{X,Y,Z}, VEL{X,Y,Z}, and/or VELPRED{X,Y,Z}
!!       are changed so as to implemwent the action of a "periodic" or "reflecting"
!!       or similar boundary.
!!
!!   It is left to the caller to take appropriate action based on these changes, such as
!!   freeing storage space for a dropped particle or assigning a particle to a different
!!   block (possibly on a different processor).
!!
!! ARGUMENTS
!!
!!   dataElement : Data structure holding the properties of one dataElement.
!!
!!   propCount : the count of properties in the dataElement data structure
!!
!!   blockID   : the id of the block number associated with the dataElement
!!   lostElements : counter for elements that go permanently missing.
!!   moved      : is true if the dataElement has moved from current block
!!   negh       : information about the neighbor
!!
!! NOTES
!!   If a setup is using user defined boundary conditions for the fluid, it is very likely
!!   to need a customized version of this routine in the setup directory.
!!
!!   The current implementation does not deal correctly with boundaries at inner boundary
!!   blocks defined with Simulation_defineDomain. Use of this GridParticles implementation
!!   together with Simulation_defineDomain may thus result in undefined behavior when a
!!   elements moves across a boundary block boundary.
!!
!!   In this implementation the following boundary conditions are specifically recognized:
!!
!!       Boundary Condition         Action
!!   ================================================================================
!!       OUTFLOW                    Drop particle
!!       DIODE
!!   --------------------------------------------------------------------------------
!!       REFLECTING                 Reflect particle (normally into the same block)
!!       HYDROSTATIC_NVREFL
!!   --------------------------------------------------------------------------------
!!       PERIODIC                   Move particle (to the opposite side of the domain)
!!   --------------------------------------------------------------------------------
!!
!!   All unrecognized boundary conditions result in dropping, as for OUTFLOW.
!!
!!***

subroutine gr_ptApplyBC(dataElement, propCount, blockID, lostElements, moved, negh, bndBox)
  
  use Grid_data,ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,&
                       gr_kmin,gr_kmax
  use gr_ptInterface, ONLY : gr_ptOneFaceBC
  use gr_ptData, ONLY : gr_ptBlk,gr_ptPosx,gr_ptPosy,gr_ptPosz
  use Grid_interface, ONLY : Grid_outsideBoundBox
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer, intent(IN) :: blockID, propCount
  real,dimension(propCount),intent(INOUT)::dataElement
  integer, intent(INOUT) :: lostElements
  logical, intent(INOUT) :: moved
  integer, dimension(MDIM),intent(INOUT) :: negh
  real, dimension(LOW:HIGH,MDIM), intent(IN) :: bndBox

  integer :: blk
  logical :: onBoundary,lostParticle

!!  lostElements=0
  onBoundary=.false.
  !! The particle can only leave through one face, which is why we have the
  !! the else if for figuring out the face. But it can be outside along multiple
  !! dimensions, which is why every dimension must be individually checked.
  !! If the particle was found to be crossing the physical boundary at any time
  !! then depending upon the boundary condition it may or may not end up in the
  !! same block. It may even be lost. That is why if the lostParticle is not true
  !! and boundary condition was applied to the particle then there is a call to
  !! Grid_outsideBoundBox to find if the final position of the particle is outside
  !! the current block
  if(dataElement(gr_ptPosx)<gr_imin) then
     call gr_ptOneFaceBC(dataElement,propCount,IAXIS,LOW,blockID,lostElements)
     onBoundary=.true.
  elseif(dataElement(gr_ptPosx)>gr_imax) then
     call gr_ptOneFaceBC(dataElement,propCount,IAXIS,HIGH,blockID,lostElements)
     onBoundary=.true.
  end if
  blk=dataElement(gr_ptBlk)
  lostParticle=((blk==NONEXISTENT).or.(blk==LOST))
  
#if(NDIM>1)
  if(.not.lostParticle) then
     if(dataElement(gr_ptPosy)<gr_jmin) then
        call gr_ptOneFaceBC(dataElement,propCount,JAXIS,LOW,blockID,lostElements)
        onBoundary=.true.
     elseif(dataElement(gr_ptPosy)>gr_jmax)then
        call gr_ptOneFaceBC(dataElement,propCount,JAXIS,HIGH,blockID,lostElements)
        onBoundary=.true.
     end if
     blk=dataElement(gr_ptBlk)
     lostParticle=((blk==NONEXISTENT).or.(blk==LOST))
  end if

#if(NDIM>2)
  if(.not.lostParticle) then
     if(dataElement(gr_ptPosz)<gr_kmin)then
        call gr_ptOneFaceBC(dataElement,propCount,KAXIS,LOW,blockID,lostElements)
        onBoundary=.true.
     elseif(dataElement(gr_ptPosz)>gr_kmax)then
        call gr_ptOneFaceBC(dataElement,propCount,KAXIS,HIGH,blockID,lostElements)
        onBoundary=.true.
     end if
     blk=dataElement(gr_ptBlk)
     lostParticle=((blk==NONEXISTENT).or.(blk==LOST))
  end if
#endif !! the endif for if(NDIM>2)

#endif !! the endif for if(NDIM>1)
  if(blk==LOST) then
     moved=.false.
  elseif((.not.lostParticle).and.onBoundary) then
     call Grid_outsideBoundBox(dataElement(gr_ptPosx:gr_ptPosz),bndBox,moved,Negh)
  end if
  return
end subroutine gr_ptApplyBC
