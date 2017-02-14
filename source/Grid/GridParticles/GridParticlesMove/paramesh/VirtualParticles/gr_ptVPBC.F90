!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/VirtualParticles/gr_ptVPBC
!!
!! NAME
!!
!!  gr_ptVPBC
!!
!! SYNOPSIS
!!
!!  gr_ptVPBC(real(INOUT)     :: dataElement(propCount),
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

subroutine gr_ptVPBC(dataElement, propCount, leftDomain, onBoundary)
  
  use Grid_data,ONLY : gr_globalDomain, gr_domainBC
  use gr_ptInterface, ONLY : gr_ptOneFaceBC
  use gr_ptData, ONLY : gr_ptBlk,gr_ptPosx,gr_ptPosy,gr_ptPosz
  use gr_ptVPData, only : gr_ptVPBufferFactor
  use Grid_interface, ONLY : Grid_getDeltas


  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer, intent(IN) ::  propCount
  real,dimension(propCount),intent(INOUT)::dataElement
  logical, intent(OUT) :: leftDomain
  integer, dimension(MDIM),intent(OUT) :: onBoundary


  integer :: blk, i, lostElements
  logical :: done, left, moved
  real,dimension(MDIM) :: pos, deltas

!!  lostElements=0


  pos(IAXIS:KAXIS)=dataElement(gr_ptPosx:gr_ptPosz)
  blk=int(dataElement(gr_ptBlk))
  onBoundary=NOBOUNDARY
  left=.false.
  call Grid_getDeltas(blk,deltas)
  do i = 1,NDIM
     if(.not.left) then
        done=.false.
        if(pos(i)<gr_globalDomain(LOW,i))then

           lostElements=0
           !! if the particle is outside the lower face of this dimension
           !! then applying boundary condition may either cause the particle
           !! to leave the domain, reflect back into the same block or
           !! go into the adjacent block for periodic boundary
           call gr_ptOneFaceBC(dataElement,propCount,i,LOW,blk,lostElements)
           !! if the particle left the domain then lostElements will be 1
           left=(lostElements>0)
           !! if the boundary conditions are periodic then the new positions
           !! of the particle will be in a block on the opposite side, so the
           !! comparison with global domain will fail. Not convinced yet that
           !! this following copy into pos is necessary
           if(gr_domainBC(LOW,i)/=PERIODIC)then
              pos(IAXIS:KAXIS)=dataElement(gr_ptPosx:gr_ptPosz)
           end if
           done=.true.
        end if
        if(.not.left) then
           !! if the particle is still in domain, and is within smearing length
           !! of the boundary then it is on boundary at the left edge
           if((pos(i)-gr_globalDomain(LOW,i))<gr_ptVPBufferFactor*deltas(i))onBoundary(i)=LEFT_EDGE

        end if
        !! do similar processing for the upper face
        if((.not.done).and.(.not.left))then
           if(pos(i)>gr_globalDomain(HIGH,i))then

              lostElements=0
              call gr_ptOneFaceBC(dataElement,propCount,i,HIGH,blk,lostElements)
              left=(lostElements>0)
              if(gr_domainBC(HIGH,i)/=PERIODIC)then
                 pos(IAXIS:KAXIS)=dataElement(gr_ptPosx:gr_ptPosz)
              end if
           end if
           if(.not.left)then
              if((gr_globalDomain(HIGH,i)-pos(i))<gr_ptVPBufferFactor*deltas(i))onBoundary(i)=RIGHT_EDGE
           end if
        end if
     end if
  end do
  leftDomain=left
  !print*,'before returning ',leftDomain
  return
end subroutine gr_ptVPBC
