!!****if* source/Grid/GridMain/Chombo/Grid_getBlkBC
!!
!! NAME
!!  Grid_getBlkBC
!!
!! SYNOPSIS
!!
!!  Grid_getBlkBC(integer(IN)  :: blockId,
!!                integer(OUT) :: faces(2,MDIM),
!!       optional,integer(OUT) :: onBoundary(2,MDIM))
!!                    
!! DESCRIPTION 
!!  Returns the boundary condition for each face of the block.
!!  
!!  this function finds out if a block face is
!!  on the physical boundary, and if so, returns the boundary condition.
!!  The boundary conditions are defined in the header file
!!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC.  NOT_BOUNDARY, also
!!  defined in constants.h is returned if a block face is not on a 
!!  physical boundary.
!!   
!! ARGUMENTS 
!!
!!  blockId - the local blockId 
!!  faces   - array returned holding boundary conditions
!!                
!!            the first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS
!! onBoundary - array returned with boundary conditions including periodic
!!
!! NOTES
!!
!!   The #define constants LOW, HIGH, IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.       
!!   instead of faces(2,3) = PERIODIC the code reads
!!   faces(HIGH, KAXIS) = PERIODIC
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkBC(blockId, faces, onBoundary)
  use Grid_data, ONLY : gr_globalDomain, gr_domainBC
  use Grid_interface, ONLY : Grid_getBlkBoundBox
  implicit none
  integer, intent(in) :: blockId
  integer, dimension(2,MDIM),intent(out):: faces
  integer, optional, dimension(LOW:HIGH,MDIM),intent(out):: onBoundary
  real, dimension(2,MDIM) :: bnd_box
  integer :: axis, face

  !NOTE: Does not behave the same as PARAMESH.  If we are
  !on a periodic boundary then this version will return "periodic"
  !in faces array.
  call Grid_getBlkBoundBox(blockId,bnd_box)

  do axis = 1,MDIM
     do face = LOW,HIGH
        if ( bnd_box(face,axis) == gr_globalDomain(face,axis) ) then
           faces(face,axis) = gr_domainBC(face,axis)
        else
           faces(face,axis) = NOT_BOUNDARY
        end if
     end do
  end do
end subroutine Grid_getBlkBC
