!!****f* source/Grid/Grid_getBlkBC
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
!!  This function finds out if a block face is
!!  on the physical boundary, and if so, returns the boundary condition.
!!
!!  This function finds out if a block face is on the physical
!!  boundary. If the face is on the boundary and the boundary
!!  conditions is not periodic, it returns the boundary condition in
!!  array "faces". If the argument "onBoundary" is present then if
!!  all boundary conditions are not periodic, it returns all values
!!  identical to faces. If one or more boundary conditions are
!!  periodic, then the array "faces" returns NOT_BOUNDARY in corresponding elements,
!!  while the array "onBoundary" returns the value
!!  PERIODIC on corresponding faces.
!!
!!  The boundary conditions are defined in the header file
!!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC.  NOT_BOUNDARY, also
!!  defined in constants.h, is returned if a block face is not on a 
!!  physical boundary.
!!   
!! ARGUMENTS 
!!
!!  blockId - the local blockId 
!!  faces   - array returned holding boundary conditions, except when periodic,
!!            if any of the faces of the block are on a physical boundary. 
!!                
!!            The first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS.
!! onBoundary - array returned with boundary conditions including periodic
!!
!! NOTES
!!
!!   The #define constants LOW, HIGH, IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.
!!   instead of faces(2,3) = PERIODIC the code reads
!!   faces(HIGH, KAXIS) = PERIODIC.
!!
!!***


subroutine Grid_getBlkBC(blockId, faces,onBoundary)
implicit none
#include "constants.h"

  integer, intent(in) :: blockId
  integer, dimension(2,MDIM),intent(out):: faces
  integer, optional, dimension(LOW:HIGH,MDIM),intent(out):: onBoundary

  faces = NOT_BOUNDARY
  
  if (present (onBoundary)) then     
     onBoundary = NOT_BOUNDARY
  endif
 

  return
end subroutine Grid_getBlkBC
