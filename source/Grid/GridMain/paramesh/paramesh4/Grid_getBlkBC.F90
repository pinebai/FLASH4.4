!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_getBlkBC
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
!!  This function finds out if a block face is on the physical
!!  boundary. If the face is on the boundary and the boundary
!!  conditions is not periodic, it returns the boundary condition in
!!  array "faces". If the argument "onBoundary" is present then if
!!  all boundary conditions are not periodic, it returns all values
!!  identical to faces. If one or more boundary conditions are
!!  periodic, then the array "faces" returns NOT_BOUNDARY,
!!  while the array "onBoundary" returns the value
!!  PERIODIC on corresponding faces.
!!
!!  The boundary conditions are defined in the header file
!!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC.  NOT_BOUNDARY, also
!!  defined in constants.h, is returned if a block face is not on a 
!!  physical boundary.
!!  
!!   
!! ARGUMENTS 
!!
!!  blockId - the local blockId 
!!  faces   - array returned holding boundary conditions (except PERIODIC),
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


subroutine Grid_getBlkBC(blockid, faces, onBoundary)

  use tree, ONLY : neigh,bnd_box
  use Grid_data,ONLY : gr_domainBC,gr_globalDomain,gr_minCellSizes
  use gr_interface, ONLY : gr_extractBCForDirection
implicit none
#include "constants.h"
  integer, intent(in) :: blockid
  integer, dimension(LOW:HIGH,MDIM),intent(out):: faces
  integer, optional, dimension(LOW:HIGH,MDIM),intent(out):: onBoundary
  integer :: idirf
  integer :: axis,face

#ifdef DEBUG_GRID
  do idirf=1,6
     if (neigh(1,idirf,blockid) <= -1000) then
        print*,'Grid_getBlkBC: neigh(1,',idirf,',',blockid,') is ',neigh(1,idirf,blockid)
     end if
  end do
#endif

  faces = NOT_BOUNDARY
  if(neigh(1,1,blockid) <= PARAMESH_PHYSICAL_BOUNDARY .AND. neigh(1,1,blockid) .NE. PERIODIC) &
       faces(LOW,IAXIS) = gr_extractBCForDirection(neigh(1,1,blockid), IAXIS,HIGH)
  if(neigh(1,2,blockid) <= PARAMESH_PHYSICAL_BOUNDARY .AND. neigh(1,2,blockid) .NE. PERIODIC) &
       faces(HIGH,IAXIS) = gr_extractBCForDirection(neigh(1,2,blockid), IAXIS,LOW)
  if(neigh(1,3,blockid) <= PARAMESH_PHYSICAL_BOUNDARY .AND. neigh(1,3,blockid) .NE. PERIODIC) &
       faces(LOW,JAXIS) = gr_extractBCForDirection(neigh(1,3,blockid), JAXIS,HIGH)
  if(neigh(1,4,blockid) <= PARAMESH_PHYSICAL_BOUNDARY .AND. neigh(1,4,blockid) .NE. PERIODIC) &
       faces(HIGH,JAXIS) = gr_extractBCForDirection(neigh(1,4,blockid), JAXIS,LOW)
  if(neigh(1,5,blockid) <= PARAMESH_PHYSICAL_BOUNDARY .AND. neigh(1,5,blockid) .NE. PERIODIC) &
       faces(LOW,KAXIS) = gr_extractBCForDirection(neigh(1,5,blockid), KAXIS,HIGH)
  if(neigh(1,6,blockid) <= PARAMESH_PHYSICAL_BOUNDARY .AND. neigh(1,6,blockid) .NE. PERIODIC) &
       faces(HIGH,KAXIS) = gr_extractBCForDirection(neigh(1,6,blockid), KAXIS,LOW)

!!  This is addition to help with MapParticlesToMesh routines. 
!!  It only works for outer boundaries, not for obstacles etc, not even a step.
!!  Idea is to return PERIODIC for the periodic BC instead of the neighbors
!!           as faces array returns.

  if(present(onBoundary)) then
     onBoundary=faces
     do axis = 1,MDIM
        do face = LOW,HIGH
           if (gr_domainBC(face,axis)==PERIODIC) then
              if (nearlySame(bnd_box(face,axis,blockid), gr_globalDomain(face,axis),  &
                             gr_minCellSizes(axis) )) &
                                                     onBoundary(face,axis)=PERIODIC
           end if
        end do
     end do
  end if

  return

contains
  logical function nearlySame(x,y,del)
    real,intent(IN) :: x,y,del

    nearlySame = (abs(x-y).LE. 0.01 * del)
  end function nearlySame

end subroutine Grid_getBlkBC





