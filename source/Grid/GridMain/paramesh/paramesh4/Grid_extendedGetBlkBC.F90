!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_extendedGetBlkBC
!!
!! NAME
!!  Grid_extendedGetBlkBC
!!
!! SYNOPSIS
!!
!!  Grid_extendedGetBlkBC(integer(IN)  :: blockId,
!!                        integer(OUT) :: faces(2,MDIM))
!!                    
!! DESCRIPTION 
!!  Returns the boundary condition for each face of the block.
!!  
!!  this function finds out if a block face is
!!  on the physical boundary, and if so, returns the boundary condition.
!!  The boundary conditions are defined in the header file
!!  constants.h, e.g., OUTFLOW, REFLECTING, PERIODIC.  NOT_BOUNDARY, also
!!  defined in constants.h, is returned if a block face is not on a 
!!  physical boundary.
!!   
!! ARGUMENTS 
!!
!!  blockId - the local blockId 
!!  faces   - array returned holding boundary conditions
!!                
!!            the first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS
!!
!! NOTES
!!
!!   The values returned by this interface differs from those
!!   returned by Grid_getBlkBC only at PERIODIC domain boundaries.
!!   Grid_extendedGetBlkBC returns PERIODIC in the components of faces that
!!   correspond to directions where the currentl block abuts a PERIODIC
!!   domain boundary.  Grid_getBlkBC returns NOT_BOUNDARY instead.
!!
!!   The #define constants LOW, HIGH, IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.       
!!   instead of faces(2,3) = REFLECTING the code reads
!!   faces(HIGH, KAXIS) = REFLECTING
!!
!! SEE ALSO
!!
!!   Grid_getBlkBC
!!***




subroutine Grid_extendedGetBlkBC(blockid, faces)

  use tree, ONLY : neigh, coord, laddress, strt_buffer, last_buffer
  use paramesh_dimensions, ONLY : ndim
  use Grid_data,ONLY : gr_meshMe,gr_anyPeriodic
  use gr_interface,ONLY : gr_extractBCForDirection
  use Driver_interface,ONLY : Driver_abortFlash
implicit none
#include "constants.h"
  integer, intent(in) :: blockid
  integer, dimension(2,MDIM),intent(out):: faces
  integer :: idirf,idir,iblk
  logical :: lfound
  integer :: pe,neighborBlockHandle
  real,dimension(MDIM) :: neighborCenter

#ifdef DEBUG_GRID
  do idirf=1,6
     if (neigh(1,idirf,blockid) <= -1000) then
        print*,'Grid_extendedGetBlkBC: neigh(1,',idirf,',',blockid,') is ',neigh(1,idirf,blockid)
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

99 format (a,1x,i1,2(1x,i4),' ->',3(1x,i4))
  do idir=1,NDIM
     if (gr_anyPeriodic(idir)) then
        if (neigh(1,2*idir-1,blockid) > 0) then
           pe = neigh(2,2*idir-1,blockid)
           if (pe == gr_meshMe) then
              neighborBlockHandle = neigh(1,2*idir-1,blockid)
           else
!!              print 99,'Grid_extendedGetBlkBC: remote LOW neighbor!',2*idir-1,blockid,gr_meshMe,neigh(:,2*idir-1,blockid)
              ! The following logic copied from Tests/amr_1blk_bcset.F90 in Paramesh4 distribution. - KW
              ! This is heavily dependent on Paramesh internals.
              neighborBlockHandle = neigh(1,2*idir-1,blockid)
              lfound = .false.
              iblk = strt_buffer
!!              print*,'laddr1,2:',strt_buffer,':',last_buffer
!!              print*,'laddr1:',laddress(1,strt_buffer:last_buffer)
!!              print*,'laddr2:',laddress(2,strt_buffer:last_buffer)
              do while(.not.lfound.and.iblk.le.last_buffer)
                 if(neighborBlockHandle.eq.laddress(1,iblk).and.pe.eq.laddress(2,iblk) ) then
                    neighborBlockHandle = iblk
                    lfound = .true.
                 endif
                 iblk = iblk+1
              enddo
              if(.not.lfound) then
                 call Driver_abortFlash('Paramesh error: Grid_extendedGetBlkBC: '// & 
                      & 'Low neighbor remote block is not in list received on this pe')
              endif
           end if
!!           print*,'Lo Nei:',neighborBlockHandle,neighborCenter
           neighborCenter(:) = coord(:,neighborBlockHandle)
           if (neighborCenter(idir) > coord(idir,blockid)) faces(LOW,idir) = PERIODIC
        else if (neigh(1,2*idir,blockid) == -1) then
           print*,'Grid_extendedGetBlkBC: LOW neighbor is -1!',2*idir,gr_meshMe,neigh(:,2*idir-1,blockid)
           call Driver_abortFlash('Grid_extendedGetBlkBC: '// & 
                & 'Low neighbor blockid is -1')
        end if

        if (neigh(1,2*idir,blockid) > 0) then
           pe = neigh(2,2*idir,blockid) 
           if (pe == gr_meshMe) then
              neighborBlockHandle = neigh(1,2*idir,blockid)
           else
!!              print 99,'Grid_extendedGetBlkBC: remote HIGH neighbor!',2*idir,blockid,gr_meshMe,neigh(:,2*idir,blockid)
              neighborBlockHandle = neigh(1,2*idir,blockid)
              lfound = .false.
              iblk = strt_buffer
              do while(.not.lfound.and.iblk.le.last_buffer)
                 if(neighborBlockHandle.eq.laddress(1,iblk).and.pe.eq.laddress(2,iblk) ) then
                    neighborBlockHandle = iblk
                    lfound = .true.
                 endif
                 iblk = iblk+1
              enddo
              if(.not.lfound) then
                 call Driver_abortFlash('Paramesh error: Grid_extendedGetBlkBC: '// & 
                      & 'High neighbor remote block is not in list received on this pe')
              endif
           end if
           neighborCenter(:) = coord(:,neighborBlockHandle)
           if (neighborCenter(idir) < coord(idir,blockid)) faces(HIGH,idir) = PERIODIC
        else if (neigh(1,2*idir,blockid) == -1) then
           print*,'Grid_extendedGetBlkBC: HIGH neighbor is -1!',2*idir,gr_meshMe,neigh(:,2*idir,blockid)
           call Driver_abortFlash('Grid_extendedGetBlkBC: '// & 
                & 'High neighbor blockid is -1')
        end if
     end if
  end do

  return
end subroutine Grid_extendedGetBlkBC





