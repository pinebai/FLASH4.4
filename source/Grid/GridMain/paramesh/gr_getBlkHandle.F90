!!****if* source/Grid/GridMain/paramesh/gr_getBlkHandle
!!
!! NAME
!!
!!  gr_getBlkHandle
!!
!! SYNOPSIS
!!
!!  call gr_getBlkHandle(integer(in)    :: remoteBlockID,
!!                       integer(in)    :: pe,
!!                       integer(INOUT) :: blockHandle)
!!
!! DESCRIPTION
!!
!!  Convert a (remoteBlockID, pe) pair to a local blockHandle, if possible.
!!  Only meaningful with PARAMESH 3 ff.
!!
!!  This routine finds where metadata for a remote block with address
!!  (remoteBlockID,pe) is 'cached' in the local arrays that hold
!!  per-block meta information (e.g., parent, child, coords, ...).
!!  If successful, it returns a handle for the information in blockHandle.
!!  This handle can be used as an index (in the blockID position, i.e.,
!!  usually the last index) for those arrays (i.e., parent, child, coords, ...).
!!
!!  If pe is the local PE, the subroutine simply returns remoteBlockID in
!!  blockHandle.
!!
!!  If pe is different from the local PE and no cached information for
!!  (remoteBlockID,pe) is found, blockHandle is left unchanged.  The caller
!!  should deposit a special value, such as -1, in blockHandle before the
!!  call and check the value after return, in order to detect this
!!  situation.
!!
!! ARGUMENTS
!!
!!   remoteBlockID :   Block ID of the block on the PE given by pe.
!!
!!   pe :              The PE on which remoteBlockID is the blocks's local ID.
!!
!!   blockHandle :     On return, changed to contain a handle for the block
!!                     if meta information is available locally;
!!                     unchanged otherwise.
!!
!! SEE ALSO
!!
!!  amr_mpi_find_blk_in_buffer
!!
!!
!!***

subroutine gr_getBlkHandle(remoteBlockID, pe, blockHandle)

#include "Flash.h"

#ifndef FLASH_GRID_PARAMESH2
  use tree, ONLY: laddress          !, strt_buffer, last_buffer
  use mpi_morton, ONLY: ladd_strt, ladd_end
#endif
  use Grid_data, ONLY : gr_meshMe

  implicit none
  integer, intent(in) :: remoteBlockID, pe
  integer, intent(INOUT) :: blockHandle

  integer :: iblk
  logical :: lfound

  if (pe == gr_meshMe) then
     blockHandle = remoteBlockID
  else

     lfound = .false.

#ifndef FLASH_GRID_PARAMESH2
     iblk = ladd_strt(pe)
     do while(.not.lfound.and.iblk.le.ladd_end(pe))
        !         iblk = strt_buffer
        !         do while(.not.lfound.and.iblk.le.last_buffer)
        if( (pe.eq.laddress(2,iblk))  .and. & 
             &         (remoteBlockID.eq.laddress(1,iblk)) ) then
           ! found the corresponding block id iblk
           blockHandle = iblk
           !          tsurrblks(3,i,j,k) = nodetype(iblk)     !?????
           lfound = .true.
#ifdef DEBUG
           if(gr_meshMe.eq.0) & 
                &       write(*,*) 'pe ',gr_meshMe,' looking for ',remoteBlockID,',',pe &
                &             ,' in slot ',iblk,' FOUND ', &
                &             laddress(:,iblk)
#endif /*DEBUG */
        else
#ifdef DEBUG
           if(gr_meshMe.eq.0) & 
                &       write(*,*) 'pe ',gr_meshMe,' looking for ',remoteBlockID,',',pe &
                &             ,' in slot ',iblk,' found ', &
                &             laddress(:,iblk)
#endif /*DEBUG */
           iblk = iblk+1  
        endif
     enddo
     !-
#endif
  endif


end subroutine gr_getBlkHandle
