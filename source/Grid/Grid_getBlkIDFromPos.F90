!!****f* source/Grid/Grid_getBlkIDFromPos
!!
!! NAME
!!  Grid_getBlkIDFromPos
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkIDFromPos(real(IN)     :: pos(:), 
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID,
!!                   optional,integer(IN)  :: comm)
!!  
!! DESCRIPTION 
!! 
!!  Returns the processor and block ID
!!  containing the cell that overlaps with the 
!!  specified position co-ordinate.
!! 
!! 
!! ARGUMENTS 
!!
!!  pos        :: co-ordinates of the point
!!  ansBlockID    :: the local blockID of the block that contains the point;
!!                   or NONEXISTENT if no matching block was found globally.
!!  ansProcID     :: the ID of the processor that contains the point;
!!                   or NONEXISTENT if no matching block was found anywhere.
!!  comm       :: if communication is necessary, a communicator must be
!!                specified here.
!!                Communication is necessary unless the setup uses BITTREE.
!!                The comm argument is ignored if BITTREE is used.
!!
!! NOTES
!!
!!  If a communicator comm is present and used, then this routine must be
!!  called collectively by all MPI tasks in the communicator. In a collective
!!  call, all tasks are expected to call with the same values for the
!!  coordinate tuple in pos. In other words, The intent(in) arguments should
!!  have the same values on all tasks.
!!
!!  On return, ansBlockID and ansBlockID will contain the same values on all tasks
!!  that provided the same input values for pos (and comm, if appropriate).
!!***

#include "constants.h"

subroutine Grid_getBlkIDFromPos(pos,ansBlockID, ansProcID,comm)

  implicit none

  real, dimension(1:MDIM), intent(IN) :: pos
  integer, intent(OUT) :: ansBlockID, ansProcID
  integer, optional, intent(IN) :: comm
  
  ansProcID=NONEXISTENT
  ansBlockID=NONEXISTENT
  
end subroutine Grid_getBlkIDFromPos
