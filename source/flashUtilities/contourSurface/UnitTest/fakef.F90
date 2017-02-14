

#include "constants.h"
module Grid_interface

!  interface Grid_getDeltas
!     subroutine Grid_getDeltas(blockId, del)
!       integer, intent(in) :: blockId
!       real, dimension(MDIM), intent(out) :: del
!     end subroutine Grid_getDeltas
!  end interface


  contains

     subroutine Grid_getDeltas(blockId, del)
       integer, intent(in) :: blockId
       real, dimension(MDIM), intent(out) :: del

       del(:) = 1.0
     end subroutine Grid_getDeltas

end module


module Driver_interface

  contains

    subroutine Driver_abortFlash (errorMessage)
      implicit none
      character(len=*), intent(in) :: errorMessage

    print *, errorMessage
    stop
    end subroutine Driver_abortFlash


end module
