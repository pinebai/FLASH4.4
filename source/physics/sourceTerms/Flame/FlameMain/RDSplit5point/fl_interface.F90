! declarations of public interface for private Flame subroutines
!
! Aaron Jackson 2010
!
! The basic function of each routine is also described here.


Module fl_interface
#include "constants.h"
#include "Flash.h"

  implicit none

  interface fl_laplacian
     subroutine fl_laplacian(lapl, flam, h, bid)
        implicit none
        real, dimension(:,:,:), intent(out) :: lapl
        real, dimension(:,:,:), intent(in) :: flam
        integer, intent(in) :: bid, h
        !  Calculate the laplacian of the flame progress variable.
        !  The block id (bid) is passed in so that we can retrieve
        !  coordinate info. The laplacian is only computed for the
        !  interior cells, although the indices of lapl also run
        !  over the guard cells to simplify indexing. h is the step size
        !  used to calculate the laplacian
     end subroutine fl_laplacian
  end interface fl_laplacian

end Module fl_interface
