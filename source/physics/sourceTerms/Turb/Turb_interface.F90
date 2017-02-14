! declarations of public interface for private Turb subroutines
!
! Aaron Jackson 2010
!
! The basic function of each routine is also described here.


Module Turb_interface
#include "constants.h"
#include "Flash.h"

  implicit none

  interface Turb_laplacian
     subroutine Turb_laplacian(lapl, scalar, h, bid)
        implicit none
        real, dimension(:,:,:), intent(out) :: lapl
        real, dimension(:,:,:), intent(in) :: scalar
        integer, intent(in) :: bid, h
        !  Calculate the laplacian of the scalar field "scalar".
        !  The block id (bid) is passed in so that we can retrieve
        !  coordinate info. The laplacian is only computed for the
        !  interior cells, although the indices of lapl also run
        !  over the guard cells to simplify indexing. h is the step size
        !  used to calculate the laplacian
     end subroutine Turb_laplacian
  end interface Turb_laplacian

  interface Turb_curlMag
     subroutine Turb_curlMag(curl, velX, velY, velZ, h, bid)
        implicit none
        real, dimension(:,:,:), intent(out) :: curl
        real, dimension(:,:,:), intent(in) :: velX, velY, velZ
        integer, intent(in) :: bid, h
        !  Calculate the curl of the velocity field.
        !  The block id (bid) is passed in so that we can retrieve
        !  coordinate info. The curl is only computed for the 
        !  interior cells, although the indices of curl also run
        !  over the guard cells to simplify indexing.
        !  For this implementation, we only need the magnitude of
        !  the curl, so we store this information as a scalar
        !  instead of a vector in curl. h is the step size used
        !  to calculate the curl
     end subroutine Turb_curlMag
  end interface Turb_curlMag

  interface Turb_calc
     subroutine Turb_calc(num_blocks, blockList)
        implicit none
        integer, intent(in)                        :: num_blocks
        integer, intent(in), dimension(num_blocks) :: blockList
        ! Calculate the laplacian of the velocity field and store
        ! in TURB_VAR, LAPY_VAR, and LAPZ_VAR for each block, then
        ! fill guard cells.
        ! Calculate the curl of the laplacian for each block, then
        ! fill guard cells again.
        ! This routine is OP2 in Colin et al. (2000)
     end subroutine Turb_calc
  end interface Turb_calc

  interface Turb_calcCompLimits
     subroutine Turb_calcCompLimits(blkLimits, compLimits, local_stepSize)
        implicit none
        integer, intent(in), dimension(LOW:HIGH,MDIM) :: blkLimits
        integer, intent(out), dimension(LOW:HIGH,MDIM) :: compLimits
        integer, intent(in) :: local_stepSize
        ! Calculate the indices over which the first operator should be 
        ! performed. If turb_stepSize is small enough compared to the number
        ! of guard cells and the type of operator being performed, then an
        ! intermediate guard cell fill is not necessary and we should
        ! calculate the first operator into the guard cells
     end subroutine Turb_calcCompLimits
  end interface Turb_calcCompLimits

  interface Turb_getFilterScale
     subroutine Turb_getFilterScale(dx, de)
        implicit none
        real, intent(in) :: dx
        real, intent(out) :: de
        ! The current implementation does not use dx, but rather the
        ! minimum dx in the simulation. This assumes that the highest
        ! refinement is used to resolve the flame.
        ! The filter scale associated with the turbulence operator is
        ! returned in de.
     end subroutine Turb_getFilterScale
  end interface Turb_getFilterScale

  interface
    subroutine Turb_finalize()
    end subroutine Turb_finalize
  end interface

  interface
    subroutine Turb_init()
    end subroutine Turb_init
  end interface

end Module Turb_interface
