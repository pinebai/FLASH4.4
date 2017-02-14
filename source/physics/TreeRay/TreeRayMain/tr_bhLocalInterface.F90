!!****ih* source/physics/TreeRay/TreeRayMain/tr_bhLocalInterface
!!
!! This is the header file for the TreeRay module that defines its
!! public interfaces.
!!***

Module tr_bhLocalInterface

#include "constants.h"
#include "Flash.h"

  interface tr_bhFinalizeCell
    subroutine tr_bhFinalizeCell(solnPoint, dl_poc, eflux, cdMaps)
      use TreeRay_data, ONLY : tr_nEb, tr_nPix, tr_nCd
      implicit none
      real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
      real, DIMENSION(tr_nEb, 0:tr_nPix-1), intent(IN) :: eflux
      real :: cdMaps(tr_nCd, 0:tr_nPix-1)
      real, DIMENSION(MDIM), intent(IN) :: dl_poc
    end subroutine tr_bhFinalizeCell
  end interface

  interface tr_bhRadToGas
    subroutine tr_bhRadToGas(solnPoint, vol_poc, area_poc, dt)
      real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
      real, intent(IN) :: vol_poc, area_poc, dt
    end subroutine tr_bhRadToGas
  end interface

  interface tr_bhGenIntersectList
    subroutine tr_bhGenIntersectList()
    end subroutine tr_bhGenIntersectList
  end interface

  interface tr_checkQuan
    subroutine tr_checkQuan(ncall)
      integer, intent(IN) :: ncall
    end subroutine tr_checkQuan
  end interface

end Module tr_bhLocalInterface

