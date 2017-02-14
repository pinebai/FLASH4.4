!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_finalize
!!
!! NAME
!!
!!  Hydro_finalize
!!
!! SYNOPSIS
!!
!!  Hydro_finalize()
!!
!! DESCRIPTION
!!
!!  Deallocates any memory that has been allocated in the Hydro Unit
!!  and prepares the unit for shutdown
!!
!!
!!***


subroutine Hydro_finalize()

  use Hydro_data, ONLY : hy_xarea,hy_xdtdx, hy_xgrav,hy_xngrav,hy_xfict
  use Hydro_data, ONLY : hy_yarea,hy_ydtdy, hy_ygrav,hy_yngrav,hy_yfict
  use Hydro_data, ONLY : hy_zarea,hy_zdtdz, hy_zgrav,hy_zngrav,hy_zfict
  use Hydro_data, ONLY : hy_fluxCorrect
  implicit none
#include "Flash.h"
#ifndef FIXEDBLOCKSIZE
  if(hy_fluxCorrect) then
!!     deallocate(hy_xarea)
!!     deallocate(hy_xngrav)
!!    deallocate(hy_xgrav)
!!     deallocate(hy_xdtdx)
!!     deallocate(hy_xfict)
!!     deallocate(hy_yarea)
!!     deallocate(hy_yngrav)
!!     deallocate(hy_ygrav)
!!     deallocate(hy_ydtdy)
!!     deallocate(hy_yfict)
!!     deallocate(hy_zarea)
!!     deallocate(hy_zngrav)
!!     deallocate(hy_zgrav)
!!     deallocate(hy_zdtdz)
!!     deallocate(hy_zfict)
  end if
#endif
end subroutine Hydro_finalize
