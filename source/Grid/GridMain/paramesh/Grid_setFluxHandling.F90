!!****if* source/Grid/GridMain/paramesh/Grid_setFluxHandling
!!
!! NAME
!!
!!  Grid_setFluxHandling
!!
!! SYNOPSIS
!!
!!  call Grid_setFluxHandling(character*(*){IN) :: handling,
!!                        OPTIONAL,integer(OUT) :: status)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  handling - A string indicating how fluxes should be handled by
!!             the Grid unit during Grid_conserveFluxes processing.
!!             Should be one of the following values:
!!               'consv_fluxes'
!!               'consv_flux_densities'
!!             These correspond to the PARAMESH 4 logicals of the
!!             same names.
!!
!!  status  -  If this optional argument is present, a nonzero value
!!             on return indicates that the Grid flux handling could
!!             not be set as requested.
!!             If this optional argument is not present, the program
!!             aborts if the requested flux handling cannot be set.
!!
!!
!!***

#include "Flash.h"

subroutine Grid_setFluxHandling(handling, status)
  use Driver_interface, ONLY: Driver_abortFlash

#ifdef FLASH_GRID_PARAMESH2
#define consv_fluxes .FALSE.
#define consv_flux_densities .TRUE.
#else
  use physicaldata, ONLY: consv_fluxes, consv_flux_densities
  use paramesh_dimensions, ONLY: ndim, red_f
#endif

  implicit none

  character(len=*),intent(IN) :: handling
  integer,intent(OUT),OPTIONAL :: status

  if (handling=='consv_fluxes') then
     if (consv_fluxes) then ! state is already as requested...
        if (present(status)) status = 0
        return
     end if
  else if (handling=='consv_flux_densities') then
     if (consv_flux_densities) then ! state is already as requested...
        if (present(status)) status = 0
        return
     end if
  else
     call Driver_abortFlash('Grid_setFluxHandling: unrecognized handling requested: ' // handling)
  end if

#ifdef LIBRARY
#define LIBRARY_OR_PM4DEV
#else
#ifdef FLASH_GRID_PARAMESH4DEV
#define LIBRARY_OR_PM4DEV
#endif
#endif

#ifndef LIBRARY_OR_PM4DEV
  ! State is different from what was requested, but we cannot change it now...
  if (present(status)) then
     status = 1                 !indicates failure
     return
  else
     call Driver_abortFlash('Grid_setFluxHandling: Handling by Grid is compiled in, cannot change at run-time to '&
          // handling)
  end if
#endif


#ifndef FLASH_GRID_PARAMESH2
#ifdef LIBRARY_OR_PM4DEV
  if (handling=='consv_fluxes') then
     consv_fluxes = .TRUE.
     consv_flux_densities = .FALSE.
     if (ndim == 3) then
        red_f = 1.0
     elseif (ndim == 2) then
        red_f = 0.5
     elseif (ndim == 1) then
        red_f = 0.25
     endif
  else if (handling=='consv_flux_densities') then
     consv_fluxes = .FALSE.
     consv_flux_densities = .TRUE.
     red_f = 0.25
  end if
#endif
#endif

  if (present(status)) status = 0
  

end subroutine Grid_setFluxHandling
