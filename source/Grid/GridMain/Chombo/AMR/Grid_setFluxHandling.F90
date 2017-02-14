!!****if* source/Grid/GridMain/Chombo/AMR/Grid_setFluxHandling
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
  use Grid_data, ONLY : gr_scaleFineFluxes
  implicit none

  character(len=*),intent(IN) :: handling
  integer,intent(OUT),OPTIONAL :: status

  if (handling=='consv_fluxes') then
     gr_scaleFineFluxes = .false.
  else if (handling=='consv_flux_densities') then
     gr_scaleFineFluxes = .true.
  end if

  if (present(status)) status = 0  

end subroutine Grid_setFluxHandling
