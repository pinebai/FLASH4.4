!!****f* source/Grid/Grid_setFluxHandling
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
!!   handling - a string indicating the handling requested;
!!              either 'consv_fluxes' or 'consv_flux_densities'.
!!   status - if present, will be set to 0 for success and 1 for failure.
!!
!!
!!***

subroutine Grid_setFluxHandling(handling, status)

  implicit none

  character(len=*),intent(IN) :: handling
  integer,intent(OUT),OPTIONAL :: status

  ! This stub pretends to always succeed:
  if (present(status)) status = 0
  

end subroutine Grid_setFluxHandling
