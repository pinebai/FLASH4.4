!!****f* source/Grid/Grid_getGeometry
!!
!! NAME
!!  Grid_getGeometry
!!
!! SYNOPSIS
!!
!!  Grid_getGeometry(integer (OUT)  :: geometry)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the global grid geometry.
!!  valid values are (CARTESIAN, POLAT, CYLINDRICAL, SPHERICAL) defined
!!  in file "constants.h"
!!
!!
!! ARGUMENTS
!!
!!  geometry - returned value
!!
!!***

subroutine Grid_getGeometry(geometry)

  implicit none

  real, intent(OUT) :: geometry

  geometry = 0

end subroutine Grid_getGeometry
