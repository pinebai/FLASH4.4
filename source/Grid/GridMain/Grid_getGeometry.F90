!!****if* source/Grid/GridMain/Grid_getGeometry
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
!!  valid values are (CARTESIAN, POLAR, CYLINDRICAL, SPHERICAL) defined
!!  in file "constants.h"
!!
!!
!! ARGUMENTS
!!
!!  geometry - returned value
!!
!!***

subroutine Grid_getGeometry(geometry)

  use Grid_data, ONLY : gr_geometry

  implicit none

  integer, intent(OUT) :: geometry
  geometry = gr_geometry

end subroutine Grid_getGeometry
