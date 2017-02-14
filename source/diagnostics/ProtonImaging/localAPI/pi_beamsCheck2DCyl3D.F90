!!****if* source/diagnostics/ProtonImaging/localAPI/pi_beamsCheck2DCyl3D
!!
!! NAME
!!
!!  pi_beamsCheck2DCyl3D
!!
!! SYNOPSIS
!!
!!  call pi_beamsCheck2DCyl3D ()
!!
!! DESCRIPTION
!!
!!  Checks the collected 3D cartesian beams data for 2D cylindrical grids (cartesian). This is
!!  a rather specialized routine, in which the dimensions where the beams are defined is different
!!  from the dimension of the underlaying grid. All checks which depend on domain grid details
!!  should go in here. Currently it contains the following:
!!
!!         1) Check, if all beam circular lens areas are completely outside the 3D cylindrical
!!            domain, defined by the 2D cylindrical grid.
!!
!!  The 3D cylindrical domain is obtained from the (x=R,y=z) 2D cylindrical grid by a 360
!!  degree revolution around the 2D cylindrical z-axis.
!!
!!  Since the target area of the proton beam is merely used to construct the direction of
!!  the individual protons, there is no need to check if this area is properly located wrt
!!  to the 3D cylindrical domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_beamsCheck2DCyl3D ()
  
  implicit none

  return
end subroutine pi_beamsCheck2DCyl3D
