!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsCheck2DCyl3D
!!
!! NAME
!!
!!  ed_beamsCheck2DCyl3D
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck2DCyl3D ()
!!
!! DESCRIPTION
!!
!!  Checks the collected 3D cartesian beams data for 2D cylindrical grids (cartesian). This is
!!  a rather specialized routine, in which the dimensions where the beams are defined is different
!!  from the dimension of the underlaying grid. All checks which depend on domain grid details
!!  should go in here. Currently it contains the following:
!!
!!     1) Check, if all beam elliptical 3D target areas are completely within the
!!        3D cylindrical domain implicated by the 2D cylindrical grid.
!!
!!     2) Check, if all beam elliptical 3D lens areas are completely outside the
!!        3D cylindrical domain implicated by the 2D cylindrical grid.
!!
!!  The 3D cylindrical domain is obtained from the (x=R,z) 2D cylindrical grid by a 360 degree
!!  revolution around the 2D cylindrical z-axis.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck2DCyl3D ()

  implicit none

  return
end subroutine ed_beamsCheck2DCyl3D
