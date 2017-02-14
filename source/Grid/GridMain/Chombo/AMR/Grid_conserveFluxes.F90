!!****if* source/Grid/GridMain/Chombo/AMR/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  Grid_conserveFluxes(integer(IN) :: axis,
!!                      integer(IN) :: level)
!!  
!! DESCRIPTION 
!!  
!!  Flux conservation is necessary when 2 blocks of differing
!!  levels (meaning having different grid spacings) border 
!!  one another. 
!!  
!!  This routine can perform flux conservation on the finest
!!  blocks, the most typical usage for the Paramesh Grid or on
!!  blocks of a certain level.
!!  
!!  The routine overwrites the flux arrays maintained by the Grid
!!  
!! ARGUMENTS 
!!
!!  axis - conserve fluxes in just one direction if 
!!         IAXIS, JAXIS, KAXIS, or in all directions if ALLDIR.
!!         These constants are defined in constants.h.
!!
!!  level - refinement level. Ignored.
!!
!!***
subroutine Grid_conserveFluxes(axis, level)
  use chombo_f_c_interface, ONLY : ch_reflux
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(in) :: axis, level

  call ch_reflux()
end subroutine Grid_conserveFluxes
