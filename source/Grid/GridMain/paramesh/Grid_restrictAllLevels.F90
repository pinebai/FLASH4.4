!!****if* source/Grid/GridMain/paramesh/Grid_restrictAllLevels
!!
!! NAME
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!! 
!!  Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid data to all refinement levels. Normally FLASH
!!  only evolves on the leaf blocks, calling this routine makes all
!!  levels have valid data.  This is mostly for visualization purposes to
!!  be able to look at different levels of resolution
!!  
!!  
!!
!!***


subroutine Grid_restrictAllLevels()

  use Timers_interface, ONLY: Timers_start, Timers_stop

  implicit none
  call Timers_start("restrictAll")
  call gr_restrictTree()
  call Timers_stop("restrictAll")
  
end subroutine Grid_restrictAllLevels
