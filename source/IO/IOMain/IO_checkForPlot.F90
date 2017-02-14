!!****if* source/IO/IOMain/IO_checkForPlot
!!
!! NAME
!!    IO_checkForPlot
!!
!! SYNOPSIS
!!
!!    IO_checkForPlot(logical(out) : wrotePlot)
!!
!! DESCRIPTION
!!
!!   This routine sets the argument to true when called on a cycle
!!   following a plot.
!!
!! ARGUMENTS
!!
!!  wrotePlot : set to true if plot written, otherwise false
!!
!!***
subroutine IO_checkForPlot(wrotePlot)
  use IO_data, ONLY: io_wrotePlot
  implicit none
  
  logical, intent(out) :: wrotePlot

  wrotePlot = io_wrotePlot

end subroutine IO_checkForPlot
