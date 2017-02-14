!!****f* source/IO/IO_checkForPlot
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
  implicit none
  
  logical, intent(out) :: wrotePlot

  wrotePlot = .false.

end subroutine IO_checkForPlot
