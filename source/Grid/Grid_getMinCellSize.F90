!!****f* source/Grid/Grid_getMinCellSize
!!
!! NAME
!!  Grid_getMinCellSize
!!
!! SYNOPSIS
!!
!!  Grid_getMinCellSize(real (OUT)  :: minCellSize)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the smallest possible cell size in a simulation in any dimension
!!  that does not represent an angle in curvilinear coordinates.
!!
!!
!!
!! ARGUMENTS
!!
!!  minCellSize - returned value
!!
!!***

subroutine Grid_getMinCellSize(minCellSize)

  implicit none

  real, intent(OUT) :: minCellSize

  minCellSize = 1.0

end subroutine Grid_getMinCellSize
