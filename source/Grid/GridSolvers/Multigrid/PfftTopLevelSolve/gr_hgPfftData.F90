!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgPfftData
!!
!! NAME
!!
!!  gr_hgPfftData
!!
!! SYNOPSIS
!!
!!  use gr_hgPfftData
!!
!!  This includes subunit scope data for the PFFT extensions to Multigrid.
!!
!!***

module gr_hgPfftData

#include "constants.h"

  implicit none

  real, allocatable, save, dimension(:) ::  gr_hgPfftInArray, gr_hgPfftOutArray, gr_hgPfftTranArray
  real, save :: gr_hgPfftPoisfact
  integer, save :: gr_hgPfftSolveFlag, gr_hgPfftLastMappedLevel
  integer, save :: gr_hgPfftMaxDirectSolveLevel

  !This is like gr_hgBndType, but initialized from runtime parameters
  !that are defined specifically for Multigrid+Pfft. May be redundant.
  !Currently used for consistency checking.
  integer, save, dimension(2*MDIM) :: gr_hgBcTypes
  !This should hold the same information as gr_hgBcTypes, but expressed in
  !symbols from the Grid namespace (defined in constants.h, e.g., PERIODIC),
  !instead of from the Multigrid namespace (defined in Multigrid.h, e.g.,
  !MG_BND_PERIODIC.)
  integer, save, dimension(2*MDIM) :: gr_hgPfftBcTypes

end module gr_hgPfftData
