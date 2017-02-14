!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/gr_mgPfftData
!!
!! NAME
!!
!!  gr_mgPfftData
!!
!! SYNOPSIS
!!
!!  use gr_mgPfftData
!!
!!  This includes subunit scope data for the PFFT extensions to Multigrid MC.
!!
!!***

module gr_mgPfftData

#include "constants.h"

  implicit none

  real, allocatable, save, dimension(:) ::  gr_mgPfftInArray, gr_mgPfftOutArray, gr_mgPfftTranArray
  real, save :: gr_mgPfftPoisfact
  integer, save :: gr_mgPfftSolveFlag, gr_mgPfftLastMappedLevel
  integer, save :: gr_mgPfftMaxDirectSolveLevel

  !This is like gr_mgBndType, but initialized from runtime parameters
  !that are defined specifically for Multigrid+Pfft. May be redundant.
  !Currently used for consistency checking.
  integer, save, dimension(2*MDIM) :: gr_mgBcTypes
  !This should hold the same information as gr_mgBcTypes, but expressed in
  !symbols from the Grid namespace (defined in constants.h, e.g., PERIODIC),
  !instead of from the Multigrid namespace (defined in Multigrid.h, e.g.,
  !MG_BND_PERIODIC.)
  integer, save, dimension(2*MDIM) :: gr_mgPfftBcTypes

end module gr_mgPfftData
