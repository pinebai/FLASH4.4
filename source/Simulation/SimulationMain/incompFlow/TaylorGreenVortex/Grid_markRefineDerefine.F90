!!****if* source/Simulation/SimulationMain/incompFlow/TaylorGreenVortex/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_myPE or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var
  use tree, ONLY : newchild, refine, derefine, stay,lrefine_max
  use Grid_interface, ONLY : Grid_markRefineSpecialized,Grid_fillGuardCells
  use Simulation_data


  implicit none
#include "constants.h"
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref

  logical :: gcMask(NUNK_VARS)

  !! Special refinement criteria -----------------
  real, dimension(7) :: specs
  integer :: lref,specsSize
  !! End of special refinement treatment ---------

  call Grid_fillGuardCells(CENTER,ALLDIR)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  sim_xMin + 2./4.*(sim_xMax - sim_xMin) +.005
  specs(2) =  sim_xMax -.005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  sim_yMin +.005
  specs(4) =  sim_yMin + 4./4.*(sim_yMax - sim_yMin) -.005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.005
  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.005
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely 
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !! Bring all qualifying blocks to this level of refinement
  lref = lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif
  return
#endif

end subroutine Grid_markRefineDerefine
