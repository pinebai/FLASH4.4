!!****if* source/Simulation/SimulationMain/unitTest/Poisson/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!
!! DESCRIPTION
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!!  This routine is normally called by the implementation of
!!  Grid_updateRefinement. It may also get called repeatedly
!!  during the initial construction of the Grid from
!!  Grid_initDomain.
!!
!! ARGUMENTS
!!
!!  none
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_initDomain
!!  gr_expandDomain
!!
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
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
  use Simulation_data, ONLY : sim_xMin, sim_xMax, sim_yMin, sim_yMax, &
                              sim_zMin, sim_zMax



#include "constants.h"

  implicit none
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

!!$  do l = 1,gr_numRefineVars
!!$     iref = gr_refine_var(l)
!!$     ref_cut = gr_refine_cutoff(l)
!!$     deref_cut = gr_derefine_cutoff(l)
!!$     ref_filter = gr_refine_filter(l)
!!$     call gr_markRefineDerefine(MyPE,iref,ref_cut,deref_cut,ref_filter)
!!$  end do

#ifdef SPECIAL_REFINEMENT

  select case(NDIM)
  case(2)
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  sim_xMin + 1./4.*(sim_xMax - sim_xMin) +.0005
  specs(2) =  sim_xMin + 3./4.*(sim_xMax - sim_xMin) -.0005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  sim_yMin + 1./4.*(sim_yMax - sim_yMin) +.0005
  specs(4) =  sim_yMin + 3./4.*(sim_yMax - sim_yMin) -.0005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0. !sim_zMin + 1./4.*(sim_zMax - sim_zMin) +.05
  specs(6) =  0. !sim_zMin + 3./4.*(sim_zMax - sim_zMin) -.05
  !! End of coordinate information -------------------------------

  case(3)
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction
  specs(1) =  sim_xMin + 1./3.*(sim_xMax - sim_xMin) +.0005
  specs(2) =  sim_xMin + 2./3.*(sim_xMax - sim_xMin) -.0005

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  sim_yMin + 1./3.*(sim_yMax - sim_yMin) +.0005
  specs(4) =  sim_yMin + 2./3.*(sim_yMax - sim_yMin) -.0005

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  sim_zMin + 0./4.*(sim_zMax - sim_zMin) +.0005
  specs(6) =  sim_zMin + 4./4.*(sim_zMax - sim_zMin) -.0005
  !! End of coordinate information -------------------------------

  end select

  !! Decide wether or not we refine only blocks completely
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 0.0

  !! Bring all qualifying blocks to this level of refinement
  lref = lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)

#endif

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif
  return
#endif

end subroutine Grid_markRefineDerefine
