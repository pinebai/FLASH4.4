!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/Grid_markRefineDerefine
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
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_markRefineDerefine()

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
       gr_refine_filter,&
       gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_max
  use Grid_interface, ONLY : Grid_markRefineSpecialized, Grid_fillGuardCells
#endif

  implicit none


  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask


  !! Special refinement criteria -----------------
  real, dimension(7) :: specs
  integer :: lref,specsSize
  !! End of special refinement treatment ---------



  !If this is not an adaptive mesh simulation we simply
  !return from this subroutine.
#ifdef FLASH_GRID_UG
  return
#else

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.false.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do


  !#ifdef SPECIAL_REFINEMENT
  !! Call for the specialized refinement
  specsSize=7
  !! Coordinate information --------------------------------------
  !! define a range of coordinates of the rectangle in x-direction

  !! |x| <= 0.05   for eta=1.e-2,   5.e-3,  2.5e-3
  !! |x| <= 0.025  for eta=1.25e-3, 6.25e-4,
  !! |x| <= 0.0125 for eta=3.125e-4, 1.5625e-4

  specs(1) =  0.25
  specs(2) =  0.75

  !! define a range of coordinates of the rectangle in y-direction
  specs(3) =  0.25
  specs(4) =  0.75 !! for alf07a
!!$  specs(4) = 1.5 !! for alf11b and alf12

  !! define a range of coordinates of the rectangle in z-direction
  specs(5) =  0.25
  specs(6) =  0.75
  !! End of coordinate information -------------------------------

  !! Decide wheather or not we refine only blocks completely
  !! contained within the rectangle (specs(7) .NE. 0.0)
  !! Otherwise, refine blocks with any overlap (specs(7) .EQ. 0.0)
  specs(7) = 1.0

  !! Bring all qualifying blocks to this level of refinement
  lref = lrefine_max

  call Grid_markRefineSpecialized (RECTANGLE,specsSize,specs,lref)
  !#endif



#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

#endif

  return
end subroutine Grid_markRefineDerefine

