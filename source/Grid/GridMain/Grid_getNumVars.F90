!!****if* source/Grid/GridMain/Grid_getNumVars
!!
!! NAME
!!  Grid_getNumVars
!! 
!! SYNOPSIS
!!  Grid_getNumVars(integer, intent(in)  :: gridStruct,
!!                  integer, intent(out) :: nVar)
!!
!! DESCRIPTION
!!  Returns the number of variables in each mesh data
!!  structure.  Allows us to write more flexible code that
!!  is not dependent on NUNK_VARS and NFACE_VARS e.t.c.
!!
!! ARGUMENTS
!!  gridStruct -- Integer value representing mesh data structure
!!  nVar -- Number of mesh variables in gridStruct
!!
!!***
#include "constants.h"
#include "Flash.h"

subroutine Grid_getNumVars(gridStruct, nVar)  
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none
  integer, intent(in) :: gridStruct
  integer, intent(out) :: nVar

  select case (gridStruct)
  case(CENTER)
     nVar = NUNK_VARS
  case(FACEX,FACEY,FACEZ)
     nVar = NFACE_VARS
  case(SCRATCH)
     nVar = NSCRATCH_GRID_VARS
  case(SCRATCH_CTR)
     nVar = NSCRATCH_CENTER_VARS
  case(SCRATCH_FACEX)
     nVar = NSCRATCH_FACEX_VARS
  case(SCRATCH_FACEY)
     nVar = NSCRATCH_FACEY_VARS
  case(SCRATCH_FACEZ)
     nVar = NSCRATCH_FACEZ_VARS
  case DEFAULT
     call Driver_abortFlash("[Grid_getNumVars]: Invalid data structure")
  end select
end subroutine Grid_getNumVars
