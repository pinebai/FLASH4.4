!!****f* source/Grid/Grid_getNumVars
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
subroutine Grid_getNumVars(gridStruct, nVar)  
  implicit none
  integer, intent(in) :: gridStruct
  integer, intent(out) :: nVar
  nVar = 0
end subroutine Grid_getNumVars
