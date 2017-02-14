!!****if* source/Grid/GridBoundaryConditions/OneRow/Grid_applyBCEdgeAllUnkVars
!!
!! NAME
!!  Grid_applyBCEdgeAllUnkVars
!!
!! SYNOPSIS
!!
!!  Grid_applyBCEdgeAllUnkVars(integer(in) :: bcType, 
!!                             integer(in) :: bcDir,
!!                             integer(in) :: guard,
!!                             real(INOUT) :: dataRow(2*guard,NUNK_VARS),
!!                             integer(in) :: face,
!!                             real(in)    :: cellCenterSweepCoord(*),
!!                             real(in)    :: secondCoord,
!!                             real(in)    :: thirdCoord,
!!                             integer(in),OPTIONAL :: blockHandle)
!!  
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Applies the boundary conditions to all unk variables at a given column
!!  of cells.
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) along a given axis (bcDir), by using and setting values
!!  for all variables (unless perhaps masked out) along a 1-dimensional
!!  column of cells perpendicular to the face.
!!     If (face=LOW)  dataRow(1:guard,:) = boundary values
!!     If (face=HIGH) dataRow(guard+1:2*guard,:) = boundary values
!!
!!  One reason why information about direction and variable is
!!  included in this interface is because Velocities need to be
!!  treated specially for REFLECTING boundary conditions. if
!!  bcDir=IAXIS, then the variable VELX_VAR is treated differently,
!!  same with VELY_VAR if bcDir=JAXIS and VELZ_VAR if
!!  bcDir=KAXIS. All supported mesh packages extract the vector passed
!!  in through the argument "dataRow" from the appropriated blocks,
!!  and send it to this routine for boundary calculation. The
!!  PERIODIC boundary is calculated by default when the blocks are
!!  exchanging data with each other.
!!  This routine currently passes handling of all other boundary condition
!!  types on to Grid_applyBCEdge, which is called for each of the variables
!!  in turn.
!!
!!  If the user wishes to apply different boundary conditions, they should
!!  make a copy of this routine in their given Simulation setup and implement 
!!  their specific boundary. 
!!
!!  
!! ARGUMENTS 
!!
!!
!!  bcType -   the type of boundary condition being applied to this face
!!  bcDir -    can take on values IAXIS,JAXIS or KAXIS. This is needed
!!             for handling the reflective boundary conditions. If bcDir=IAXIS,
!!             and boundary conditions are reflective, and X velocity is
!!             treated differently from all other variables. similarly if bcDir
!!             is JAXIS, then Y velocity is different.
!!  guard -    number of guardcells 
!!  dataRow -  storage for the data being operated upon.
!!  face    -  can take values LOW and HIGH, defined in constants.h
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  cellCenterSweepCoord - vector of (at least) 2*guard cell center coordinate
!!                         values in the sweep direction. This is not needed
!!                         for simple boundary condition types such as REFLECTIVE
!!                         or OUTFLOW, but is provided so that more complex
!!                         boundary conditions can make use of it.
!!  secondCoord,thirdCoord - scalar cell center coordinate values in the coordinate
!!                         directions perpendicular to the sweep direction.
!!                         This is not needed for simple boundary condition types
!!                         such as REFLECTIVE or OUTFLOW, but is provided so that
!!                         more complex boundary conditions can make use of it.
!!                         The meaning depends on the sweep direction bcDir as
!!                         follows
!!                          bcDir   |    secondCoord       thirdCoord
!!                          ------------------------------------------
!!                          IAXIS   |    Y(j) *            Z(k) **
!!                          JAXIS   |    X(i)              Z(k) **
!!                          KAXIS   |    X(i)              Y(j)
!!                         *)  if NDIM > 1
!!                         **) if NDIM > 2
!!                         These dummy arguments are ignored (and an implementation
!!                         of this interface should not attempt to access them) if
!!                         they do not make sense based on the dimensionality of
!!                         the problem.
!!  blockHandle - available for passing on to some handlers for boundary conditions
!!                that may need it, ignored in the default BC handling implementations
!!
!!
!!  NOTES 
!!
!!            This routine is common to all the mesh packages supported
!!            The mesh packages extract the small vector relevant to
!!            boundary conditions calculations from their Grid data 
!!            structures. If users wish to apply a different but simple
!!            directional boundary condition, they can replace this routine.
!!            More complex boundary conditions support will come in 
!!            future releases. In the meantime, users could use the
!!            Grid_applyBC.F90 file as a template. That file is currently
!!            used only for some Grid implementations (UG and Paramesh2, not
!!            for Paramesh 4), but gives a good outline for handling complex
!!            boundary conditions.
!!            
!!
!!***
#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_applyBCEdgeAllUnkVars(bcType,bcDir,guard,dataRow,face,&
     cellCenterSweepCoord, secondCoord,thirdCoord, blockHandle)
  use Grid_data, ONLY :gr_meshMe
  use Grid_interface, ONLY : Grid_applyBCEdge
  use Driver_interface, ONLY : Driver_abortFlash
#include "Flash.h"
#ifdef FLASH_GRID_PARAMESH3OR4
  use physicaldata, ONLY : gcell_on_cc
#define MASKED_IN(ivar) gcell_on_cc(ivar)
#else
#define MASKED_IN(ivar) .TRUE.
#endif
  implicit none
#include "constants.h"

  integer,intent(IN):: bcType,bcDir,guard,face
  real,dimension(2*guard,NUNK_VARS),intent(INOUT)::dataRow
  real,intent(IN):: cellCenterSweepCoord(*), secondCoord,thirdCoord
  integer,intent(IN),OPTIONAL:: blockHandle ! passed on in the default implementation - KW
  integer :: ivar


  if(bcType==PERIODIC) return

  do ivar = 1,NUNK_VARS
     if(MASKED_IN(ivar)) then
        call Grid_applyBCEdge(bcType,bcDir,guard,ivar,dataRow(:,ivar),face,CENTER, blockHandle, secondCoord, thirdCoord)
     endif
  end do

  return
end subroutine Grid_applyBCEdgeAllUnkVars
