!!****f* source/Grid/Grid_applyBCEdge
!!
!! NAME
!!  Grid_applyBCEdge
!!
!! SYNOPSIS
!!
!!  Grid_applyBCEdge(integer(IN)              :: bcType, 
!!                   integer(IN)              :: bcDir,
!!                   integer(IN)              :: guard,
!!                   integer(IN)              :: var,
!!                   real(INOUT),dimension(:) :: dataRow(2*guard),
!!                   integer(IN)              :: face,
!!                   integer(IN)              :: gridDataStruct,
!!                   integer(IN),OPTIONAL     :: blockHandle,
!!                   real(in),OPTIONAL    :: secondCoord,
!!                   real(in),OPTIONAL    :: thirdCoord)
!!  
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Applies the boundary conditions to a given vector.
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) of a given vector. 
!!     If (face=LOW)dataRow(1:guard) = boundary values
!!     If (face=HIGH) dataRow(guard+1:2*guard) = boundary values
!!
!!  The reason why information about direction and variable is
!!  included in this interface is because velocities need to be
!!  treated specially for REFLECTING boundary conditions. if
!!  bcDir=IAXIS, then the variable VELX_VAR is treated differently,
!!  same with VELY_VAR if bcDir=JAXIS and VELZ_VAR if
!!  bcDir=KAXIS. This special treatment is currently only done if
!!  gridDataStruct=CENTER (face variables are assumed not to be
!!  velocities). All supported mesh packages extract the vector passed
!!  in through the argument "dataRow" from the appropriated blocks,
!!  and send it to this routine for boundary calculation. This routine
!!  currently calculates REFLECTING and OUTFLOW boundaries. The
!!  PERIODIC boundary is calculated by default when the blocks are
!!  exchanging data with each other.
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
!!  var   -    The variable on which boundary conditions are applied
!!             It is used with bcDir for reflective boundary conditions
!!             to correctly handle velocities. This argument is redundant 
!!             for all other variables.
!!  dataRow -  storage for the data being operated upon.
!!  face    -  can take values LOW and HIGH, defined in constants.h
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (unk or work for PM) (default)
!!                   FACEX  face centered variable on faces normal to IAXIS
!!                   FACEY  face centered variable on faces normal to JAXIS
!!                   FACEZ  face centered variable on faces normal to KAXIS
!!   blockHandle - the identity of the block under consideration
!!  secondCoord,thirdCoord - scalar coordinate values in the coordinate
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
!!                         These are either coordinates of cell centers or of
!!                         face centers, as appropriate depending on gridDataStruct.
!!
!!  NOTES 
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


subroutine Grid_applyBCEdge(bcType,bcDir,guard,var,dataRow,face,gridDataStruct, blockHandle,&
     secondCoord, thirdCoord)
  implicit none

  integer, intent(in):: bcType
  integer,intent(IN) :: var,guard,face,bcDir,gridDataStruct
  real,dimension(:),intent(INOUT)::dataRow
  integer,intent(IN),OPTIONAL :: blockHandle
  real,intent(IN),OPTIONAL :: secondCoord,thirdCoord
  
  return
end subroutine Grid_applyBCEdge
