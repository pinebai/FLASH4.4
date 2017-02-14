!!****if* source/Grid/localAPI/gr_bicgMapBcType
!!
!! NAME
!!  gr_bicgMapBcType
!!
!! SYNOPSIS
!!
!!  call gr_bicgMapBcType(integer(OUT)  :: bcTypeToApply,
!!                        integer(IN)   :: bcTypeFromGrid,
!!                        integer(IN)   :: varIndex,
!!                        integer(IN)   :: gridDataStruct,
!!                        integer(IN)   :: axis,
!!                        integer(IN)   :: face,
!!                        integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine returns a boundary condition type (BC type) for
!!  a specific variable in a specific grid datastructure, based on
!!  the BC type that is stored in the structures of the Grid that
!!  keep track of block neighbors.
!!
!!  Normally, this routine just returns the bcTypeFromGrid as bcTypeToApply.
!!  This routine is adapted from the HG solver. It gives the solver the opportunity to
!!  apply different types of boundary conditions from what the
!!  Grid would normally use.
!!
!! ARGUMENTS
!!
!!  bcTypeToApply :  Returns the BC type that the caller should
!!                   (perhaps after further mapping) apply.
!!                   The caller should use this to determine the actual
!!                   handling of the boundary condition, for "this"
!!                   variable (given by varIndex), data structure (given
!!                   by gridDataStruct), and PARAMESH array slot (given
!!                   by idest).
!!
!!  bcTypeFromGrid:  the BC type as known to the parts of the Grid unit
!!                   that keep track of boundary conditions on a per-direction
!!                   and/or per-block (but not, per-variable) basis
!!
!!  varIndex :       Index that selects a variable within the data
!!                   structure given by the gridDataStruct argument (below).
!!                   For gridDataStruct==CENTER, this is an UNK variable
!!                   such as DENS_VAR. For WORK, this is meaningless
!!                   and should be set to 1.
!!
!!  gridDataStruct : Integer value specifying the data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   WORK   work array specific to paramesh
!!
!!  idest:           Index that selects a slot within the UNK1 / WORK1 / ..
!!                   Paramesh 4 array for which guard cells are being filled.
!!                   Meaningless for Paramesh2.
!!                   When guard cells are being filled as an effect of
!!                   calling the Paramesh 4 implementation of
!!                   amr_1blk_guardcell, then the value of idest can be
!!                   used to distinguish between the following -
!!                   idest == 1 : currently filling guard cells in a block that
!!                                are actually to be returned to the user.
!!                   idest == 2 : currently filling guards of a parent block of
!!                                another block, for the purpose of serving as
!!                                input to coarse-to-fine interpolation.
!!
!!  axis :           the dimension for applying BC
!!
!!  face :           the face of for applying BC
!!
!! SEE ALSO
!!  gr_bcMapBcType
!!  gr_bcApplyToOneFace
!!***

subroutine gr_bicgMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
     axis,face,idest)
  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
  integer,intent(IN),OPTIONAL:: idest

  bcTypeToApply = 0
  return
end subroutine gr_bicgMapBcType
