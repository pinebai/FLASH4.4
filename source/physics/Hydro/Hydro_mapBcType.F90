!!****f* source/physics/Hydro/Hydro_mapBcType
!!
!! NAME
!!  Hydro_mapBcType
!!
!! SYNOPSIS
!!
!!  call Hydro_mapBcType(integer(OUT)  :: bcTypeToApply,
!!                      integer(IN)   :: bcTypeFromGrid,
!!                      integer(IN)   :: varIndex,
!!                      integer(IN)   :: gridDataStruct,
!!                      integer(IN)   :: axis,
!!                      integer(IN)   :: face,
!!             OPTIONAL,integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine returns a boundary condition type (BC type) for
!!  a specific variable in a specific grid datastructure, based on
!!  the BC type that is stored in the structures of the Grid that
!!  keep track of block neighbors.
!!
!!  Normally, this routine just returns the bcTypeFromGrid as bcTypeToApply,
!!  or 0 (which basically has the same effect).
!!  However, this routines gives the Hydro unit the opportunity to
!!  apply different types of boundary conditions from what the
!!  Grid would normally use.
!!
!! ARGUMENTS
!!
!!  bcTypeToApply :  The boundary condition type returned to the caller.
!!                   The caller should use this to determine the actual
!!                   handling of the boundary condition, for "this"
!!                   variable (given by varIndex), data structure (given
!!                   by gridDataStruct), and PARAMESH array slot (given
!!                   by idest).
!!  bcTypeFromGrid : The BC type that is stored in the structures of the
!!                   Grid that keep track of block neighbors - in particular,
!!                   the "neigh" (Paramesh2) and "surr_blks" (Paramesh3 ff)
!!                   arrays. Typically this is the BC type set with the
!!                   {x,y,z}{l,u}_boundary_type runtime parameters.
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
!!  face           : The face on which to apply BC. Can take values LOW and HIGH,
!!                   for uper and lower (left and right) side of the block,
!!                   respectively. Note that LOW and HIGH are defined in constants.h
!!                   to have the values 1 and 2, respectively.
!! 
!!  axis           : the dimension along which to apply BC
!!
!!  idest:           Can usually be ignored; see Grid_bcApplyToRegion and gr_bcMapBcType
!!                   if you really want to know more.
!!
!!  HISTORY
!!    2013       Created                                                - Klaus Weide
!!
!! SEE ALSO
!!  Grid_bcApplyToRegion
!!  gr_bcMapBcType
!!  gr_bcApplyToOneFace
!!***

subroutine Hydro_mapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
     axis,face,idest)

  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
  integer,intent(IN),OPTIONAL:: idest

  bcTypeToApply = 0

  return
end subroutine Hydro_mapBcType








