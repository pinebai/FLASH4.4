!!****if* source/Grid/GridBoundaryConditions/gr_bcMapBcType
!!
!! NAME
!!  gr_bcMapBcType
!!
!! SYNOPSIS
!!
!!  call gr_bcMapBcType(integer(OUT)  :: bcTypeToApply,
!!                      integer(IN)   :: bcTypeFromGrid,
!!                      integer(IN)   :: varIndex,
!!                      integer(IN)   :: gridDataStruct,
!!                      integer(IN)   :: axis,
!!                      integer(IN)   :: face,
!!                      integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine returns a boundary condition type (BC type) for
!!  a specific variable in a specific grid datastructure, based on
!!  the BC type that is stored in the structures of the Grid that
!!  keep track of block neighbors.
!!
!!  Normally, this routine just returns the bcTypeFromGrid as bcTypeToApply.
!!  However, this routines gives the opportunity for handling certain
!!  variables differently by returning something else.
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
!!  face           : The face on which to apply BC
!! 
!!  axis           : the dimension along which to apply BC
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
!! SEE ALSO
!!  Grid_bcApplyToRegion
!!***
#define DEBUG_GRID

#include "Flash.h"

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_bcMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
     axis,face,idest)

#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_interface, ONLY : Hydro_mapBcType
  use gr_hgInterface, ONLY : gr_hgMapBcType
  use gr_mgInterface, ONLY : gr_mgMapBcType
  use gr_bicgInterface, ONLY : gr_bicgMapBcType  

  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
  integer,intent(IN),OPTIONAL:: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  integer :: bcType, bcTypeMapped
  logical :: validGridDataStruct


#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     print *, "gr_bcMapBcType: gridDataStruct set to improper value"
     print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abortFlash("gr_bcMapBcType gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varIndex/=1)) &
       call Driver_abortFlash("gr_bcMapBcType: varCount be 1 for work array")
  if((gridDataStruct==CENTER).and.(varIndex > NUNK_VARS)) &
       call Driver_abortFlash("gr_bcMapBcType: varIndex be <= NUNK_VARS for unk array")
       
#endif

  

  bcType = bcTypeFromGrid

!  write(*,*) 'bcType before hg=',bcType

  call Hydro_mapBcType(bcTypeMapped,bcType,varIndex,gridDataStruct,axis,face,idest)
  if (bcTypeMapped .NE. 0) bcType = bcTypeMapped

  call gr_hgMapBcType(bcTypeMapped,bcType,varIndex,gridDataStruct,axis,face,idest)
  if (bcTypeMapped .NE. 0) bcType = bcTypeMapped

  call gr_mgMapBcType(bcTypeMapped,bcType,varIndex,gridDataStruct,axis,face,idest)
  if (bcTypeMapped .NE. 0) bcType = bcTypeMapped

!  write(*,*) 'bcType before bicg=',bcType

  call gr_bicgMapBcType(bcTypeMapped,bcType,varIndex,gridDataStruct,axis,face,idest)
  if (bcTypeMapped .NE. 0) bcType = bcTypeMapped


!!  Repeat for other subunits that should get a chance to further map BC types:
!!$  call gr_SUBUNIT1MapBcType(bcTypeMapped,bcType,varIndex,gridDataStruct,axis,face,idest)
!!$  if (bcTypeMapped .NE. 0) bcType = bcTypeMapped
!!$
!!$  call gr_SUBUNIT2MapBcType(bcTypeMapped,bcType,varIndex,gridDataStruct,axis,face,idest)
!!$  if (bcTypeMapped .NE. 0) bcType = bcTypeMapped
!!  etc....


  if(bcType==DIODE) then
     if (gridDataStruct.NE.CENTER) bcType = OUTFLOW           ! Treat as zero-gradient.

  else if(bcType==HYDROSTATIC_F2_NVDIODE .OR. bcType==HYDROSTATIC_F2_NVOUT .OR. bcType==HYDROSTATIC_F2_NVREFL) then
     if (gridDataStruct.NE.CENTER) bcType = OUTFLOW           ! Treat as zero-gradient.

  else if(bcType==HYDROSTATIC_NVDIODE .OR. bcType==HYDROSTATIC_NVOUT .OR. bcType==HYDROSTATIC_NVREFL) then
     if (gridDataStruct.NE.CENTER) bcType = OUTFLOW           ! Treat as zero-gradient.
  end if

  bcTypeToApply = bcType


  return
end subroutine gr_bcMapBcType








