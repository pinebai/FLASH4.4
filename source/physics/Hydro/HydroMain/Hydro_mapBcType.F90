!!****if* source/physics/Hydro/HydroMain/Hydro_mapBcType
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
!!                   to have the values 1 and 2, respectively, and this implementation
!!                   of Hydro_mapBcType assumes these definitions.
!! 
!!  axis           : the dimension along which to apply BC
!!
!!  idest:           Can usually be ignored; see Grid_bcApplyToRegion and gr_bcMapBcType
!!                   if you really want to know more.
!!
!!  HISTORY
!!
!!    2013       Created                                                - Klaus Weide
!!    2013       Correct treatment of face-centered B fields outflow BC - Dongwook Lee
!!
!! SEE ALSO
!!  Grid_bcApplyToRegion
!!  gr_bcMapBcType
!!  gr_bcApplyToOneFace
!!***

#include "Flash.h"

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif
#define DEBUG_GRID

subroutine Hydro_mapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
     axis,face,idest)

#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_data, ONLY: hy_meshMe
  

  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
  integer,intent(IN),OPTIONAL:: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct
  logical :: isMagComponent


#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     if (hy_meshMe .EQ. MASTER_PE) print *, "Hydro_mapBcType: gridDataStruct set to improper value"
     if (hy_meshMe .EQ. MASTER_PE) print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abortFlash("Hydro_mapBcType gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varIndex/=1)) &
       call Driver_abortFlash("Hydro_mapBcType: varCount be 1 for work array")
  if((gridDataStruct==CENTER).and.(varIndex > NUNK_VARS)) &
       call Driver_abortFlash("Hydro_mapBcType: varIndex be <= NUNK_VARS for unk array")
       

#endif
  
  bcTypeToApply = 0
  return  
  isMagComponent = .FALSE.
#ifdef MAG_FACE_VAR
  if (varIndex==MAG_FACE_VAR) isMagComponent = .TRUE.
#endif
#ifdef MAGI_FACE_VAR
  if (varIndex==MAGI_FACE_VAR) isMagComponent = .TRUE.
#endif


  if (isMagComponent) then
                 !! Note: A special consideration is applied for face variables using outflow BC.
                 !! The conventional outflow treatment breaks divB=0 condition that is subject to
                 !! satisfy for the constrained-transport (CT) facevariables. 
                 !! There are two ways to do this depending on to which variables Neumann BC applies:
                 !! (1) apply Neumann condition to the two transverse (to the boundary or asix) fields,
                 !!     then calculate the normal field using the divB=0 relationship as a closure,
                 !! (2) apply Neumann condition to the normal field, then seek for the transverse fields
                 !!     that satisfy divB=0. A simple way is to set them as constants.
                 !!     For example, if we set Neumann for Bx along x-direction, then assigning 0 to
                 !!     By and Bz will satisfy divB=0 automatically.
                 !! We adopt the approach (2) for FLASH implementation as there are two advantages:
                 !! (i)  BC handlings in Paramesh can only work for each gridDataStruct at a time.
                 !!      This becomes a limitation for the approach (1) because we need
                 !!      FACEX, FACEY and FACEZ all at the same time to close Bx using Neumann By & Bz.
                 !!      On the other hand, the approach (2) only requires each FACE* gridDataStruct
                 !!      without needing the other two for closure.
                 !! (ii) The normal fields have more meaningful usages in calculations. For instance,
                 !!      They are directly used for n+1/2 Riemann states, whereas the corresponding
                 !!      states for transverse fields are always reconstructed. In this sense, it is
                 !!      better to apply Neumann for the normal fields, while any constant values can
                 !!      be assigned to the transverse fields. Such constant values are then filled,
                 !!      satisfying divB=0, without affecting Riemann state reconstructions because
                 !!      it is the (Neumann-treated) cell-centered B-fields that are used for 
                 !!      the transverse fields reconstructions.
     if (gridDataStruct == FACEX .OR. gridDataStruct == FACEY .OR. gridDataStruct == FACEZ) then
        if ((gridDataStruct == FACEX .and. axis == IAXIS) .or. &
            (gridDataStruct == FACEY .and. axis == JAXIS) .or. &
            (gridDataStruct == FACEZ .and. axis == KAXIS)) then
           select case (bcTypeFromGrid)
           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
              bcTypeToApply = REFLECTING
           end select
        end if
        if ((gridDataStruct == FACEX .and. axis .NE. IAXIS) .or. &
            (gridDataStruct == FACEY .and. axis .NE. JAXIS) .or. &
            (gridDataStruct == FACEZ .and. axis .NE. KAXIS)) then
           select case (bcTypeFromGrid)
           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
              bcTypeToApply = GRIDBC_ZERO
           end select
        end if
     end if
  end if

  return
end subroutine Hydro_mapBcType








