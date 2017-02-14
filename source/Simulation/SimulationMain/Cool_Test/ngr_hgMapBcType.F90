!!****if* source/Simulation/SimulationMain/Cool_Test/ngr_hgMapBcType
!!
!! NAME
!!  gr_hgMapBcType
!!
!! SYNOPSIS
!!
!!  call gr_hgMapBcType(integer(OUT)  :: bcTypeToApply,
!!                      integer(IN)   :: bcTypeFromGrid,
!!                      integer(IN)   :: varIndex,
!!                      integer(IN)   :: gridDataStruct,
!!                      integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine returns a boundary condition type (BC type) for
!!  a specific variable in a specific grid datastructure, based on
!!  the BC type that is stored in the structures of the Grid that
!!  keep track of block neighbors.
!!
!!  Normally, this routine just returns the bcTypeFromGrid as bcTypeToApply.
!!  However, this routines gives the HG solver the opportunity to
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
!!  idest:           Index that selects a slot within the UNK1 / WORK1 / ..
!!                   Paramesh3 array for which guard cells are being filled.
!!                   Meaningless for Paramesh2.
!!                   When guard cells are being filled as an effect of
!!                   calling the Paramesh3 / Paramesh4 implementation of
!!                   amr_1blk_guardcell, then the value of idest can be
!!                   used to distinguish between the following -
!!                   idest == 1 : currently filling guard cells in a block that
!!                                are actually to be returned to the user.
!!                   idest == 2 : currently filling guards of a parent block of
!!                                another block, for the purpose of serving as
!!                                input to coarse-to-fine interpolation.
!!
!! SEE ALSO
!!  gr_bcMapBcType
!!  gr_bcApplyToOneFace
!!***

#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH4DEV
#define FLASH_GRID_PARAMESH3
#endif
#ifdef FLASH_GRID_PARAMESH3
!!REORDER(5): unk1, facevar1[xyz]
#endif
#ifdef FLASH_GRID_UG
!!REORDER(5): unk, facevar[xyz]
#endif
#ifdef FLASH_GRID_PARAMESH2
!!REORDER(5) : unk
#endif

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_hgMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct,idest)

#include "constants.h"
#include "Multigrid.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY: gr_meshMe
  
#ifdef FLASH_GRID_UG
  use physicaldata, ONLY: unk,facevarx,facevary,facevarz
#endif
#ifdef FLASH_GRID_PARAMESH3
  use physicaldata, ONLY: unk1,facevarx1, facevary1, facevarz1, gcell_on_cc, gcell_on_fc
  use workspace, ONLY : work1
#endif

#ifdef FLASH_GRID_PARAMESH2
  use workspace, ONLY : work
  use physicaldata, ONLY : unk
#endif

  use gr_hgData, ONLY: gr_hgNowActive, gr_hgCurrentGcReq_extrap, gr_hgBndConditions

  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct
  integer,intent(IN),OPTIONAL:: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct


#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     if (gr_meshMe .EQ. MASTER_PE) print *, "gr_hgMapBcType: gridDataStruct set to improper value"
     if (gr_meshMe .EQ. MASTER_PE) print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abortFlash("gr_hgMapBcType gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varIndex/=1)) &
       call Driver_abortFlash("gr_hgMapBcType: varCount be 1 for work array")
  if((gridDataStruct==CENTER).and.(varIndex > NUNK_VARS)) &
       call Driver_abortFlash("gr_hgMapBcType: varIndex be <= NUNK_VARS for unk array")
       

#endif
  
  bcTypeToApply = bcTypeFromGrid

  if (gr_hgNowActive) then
     if (gridDataStruct == WORK) then
        select case (bcTypeFromGrid)
        case(OUTFLOW,DIODE,REFLECTING)
!!$           bcTypeToApply = GRIDBC_GIMME_WORK
!!$           return

           select case (gr_hgBndConditions)
           case(MG_BND_DIRICHLET,MG_BND_GIVENVAL)
              if (present(idest)) then
                 select case (idest)
                 case(1)
                    if (gr_hgCurrentGcReq_extrap) then
                       bcTypeToApply = GRIDBC_MG_EXTRAPOLATE
                    else
                       bcTypeToApply = DIRICHLET
                    end if
                 case(2)
                    bcTypeToApply = GRIDBC_MG_EXTRAPOLATE
                 end select
#ifdef DEBUG_GRID
           if (idest .NE. 1) &
           print*,'gr_hgMapBcType:',bcTypeToApply,'<-',bcTypeFromGrid,varIndex,gridDataStruct,idest   ! within DEBUG 
#endif DEBUG_GRID
              end if
           case(MG_BND_NEUMANN)
              bcTypeToApply = DIRICHLET
           case default
              if (gr_meshMe .EQ. MASTER_PE) print*,  &
                   "gr_hgMapBcType: An unexpected gr_hgBndConditions was requested for the WORK array:", &
                 gr_hgBndConditions, bcTypeFromGrid, idest
              call Driver_abortFlash("gr_hgMapBcType: Unexpected boundary condition requested for WORK array!")
           end select

        case default
           bcTypeToApply = OUTFLOW
        !   if (gr_meshMe .EQ. MASTER_PE) print*, &    
!   "gr_hgMapBcType: An unexpected bcTypeFromGrid was requested for the WORK array:",bcTypeFromGrid,idest
        !   call Driver_abortFlash("gr_hgMapBcType: An unexpected bcTypeFromGrid was requested for the WORK array!")
        end select
     end if
  end if

  return
end subroutine gr_hgMapBcType








