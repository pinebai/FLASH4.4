!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgMapBcType
!!
!! NAME
!!  gr_mgMapBcType
!!
!! SYNOPSIS
!!
!!  call gr_mgMapBcType(integer(OUT)  :: bcTypeToApply,
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
!!  This routine is adapted from the HG solver. It allows to
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
!!  gr_bcMapBcType
!!  gr_bcApplyToOneFace
!!***

#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH3OR4
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

subroutine gr_mgMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
     axis,face,idest)

#include "constants.h"
#include "Multigrid.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY: gr_meshMe
  
#ifdef FLASH_GRID_UG
  use physicaldata, ONLY: unk,facevarx,facevary,facevarz
#endif
#ifdef FLASH_GRID_PARAMESH3OR4
  use physicaldata, ONLY: unk1,facevarx1, facevary1, facevarz1, gcell_on_cc, gcell_on_fc
  use workspace, ONLY : work1
#endif

#ifdef FLASH_GRID_PARAMESH2
  use workspace, ONLY : work
  use physicaldata, ONLY : unk
#endif

  use gr_mgData, ONLY: gr_mgBndTypes

  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
  integer,intent(IN),OPTIONAL:: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct

  bcTypeToApply = bcTypeFromGrid

  if (gridDataStruct == WORK) then
     select case (bcTypeFromGrid)
                                             !These are added for the Incomp Navier Stokes Solver
     case(OUTFLOW,DIODE,REFLECTING,DIRICHLET,SLIP_INS,NOSLIP_INS,MOVLID_INS,INFLOW_INS,OUTFLOW_INS)
           
           select case (gr_mgBndTypes(2*axis-1))

           case(MG_BND_DIRICHLET,MG_BND_GIVENVAL)
              if (present(idest)) then

                 bcTypeToApply = DIRICHLET

#ifdef DEBUG_GRID
                 if (idest .NE. 1) &
                      print*,'gr_mgMapBcType:',bcTypeToApply,'<-',bcTypeFromGrid,varIndex,gridDataStruct,idest   ! within DEBUG
#endif DEBUG_GRID
              end if

           case(MG_BND_NEUMANN)
              if (present(idest)) then

                 bcTypeToApply = NEUMANN_INS

#ifdef DEBUG_GRID
                 if (idest .NE. 1) &
                      print*,'gr_mgMapBcType:',bcTypeToApply,'<-',bcTypeFromGrid,varIndex,gridDataStruct,idest   ! within DEBUG 
#endif DEBUG_GRID
              end if
           end select

     case default
           if (gr_meshMe .EQ. MASTER_PE) print*, &
                "gr_mgMapBcType: An unexpected bcTypeFromGrid was requested for the WORK array:",bcTypeFromGrid,idest
           call Driver_abortFlash("gr_mgMapBcType: An unexpected bcTypeFromGrid was requested for the WORK array!")
     end select
  end if

  return
end subroutine gr_mgMapBcType








