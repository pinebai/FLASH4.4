!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/gr_mgInitSlv
!!
!! NAME
!!
!!  gr_mgInitSlv
!!
!! SYNOPSIS
!!
!!  call gr_mgInitSlv(integer(in) :: bndtypes)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   bndtypes : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_initSlv()

!  Description: Multigrid initialization routine.


subroutine gr_mgInitSlv(bndTypes)

!===============================================================================

#include "Flash.h"

use gr_mgData

use amr_mg_common, ONLY : amr_mg_min_required_level

use tree, only : nodetype,newchild,lrefine,maxblocks_tr,grid_changed


use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                              Grid_getListOfBlocks, &
                              Grid_getMaxCommonRefinement

use Grid_data, ONLY : gr_geometry, gr_meshMe, gr_meshComm

use gr_mgPfftdata, ONLY : gr_mgPfftMaxDirectSolveLevel,gr_mgbcTypes

use Driver_interface, ONLY : Driver_abortFlash
use RuntimeParameters_interface, ONLY : RuntimeParameters_get
use Timers_interface, ONLY: Timers_start, Timers_stop

use physicaldata, only : diagonals

implicit none
#include "Multigrid.h"
#include "constants.h"
#include "Flash_mpi.h"

integer, intent(in) :: bndTypes(6)

integer            :: lb, i, lmin, lmax, ierr, lnblocks2
integer, parameter :: MAXDIM2 = 3
integer            :: nbr_blks(2*MAXDIM2)

logical, save :: FirstCall = .true.
integer       :: igeom
logical       :: bnd_is_valid

real, pointer, dimension(:,:,:,:) :: unkt
integer, pointer, dimension(:)  :: nodetype2
logical, pointer, dimension(:)  :: newchild2

integer, parameter :: nxb = NXB
integer, parameter :: nyb = NYB
integer, parameter :: nzb = NZB
integer, parameter :: nguard = NGUARD
integer, parameter :: ndim = NDIM
integer, parameter :: k2d = K2D
integer, parameter :: k3d = K3D

integer             :: eachboundary

character*70        :: internalFile
integer             :: comm, globalTopLevel
logical             :: requestMap
logical             :: suppressPfft
integer             :: gridChanged
integer,save        :: prevSolveLevel = -1

!===============================================================================

comm = gr_meshComm

!The PARAMESH variable "grid_changed" takes the value 1 when the
!grid has changed.  Otherwise, it has a value 0 when the grid is the same.
!We store its value locally because if it has value 1, it will be reset to
!value 0 in amr_mg_morton_process (called by amr_mg_init).
!gr_hgPfftInitGrid needs to know the original value of grid_changed.
gridChanged = grid_changed
!==========================================================================

! Determine the maximum common refinement level in the grid.
call Grid_getMaxCommonRefinement(comm, globalTopLevel)

!Choose the direct solve level: gr_hgPfftMaxDirectSolveLevel is used to
!force a solve on a coarser level than globalTopLevel.  We will ignore
!this value if it is greater than globalTopLevel.

solveLevel = min(globalTopLevel, gr_mgPfftMaxDirectSolveLevel)

if (any(bndTypes(1:2*NDIM).GT.MG_BND_GIVENVAL)) then
   solveLevel = 1
   suppressPfft = .TRUE.
else
   suppressPfft = .FALSE.
end if

if (solveLevel < prevSolveLevel) then
   print*,'Lowering solveLevel from',prevSolveLevel,' to',solveLevel,' grid_changed was',grid_changed
   grid_changed = 1
   prevSolveLevel = solveLevel
end if
if (prevSolveLevel == -1) prevSolveLevel = solveLevel

!This routine initializes the multigrid and sets up data
!arrays for data exchange operations for the different levels
!in the PARAMESH tree hierarchy.
call Timers_start("amr_morton_process")
amr_mg_min_required_level = solveLevel !Set value for amr_mg_init.
diagonals = .false.                    !Set it S.t. doesn't fill diagonals.
call amr_mg_init()
diagonals = .true.
call Timers_stop("amr_morton_process")

! Initialize source data.
! Initialize PFFT extensions.
!DEV: CD.  I think we are passing a factor of 1.0 because we are multiplying
!by the Poisson factor in gr_hgLevelMultiplyScalar, so we don't need to do it twice.
if (.NOT. suppressPfft) call gr_mgPfftInitGrid(solvelevel, gridChanged, 1.)

call Grid_getLocalNumBlks(lnblocks2)

do lb = 1, lnblocks2
   nodetype_save(lb) = nodetype(lb)
   newchild_save(lb) = newchild(lb)
end do
! Determine the smallest and largest levels of refinement in the mesh.
lmax = -1
lmin = 10000000
do i = 1, lnblocks2
  if (nodetype(i) == 1) then
    lmin = min( lrefine(i), lmin )
    lmax = max( lrefine(i), lmax )
  endif
enddo
call mpi_allreduce ( lmin, mesh_lrefmin, 1, MPI_INTEGER, MPI_MIN, &
                     MPI_COMM_WORLD, ierr )
call mpi_allreduce ( lmax, mesh_lrefmax, 1, MPI_INTEGER, MPI_MAX, &
                     MPI_COMM_WORLD, ierr )
! Assign Boundary condition types to gr_mgBndTypes, same as hg solver:
do eachBoundary = 1, 2*NDIM
   if ((bndTypes(eachBoundary) == gr_mgbcTypes(eachBoundary)) .OR. suppressPfft .OR. &
      (bndTypes(eachBoundary) == MG_BND_GIVENVAL .AND.                             &
       gr_mgbcTypes(eachBoundary)==MG_BND_DIRICHLET) ) then

      gr_mgBndTypes(eachBoundary) = bndTypes(eachBoundary)

   else

      if (gr_meshMe .eq. 0) then
         write(*,*) 'gr_mgSolve Error: Boundary Conditions for Poisson Solver is inconsistent.'
         write(*,*) 'gr_mgSolve Error: direction=',eachBoundary
         write(*,*) 'gr_mgSolve Error: gr_mgbcTypes(direction) =',gr_mgbcTypes(eachBoundary)
         write(*,*) 'gr_mgSolve Error: bndTypes(direction)     =',bndTypes(eachBoundary)
      endif
      call Driver_abortFlash('gr_mgSolve Error: BC type argument inconsistent')
   end if
end do

! Assign value to mg_bnd_cond: case 0, substract mean from source,
!                              case 1, subtract given value of solution as
!                                      in Dirichlet BCs.

mg_bnd_cond = 0 !Assume all BCs are Periodic or Neumann
do  eachBoundary = 1, 2*NDIM
  if ((bndTypes(eachBoundary)==MG_BND_GIVENVAL) .or. &
      (bndTypes(eachBoundary)==MG_BND_DIRICHLET)   ) &
      mg_bnd_cond = 1 ! Case Some BC is DIRICHLET-GIVENVAL
enddo
!===============================================================================

return
end
