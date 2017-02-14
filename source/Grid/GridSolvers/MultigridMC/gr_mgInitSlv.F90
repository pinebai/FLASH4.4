!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgInitSlv
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

use tree, only : nodetype,newchild,lrefine,maxblocks_tr

use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                              Grid_getListOfBlocks

use Grid_data, ONLY : gr_geometry

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

integer :: eachboundary

!===============================================================================

!This routine initializes the multigrid and sets up data
!arrays for data exchange operations for the different levels
!in the PARAMESH tree hierarchy.
call Timers_start("amr_morton_process")
diagonals = .false.                    !Set it S.t. doesn't fill diagonals.
call amr_mg_init()
diagonals = .true.
call Timers_stop("amr_morton_process")

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
   gr_mgBndTypes(eachBoundary) = bndTypes(eachBoundary)
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
