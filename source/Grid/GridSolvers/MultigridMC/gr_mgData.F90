!*******************************************************************************

!  Module:      gr_mgData()

!  Description: Global type and variable declarations for the multigrid solver.

module gr_mgData



implicit none

#include "constants.h"  ! for MDIM

!===============================================================================


!logical, save ::  writeflag_mg

! Coarse Solution level
integer,save :: Solvelevel

! Largest and smallest levels of refinement in the mesh.
integer,save :: mesh_lrefmax, mesh_lrefmin

! A place to save nodetypes as paramesh has them before we change
! them to perform various one level operations.
integer,save, allocatable, dimension(:) :: nodetype_save !(maxblocks_tr)

! A place to save child data before the prolongation messes with it
logical,save, allocatable, dimension(:) :: newchild_save !(maxblocks_tr)

! Ranges of interior indices for blocks.

integer,save :: ili, iui, jli, jui, kli, kui

! Ranges of exterior indices for blocks.

integer,save :: ile, iue, jle, jue, kle, kue

! Boundary conditions, flag to substract mean on the source.
integer,save, dimension(2*MDIM) :: gr_mgBndTypes
integer,save :: mg_bnd_cond

! Grid geometry.

integer,save :: mg_geometry

! Just doing a quadrant?

logical,save :: quadrant

!     Second order interpolation for the work array:
integer, save :: interp_work = 2
integer, save, allocatable, dimension(:) :: interp_mask_work_mg
integer, save, allocatable, dimension(:) :: interp_mask_work_save

integer, save :: gr_mgDiffOpDiscretize

! Smoother types:
integer, parameter :: rbgs_smoother = 1
integer, parameter :: zebra_smoother= 2

!===============================================================================

end module gr_mgData
