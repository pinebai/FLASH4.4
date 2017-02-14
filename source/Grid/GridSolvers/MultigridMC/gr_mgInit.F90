!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgInit
!!
!! NAME
!!
!!  gr_mgInit
!!
!! SYNOPSIS
!!
!!  call gr_mgInit()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_init()

!  Description: Multigrid initialization routine.


  subroutine gr_mgInit()

!===============================================================================

#include "Flash.h"

  use gr_mgData

  use tree, only : nodetype,newchild,lrefine,maxblocks_tr

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks 

  use Grid_data, ONLY : gr_geometry

  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "Multigrid.h"
#include "constants.h"
include "mpif.h"

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


!===============================================================================


! INITILIAZE PARAMESH multigrid specific data

! Determine index ranges for interior zones.
ili = 1 + nguard
iui = nxb + nguard
jli = 1 + k2d*nguard
jui = nyb + k2d*nguard
kli = 1 + k3d*nguard
kui = nzb + k3d*nguard

! Determine index ranges for exterior zones.
ile = 1
iue = nxb + 2*nguard
jle = 1
jue = nyb + 2*nguard*k2d
kle = 1
kue = nzb + 2*nguard*k3d

! Determine mesh geometry (only on the first call).
if (FirstCall) then

  igeom = gr_geometry
  call RuntimeParameters_get('quadrant',    quadrant)

  select case (ndim)
  case (1)
    if (igeom == CARTESIAN) then
      mg_geometry = MG_GEOM_1DCARTESIAN
!    elseif (igeomx == geom_cylrad) then
!      mg_geometry = MG_GEOM_1DCYLINDRICAL
!    elseif (igeomx == geom_sphrad) then
!      mg_geometry = MG_GEOM_1DSPHERICAL
    else
      mg_geometry = MG_GEOM_INVALID
    endif
  case (2)
    if ((igeom == CARTESIAN)) then
      mg_geometry = MG_GEOM_2DCARTESIAN
!    elseif ((igeomx == geom_cylrad) .and. &
!            (igeomy == geom_planar)) then
!      mg_geometry = MG_GEOM_2DCYLAXISYM
!    elseif ((igeomx == geom_cylrad) .and. &
!            (igeomy == geom_cylang)) then
!      mg_geometry = MG_GEOM_2DCYLPOLAR
!    elseif ((igeomx == geom_sphrad) .and. &
!            (igeomy == geom_sphtheta)) then
!      mg_geometry = MG_GEOM_2DSPHAXISYM
    else
      mg_geometry = MG_GEOM_INVALID
    endif
  case (3)
    if ((igeom == CARTESIAN)) then
      mg_geometry = MG_GEOM_3DCARTESIAN
!    elseif ((igeomx == geom_cylrad) .and. &
!            (igeomy == geom_cylang) .and. &
!            (igeomz == geom_planar)) then
!      mg_geometry = MG_GEOM_3DCYLINDRICAL
!    elseif ((igeomx == geom_sphrad) .and. &
!            (igeomy == geom_sphtheta) .and. &
!            (igeomz == geom_sphphi)) then
!      mg_geometry = MG_GEOM_3DSPHERICAL
    else
      mg_geometry = MG_GEOM_INVALID
    endif
  end select
  if (mg_geometry == MG_GEOM_INVALID) then

    call Driver_abortFlash("gr_mgInit:  invalid geometry type!")

  endif
  FirstCall = .false.
endif

!===============================================================================

return
end
