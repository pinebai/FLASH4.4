!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bicgInitSlv
!!
!! NAME
!!
!!  gr_bicgInitSlv
!!
!! SYNOPSIS
!!
!!  call gr_bicgInitSlv(integer(in) :: bndtypes)
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

!  Routine:     bicg_initSlv()

!  Description: BiPCGStab initialization routine.


  subroutine gr_bicgInitSlv(bndTypes)

!===============================================================================

#include "Flash.h"

  use gr_bicgData

!!$  use tree, only : nodetype,newchild,lrefine,maxblocks_tr
!!$
  use Grid_interface,    ONLY : GRID_PDE_BND_DIRICHLET, &
                                GRID_PDE_BND_GIVENVAL
!!$
!!$  use Grid_data, ONLY : gr_geometry
!!$
!!$  use Driver_interface, ONLY : Driver_abortFlash
!!$  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "constants.h"

  integer, intent(in) :: bndTypes(6)

  integer            :: lb, i

!!$, lmin, lmax, ierr, lnblocks2
!!$  integer, parameter :: MAXDIM2 = 3
!!$  integer            :: nbr_blks(2*MAXDIM2)

!!$  logical, save :: FirstCall = .true.
!!$  integer       :: igeom
!!$  logical       :: bnd_is_valid

!!$  real, pointer, dimension(:,:,:,:) :: unkt
!!$  integer, pointer, dimension(:)  :: nodetype2
!!$  logical, pointer, dimension(:)  :: newchild2

  integer :: eachboundary


!===============================================================================


! Assign Boundary condition types to gr_bicgBndTypes:
  do eachBoundary = 1, 2*NDIM

     gr_bicgBndTypes(eachBoundary) = bndTypes(eachBoundary)

  end do

! Assign value to bicg_bnd_cond: case 0, substract mean from source, 
!                                case 1, subtract given value of solution as
!                                        in Dirichlet BCs.
  bicg_bnd_cond = 0 !Assume all BCs are Periodic or Neumann
  do  eachBoundary = 1, 2*NDIM
    if ((bndTypes(eachBoundary)==GRID_PDE_BND_DIRICHLET) .or. &
        (bndTypes(eachBoundary)==GRID_PDE_BND_GIVENVAL) ) &
        bicg_bnd_cond = 1 ! Case Some BC is DIRICHLET
  enddo

!===============================================================================

return
end subroutine gr_bicgInitSlv
