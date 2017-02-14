!*******************************************************************************

!  Routine:     Grid_SolvePoisson()

!  Description: Driver routine for the multigrid Poisson solver.  This routine
!               interprets the supplied boundary conditions and either calls
!               the Poisson solver directly (in the case of periodic or
!               Dirichlet boundaries) or uses James' image-mass method to
!               handle isolated boundaries.  The primary purpose of this
!               routine is to provide an interface to the Poisson solver that
!               understands the boundary condition logic, and as such it is not
!               particularly important what Poisson solver algorithm is used,
!               as long as the routine has the same interface as multigrid().

!  Parameters:  isrc            Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               isoln           Index for solution array.  The solution is
!                               written directly into this variable.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: The following have been added to the interface only to get this
! version of the paramesh3.x poisson solver to work.  These have no effect.
!               bcTypes(6)     Boundary condition types array:
!                                 MG_BND_PERIODIC
!                                 MG_BND_DIRICHLET
!                                 MG_BND_NEUMANN
!
!                                 index 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z  6 = +z
!               bcValues(2,6)  Values for dirichlet and neumann boundary
!                               conditions.  If Robins, then bcValues(1,*) holds
!                               dirichlet and bcValues(2,*) neumann conditions; if
!                               neumann or dirichlet, then bcValues(1,*) holds
!                               neumann or dirichlet values and bcValues(2,*) goes unused
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries
!               poisfact        Constant Poisson factor.  Used to scale the
!                               source function.  For example, for gravity,
!                               poisfact = 4*pi*G.


subroutine Grid_SolvePoisson (isoln, isrc, bcTypes, bcValues, poisFact)

!===============================================================================
#include "Flash.h"

use gr_mgData

use Grid_interface,    ONLY : Grid_getListOfBlocks, &
                              Grid_getBlkPtr, &
                              Grid_releaseBlkPtr

use gr_mgInterface,  ONLY: gr_mgSolve
use gr_isoInterface, ONLY: gr_isoImageMass, gr_isoImageBdry
use Timers_interface, ONLY : Timers_start, Timers_stop
use Driver_interface, ONLY : Driver_abortFlash
use RuntimeParameters_interface, ONLY : RuntimeParameters_get


implicit none

integer,intent(in) :: isoln, isrc
integer,intent(in), dimension(6) :: bcTypes
real,intent(in), dimension(2,6) :: bcValues
real,intent(inout)    :: poisFact

integer       :: lb
integer, save :: imgw1, imgw2, imgw3, imgw4, imgw5, imgw6, imgw7

#ifdef IMGP_VAR
integer :: blockCount
integer :: blockList(MAXBLOCKS)
integer :: blockID
integer, dimension(6) :: bcTypes2
#endif

real, pointer, dimension(:,:,:,:) :: unkt

integer :: mgrid_smoother

external poisson_mg_solve, poisson_mg_residual, poisson_mg_residual_leafs, &
         poisson_mg_relax, poisson_mg_relax_ZEBRA

#ifdef IMGP_VAR
procedure(),POINTER :: relaxerPtr
#endif

!===============================================================================

! Get key numbers from the database for the temporary variables we need.

imgw1 = MGW1_VAR
imgw2 = MGW2_VAR
imgw3 = MGW3_VAR
imgw4 = MGW4_VAR
imgw5 = MGW5_VAR
imgw6 = MGW6_VAR
imgw7 = MGW7_VAR


call RuntimeParameters_get("mgrid_smoother",mgrid_smoother)


!===============================================================================

! Call the multigrid Poisson solver for different types of boundary conditions.

#ifdef IMGP_VAR

#include "Multigrid.h"
#include "constants.h"

! Set procedure pointer to the desired smoother implementation.
if (mgrid_smoother .eq. rbgs_smoother) then
   relaxerPtr => poisson_mg_relax
elseif (mgrid_smoother .eq. zebra_smoother) then
   relaxerPtr => poisson_mg_relax_ZEBRA
else
   NULLIFY(relaxerPtr)          !provoke an error
!!$   relaxerPtr => poisson_mg_relax ! .. or maybe rather fall back?
endif

select case (bcTypes(1))

!-------------------------------------------------------------------------------
  case (MG_BND_ISOLATED) ! isolated boundary conditions
#ifdef DEBUG_GRID
  print*,'isolated'
#endif


     bcTypes2(:) = MG_BND_DIRICHLET

     ! First get Dirichlet solution
     call gr_mgSolve (isrc, isoln, poisFact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                     bcTypes2, poisson_mg_solve, poisson_mg_residual, &
                     poisson_mg_residual_leafs,relaxerPtr)


    ! Construct the image mass distribution
    call gr_isoImageMass(isoln, imgw7)


    call Grid_getListOfBlocks(LEAF,blockList,blockCount)
    ! Compute the boundary values of the image mass potential
    do lb = 1,blockCount

        blockID = blockList(lb)

        ! Point to blocks center vars:
        call Grid_getBlkPtr(blockID,unkt,CENTER)
        unkt(imgw6,:,:,:) = 0.
        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,unkt,CENTER)

    enddo

    call gr_isoImageBdry(imgw7, imgw6, 1.)

    bcTypes2(:) = MG_BND_GIVENVAL

    ! Get the isolated potential of the image mass distribution
    call gr_mgSolve (imgw7, imgw6, 1., imgw1, imgw2, imgw3, imgw4, imgw5, &
                    bcTypes2, poisson_mg_solve, poisson_mg_residual, &
                    poisson_mg_residual_leafs,relaxerPtr)

    ! Finally subtract the image potential to obtain the isolated potential
    do lb = 1,blockCount

        blockID = blockList(lb)

        ! Point to blocks center vars:
        call Grid_getBlkPtr(blockID,unkt,CENTER)
        unkt(isoln,:,:,:) = unkt(isoln,:,:,:) - unkt(imgw6,:,:,:)
        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,unkt,CENTER)

    enddo

!-------------------------------------------------------------------------------

  case (MG_BND_PERIODIC) ! periodic boundary conditions
    call gr_mgSolve (isrc, isoln, poisFact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                    bcTypes, poisson_mg_solve, poisson_mg_residual, &
                    poisson_mg_residual_leafs,relaxerPtr)
!-------------------------------------------------------------------------------

  case (MG_BND_DIRICHLET) ! Dirichlet boundary conditions

    call gr_mgSolve (isrc, isoln, poisFact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                    bcTypes, poisson_mg_solve, poisson_mg_residual, &
                    poisson_mg_residual_leafs,relaxerPtr)

!-------------------------------------------------------------------------------

  case (MG_BND_NEUMANN) ! Neumann boundary conditions

    call gr_mgSolve (isrc, isoln, poisFact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                    bcTypes, poisson_mg_solve, poisson_mg_residual, &
                    poisson_mg_residual_leafs,relaxerPtr)

!-------------------------------------------------------------------------------

  case default
    call Driver_abortFlash("[poisson]  invalid boundary condition type!")


!-------------------------------------------------------------------------------

end select


#else

    call Timers_start("Multigrid_solve")

    call doSolve ( bcTypes )

    call Timers_stop("Multigrid_solve")

#endif


!===============================================================================

return

contains
  ! Auxiliary subroutine, introduced because procedure pointers are not available
  ! in compiliers that do not support this Fortran 2003 feature
  subroutine doSolve(bcTypesLoc)
    integer,intent(in), dimension(:) :: bcTypesLoc

   if (mgrid_smoother .eq. rbgs_smoother) then
      call gr_mgSolve (isrc, isoln, poisFact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                       bcTypesLoc, poisson_mg_solve, poisson_mg_residual, &
                       poisson_mg_residual_leafs,poisson_mg_relax)
   else
      call gr_mgSolve (isrc, isoln, poisFact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                       bcTypesLoc, poisson_mg_solve, poisson_mg_residual, &
                       poisson_mg_residual_leafs,poisson_mg_relax_ZEBRA)
   end if
  end subroutine DoSolve
end
