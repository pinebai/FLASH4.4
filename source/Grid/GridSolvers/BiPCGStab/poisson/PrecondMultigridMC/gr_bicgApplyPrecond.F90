!!****if* source/Grid/GridSolvers/BiPCGStab/poisson/PrecondMultigridMC/gr_bicgApplyPrecond
!!
!! NAME
!!
!!  gr_bicgApplyPrecond
!!
!! SYNOPSIS
!!
!!  call gr_bicgApplyPrecond(integer, intent(IN)  :: isoln,
!!                           integer, intent(IN)  :: isrc,
!!                           integer, dimension(6), intent(IN)  :: bctypes,
!!                           real, dimension(2,6), intent(IN)  :: bcvalues)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isoln : 
!!
!!   isrc : 
!!
!!   bctypes : 
!!
!!   bcvalues : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     poisson()

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
!               bc_types(6)     Boundary condition types array:
!                                 MG_BND_PERIODIC
!                                 MG_BND_DIRICHLET
!                                 MG_BND_NEUMANN
!
!                                 index 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z  6 = +z
!               bc_values(2,6)  Values for dirichlet and neumann boundary
!                               conditions.  If Robins, then bc_values(1,*) holds
!                               dirichlet and bc_values(2,*) neumann conditions; if
!                               neumann or dirichlet, then bc_values(1,*) holds
!                               neumann or dirichlet values and bc_values(2,*)goes unused
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries
!               poisfact        Constant Poisson factor.  Used to scale the
!                               source function.  For example, for gravity,
!                               poisfact = 4*pi*G.


subroutine gr_bicgApplyPrecond (isoln, isrc, bcTypes, bcValues)

!===============================================================================
#include "Flash.h"

use gr_mgData

use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                              Grid_getListOfBlocks, &
                              Grid_getBlkPtr,       &
                              Grid_releaseBlkPtr,   &
                              GRID_PDE_BND_DIRICHLET, &
                              GRID_PDE_BND_NEUMANN,   &  
                              GRID_PDE_BND_PERIODIC

use gr_mgInterface,   ONLY: gr_mgSolve
!use gr_isoInterface, ONLY: gr_isoImageMass, gr_isoImageBdry
use Timers_interface, ONLY : Timers_start, Timers_stop
use Driver_interface, ONLY : Driver_abortFlash

use gr_bicgData, ONLY : gcellflg

implicit none
#include "Multigrid.h"

integer, intent(IN) :: isoln, isrc
integer, dimension(6), intent(IN) :: bcTypes
real, dimension(2,6), intent(IN) :: bcValues



integer       :: lb, lnblocks2, MyPE2, MasterPE2
integer, save :: imgw1, imgw2, imgw3, imgw4, imgw5, imgw6, imgw7
integer       :: i, j, k

!!$#ifdef IMGP_VAR
!!$integer :: blockCount
!!$integer :: blockList(MAXBLOCKS)
!!$#endif

  integer               :: eachboundary
  real                  :: poisfact 
  logical, save :: first_call = .true.

  integer, dimension(6), save :: bc_types

real, pointer, dimension(:,:,:,:) :: unkt

external poisson_mg_solve, poisson_mg_residual, poisson_mg_residual_leafs, &
         poisson_mg_relax

!===============================================================================

! Get key numbers from the database for the temporary variables we need.

imgw1 = MGW1_VAR
imgw2 = MGW2_VAR
imgw3 = MGW3_VAR
imgw4 = MGW4_VAR
imgw5 = MGW5_VAR
imgw6 = MGW6_VAR
imgw7 = MGW7_VAR

!===============================================================================

! Call the multigrid Poisson solver for different types of boundary conditions.
! Map GRID_PDE_ boundary conditions to MG_BND_ ones:
  if (first_call) then
  do eachboundary = 1,2*NDIM
   select case(bcTypes(eachBoundary))
     case (GRID_PDE_BND_PERIODIC)
        bc_types(eachBoundary)  = MG_BND_PERIODIC
     case (GRID_PDE_BND_NEUMANN)
        bc_types(eachBoundary)  = MG_BND_NEUMANN
     case (GRID_PDE_BND_DIRICHLET)
        bc_types(eachBoundary)  = MG_BND_DIRICHLET
     case default
           write(*,*) 'In MC Multigrid gr_bicgApplyPrecond: Boundary Condition not supported'
           call Driver_abortFlash('BCs unsupported')
     end select
  enddo

  first_call = .false.
  endif


poisfact = 1.0
call gr_mgSolve (isrc, isoln, poisfact, imgw1, imgw2, imgw3, imgw4, imgw5, &
                    bc_types, poisson_mg_solve, poisson_mg_residual, poisson_mg_residual_leafs, &
                    poisson_mg_relax)

! No need to fill guardcells before A*x operation in BiPCGStab
gcellflg = .false.


!!$#ifdef IMGP_VAR
!!$
!!$! This is a hack to make this version of poisson.F90's interface
!!$! consistent with 'blessed' version.
!!$
!!$
!!$select case (bc_types(1))
!!$
!!$!-------------------------------------------------------------------------------
!!$  print*,'isolated'
!!$  case (MG_BND_ISOLATED) ! isolated boundary conditions
!!$
!!$    
!!$     call Grid_getListOfBlocks(LEAF,blockList,blockCount)
!!$    
!!$    
!!$     bc_types2 = MG_BND_DIRICHLET
!!$
!!$     ! First get Dirichlet solution
!!$     call multigrid (isrc, isoln, poisfact, imgw1, imgw2, imgw3, imgw4, imgw5, &
!!$                     bc_types2, poisson_mg_solve, poisson_mg_residual, &
!!$                     poisson_mg_relax)
!!$
!!$
!!$    ! Construct the image mass distribution 
!!$    call call gr_isoImageMass(isoln, imgw7)
!!$
!!$
!!$    call Grid_getListOfBlocks(LEAF,blockList,blockCount)
!!$    ! Compute the boundary values of the image mass potential
!!$    do lb = 1,blockCount
!!$
!!$        blockID = blockList(lb)
!!$
!!$        ! Point to blocks center vars:
!!$        call Grid_getBlkPtr(blockID,unkt,CENTER)
!!$        unkt(imgw6,:,:,:) = 0.
!!$        ! Release pointers:
!!$        call Grid_releaseBlkPtr(blockID,unkt,CENTER)
!!$
!!$    enddo
!!$
!!$    call gr_isoImageBdry(imgw7, imgw6, 1.)
!!$
!!$    bc_types2 = MG_BND_GIVENVAL
!!$
!!$    ! Get the isolated potential of the image mass distribution
!!$    call multigrid (imgw7, imgw6, 1., imgw1, imgw2, imgw3, imgw4, imgw5, &
!!$                    bc_types2, poisson_mg_solve, poisson_mg_residual, &
!!$                    poisson_mg_relax)
!!$
!!$    ! Finally subtract the image potential to obtain the isolated potential
!!$    do lb = 1,blockCount
!!$
!!$        blockID = blockList(lb)
!!$
!!$        ! Point to blocks center vars:
!!$        call Grid_getBlkPtr(blockID,unkt,CENTER)
!!$        unkt(isoln,:,:,:) = unkt(isoln,:,:,:) - unkt(imgw6,:,:,:)
!!$        ! Release pointers:
!!$        call Grid_releaseBlkPtr(blockID,unkt,CENTER)
!!$
!!$    enddo
!!$
!!$!-------------------------------------------------------------------------------
!!$
!!$  case (MG_BND_PERIODIC) ! periodic boundary conditions
!!$    call multigrid (isrc, isoln, poisfact, imgw1, imgw2, imgw3, imgw4, imgw5, &
!!$                    bc_types, poisson_mg_solve, poisson_mg_residual, &
!!$                    poisson_mg_relax)
!!$!-------------------------------------------------------------------------------
!!$
!!$  case (MG_BND_DIRICHLET) ! Dirichlet boundary conditions
!!$
!!$    call multigrid (isrc, isoln, poisfact, imgw1, imgw2, imgw3, imgw4, imgw5, &
!!$                    bc_types, poisson_mg_solve, poisson_mg_residual, &
!!$                    poisson_mg_relax)
!!$
!!$!-------------------------------------------------------------------------------
!!$
!!$  case (MG_BND_NEUMANN) ! Neumann boundary conditions
!!$
!!$    call multigrid (isrc, isoln, poisfact, imgw1, imgw2, imgw3, imgw4, imgw5, &
!!$                    bc_types, poisson_mg_solve, poisson_mg_residual, &
!!$                    poisson_mg_relax)
!!$
!!$!-------------------------------------------------------------------------------
!!$
!!$  case default
!!$    call Driver_abortFlash("[poisson]  invalid boundary condition type!")
!!$
!!$
!!$!-------------------------------------------------------------------------------
!!$
!!$end select
!!$
!!$
!!$#else
!!$
!!$    call Timers_start("Multigrid_solve")
!!$
!!$    call multigrid (isrc, isoln, poisfact, imgw1, imgw2, imgw3, imgw4, imgw5, &
!!$                    bc_types, poisson_mg_solve, poisson_mg_residual, &
!!$                    poisson_mg_relax)
!!$
!!$    call Timers_stop("Multigrid_solve")
!!$    
!!$#endif
!!$
!!$
!!$!===============================================================================

return
end subroutine gr_bicgApplyPrecond
