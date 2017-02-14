!!****if* source/Grid/GridSolvers/Multigrid/gr_hgSetZeroBoundary
!!
!! NAME
!!  gr_hgSetZeroBoundary
!!
!! SYNOPSIS
!!
!!  gr_hgSetZeroBoundary(integer, intent(in) :: level,
!!                       integer, intent(in) :: gr_iSoln)
!! DESCRIPTION
!!  
!!  Set external values on the coarsest grid level for the Poisson solver.
!!  Currently this only effects given-value and dirichlet conditions, and
!!  only works with the coarsest-level block or the chosen level.
!!
!! ARGUMENTS
!! 
!!  level      - the level to fill, or 0 to fill all levels
!!  gr_iSoln - the grid variable to be filled
!!
!! NOTE
!!  This routine used to have an additional argument "bndtype", whose definition
!!  was
!!               bndtype     If 0, do regardless of boundary condition;
!!                           if 1, do only if Dirichlet/given value boundaries
!!  However!  This routine was only ever called with bndtype=0 so it was removed
!!
!!***

!!REORDER(5): unk

#include "Flash.h"
#include "constants.h"
#include "Multigrid.h"

subroutine gr_hgSetZeroBoundary(level, gr_iSoln)

!==============================================================================
  use physicaldata, ONLY : unk
  use tree, ONLY : lnblocks,neigh,lrefine
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use gr_hgdata, ONLY:  gr_hgBndTypes
 
  implicit none
  
  integer, intent(in) :: gr_iSoln, level
  
  integer             :: b, i, j, k

!==============================================================================
  
  ! This routine doesn't apply for periodic or Neuman boundary conditions.
  ! DEV What routine DOES?
  if (ALL((gr_hgBndTypes(1:2*NDIM) .NE. MG_BND_DIRICHLET)  &
     &   .AND. (gr_hgBndTypes(1:2*NDIM) .NE. MG_BND_GIVENVAL))) return

  call Timers_start("gr_hgSetZeroBoundary")
  
  do b = 1, lnblocks
     if ((level == 0) .or. (lrefine(b) == level)) then

        if (gr_hgBndTypes(1)==MG_BND_DIRICHLET.OR.gr_hgBndTypes(1)==MG_BND_GIVENVAL) then
           if (neigh(1,ILO_FACE,b) <= PARAMESH_PHYSICAL_BOUNDARY) then
              unk(gr_iSoln,NGUARD,:,:,b) = 0.
           endif
        end if
        if (gr_hgBndTypes(2)==MG_BND_DIRICHLET.OR.gr_hgBndTypes(2)==MG_BND_GIVENVAL) then
           if (neigh(1,IHI_FACE,b) <= PARAMESH_PHYSICAL_BOUNDARY) then
              unk(gr_iSoln,NGUARD+NXB+1,:,:,b) = 0.
           endif
        endif

        if (NDIM >= 2) then
           if (gr_hgBndTypes(3)==MG_BND_DIRICHLET.OR.gr_hgBndTypes(3)==MG_BND_GIVENVAL) then
              if (neigh(1,JLO_FACE,b) <= PARAMESH_PHYSICAL_BOUNDARY) then
                 unk(gr_iSoln,:,((NGUARD-1)*K2D)+1,:,b) = 0.
              endif
           endif
           if (gr_hgBndTypes(4)==MG_BND_DIRICHLET.OR.gr_hgBndTypes(4)==MG_BND_GIVENVAL) then
              if (neigh(1,JHI_FACE,b) <= PARAMESH_PHYSICAL_BOUNDARY) then
                 unk(gr_iSoln,:,((NGUARD+NYB)*K2D)+1,:,b) = 0.
              endif
           endif
        endif

        if (NDIM == 3) then
           if (gr_hgBndTypes(5)==MG_BND_DIRICHLET.OR.gr_hgBndTypes(5)==MG_BND_GIVENVAL) then
              if (neigh(1,KLO_FACE,b) <= PARAMESH_PHYSICAL_BOUNDARY) then
                 unk(gr_iSoln,:,:,((NGUARD-1)*K3D)+1,b) = 0.
              endif
           endif
           if (gr_hgBndTypes(6)==MG_BND_DIRICHLET.OR.gr_hgBndTypes(6)==MG_BND_GIVENVAL) then
              if (neigh(1,KHI_FACE,b) <= PARAMESH_PHYSICAL_BOUNDARY) then
                 unk(gr_iSoln,:,:,((NGUARD+NZB)*K3D)+1,b) = 0.
              endif
           endif
        endif
     endif

  enddo

  call Timers_stop("gr_hgSetZeroBoundary")
  
  !=====================================================================
  
  return
end subroutine gr_hgSetZeroBoundary
