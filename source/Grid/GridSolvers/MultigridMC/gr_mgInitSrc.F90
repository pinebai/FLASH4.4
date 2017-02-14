!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgInitSrc
!!
!! NAME
!!
!!  gr_mgInitSrc
!!
!! SYNOPSIS
!!
!!  call gr_mgInitSrc(integer(in) :: isrc_dens,
!!                    real(in) :: poisfact,
!!                    integer(in) :: img_src,
!!                    integer(in) :: img_soln)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isrc_dens : 
!!
!!   poisfact : 
!!
!!   img_src : 
!!
!!   img_soln : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_init_src()

!  Description: Initializes the right-hand side of the equation to be solved
!               with multigrid.  This routine takes the density and sets up
!               the source term for the Poisson equation.  If periodic
!               boundaries are to be used, then the average is subtracted, so
!               that the source term averages to zero.


subroutine gr_mgInitSrc (isrc_dens, poisfact, img_src, img_soln)

!===============================================================================
#include "Flash.h"

use gr_mgData

  use tree, only : nodetype,neigh

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getDeltas,       &
                                Grid_getBlkPhysicalSize

  use gr_interface, ONLY : gr_findMean

implicit none
#include "Multigrid.h"
#include "constants.h"
include 'mpif.h'

integer, intent(in) :: isrc_dens, img_src, img_soln
real, intent(in)    :: poisfact

integer       :: lb, i, j, k, ierr, lnblocks
real          :: lvol, vol, lsum, bsum, sum
real          :: nbinv, cvol, bvol, delx, dely, delz

integer, parameter :: MAXDIMS = 3
real               :: size(MAXDIMS), nbrs(2*MAXDIMS)

real, pointer, dimension(:,:,:,:) :: solnData

integer :: blockCount,blockID
integer :: blockList(MAXBLOCKS)

integer, parameter :: nxb = NXB
integer, parameter :: nyb = NYB
integer, parameter :: nzb = NZB
integer, parameter :: nguard = NGUARD
integer, parameter :: ndim = NDIM
integer, parameter :: k2d = K2D
integer, parameter :: k3d = K3D

!===============================================================================


! Copy density into source array and multiply by Poisson source factor.
! We set values on all leaf-node blocks, then restrict to propagate the
! source term to all levels.

call Grid_getListOfBlocks(LEAF,blockList,blockCount)

do lb = 1, blockCount
   blockID = blockList(lb)
   ! Point to blocks center vars:
   call Grid_getBlkPtr(blockID,solnData,CENTER)
   solnData(img_src,:,:,:) = solnData(isrc_dens,:,:,:) * poisfact
   ! Release pointers:
   call Grid_releaseBlkPtr(blockID,solnData,CENTER)
enddo


! Restricting the source is not necessary for the Martin & Cartwright algorithm.
!!do i = mesh_lrefmax, 2, -1
!!  call mg_restrict (i, img_src, img_src)
!!enddo

! For periodic or Neumann boundary conditions, we must subtract off the
! average value of the source.

if (mg_bnd_cond == 0) then ! Substract source mean.

!!$  lvol = 0.
!!$  lsum = 0.
!!$  nbinv = 1. / real(nxb)
!!$  if (ndim >= 2) nbinv = nbinv / real(nyb)
!!$  if (ndim == 3) nbinv = nbinv / real(nzb)
!!$
!!$  do lb = 1, blockCount
!!$
!!$      blockID = blockList(lb)
!!$
!!$      ! Get BlockSize:
!!$      call Grid_getBlkPhysicalSize(blockID,size)
!!$      ! Point to blocks center vars:
!!$      call Grid_getBlkPtr(blockID,solnData,CENTER)      
!!$
!!$      bvol = size(1)
!!$      if (ndim >= 2) bvol = bvol * size(2)
!!$      if (ndim == 3) bvol = bvol * size(3)
!!$      cvol = bvol * nbinv
!!$      lvol = lvol + bvol
!!$      bsum = 0.
!!$      do k = kli, kui
!!$        do j = jli, jui
!!$          do i = ili, iui
!!$            bsum = bsum + solnData(img_src,i,j,k)
!!$          enddo
!!$        enddo
!!$      enddo
!!$      lsum = lsum + bsum * cvol
!!$
!!$      ! Release pointers:
!!$      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
!!$    
!!$  enddo
!!$
!!$  call mpi_allreduce ( lsum, sum, 1, MPI_DOUBLE_PRECISION, & 
!!$                       MPI_SUM, MPI_COMM_WORLD, ierr )
!!$  call mpi_allreduce ( lvol, vol, 1, MPI_DOUBLE_PRECISION, & 
!!$                       MPI_SUM, MPI_COMM_WORLD, ierr )
!!$
!!$  sum = sum / vol

  call gr_findMean(img_src,2,.false.,sum)

  do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,solnData,CENTER) 
      solnData(img_src,:,:,:) = solnData(img_src,:,:,:) - sum
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  enddo

! For given-value boundary conditions, we must subtract exterior boundary
! values (assumed to be stored in the first layer of boundary zones for
! the solution variable).  For cylindrical or spherical coordinates, r=0
! is always a Neumann boundary, and for 2D axisymmetric coordinates with
! one quadrant suppressed, z=0 is also a Neumann boundary.

else if (mg_bnd_cond == 1) then

  do lb = 1, blockCount

      blockID = blockList(lb) 

      nbrs(1:2) = neigh(:,IAXIS,blockID)
      nbrs(3:4) = neigh(:,JAXIS,blockID)
      nbrs(5:6) = neigh(:,KAXIS,blockID)

    
      ! Get BlockSize:
      call Grid_getBlkPhysicalSize(blockID,size)

      delx = size(1) / nxb
      if (ndim >= 2) dely = size(2) / nyb
      if (ndim == 3) delz = size(3) / nzb

      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,solnData,CENTER) 


      if ((mg_geometry == MG_GEOM_1DCARTESIAN) .or. &
          (mg_geometry == MG_GEOM_2DCARTESIAN) .or. &
          (mg_geometry == MG_GEOM_3DCARTESIAN)) then

        if (nbrs(1) <= -20) then
          solnData(img_src,nguard+1,:,:) = &
            solnData(img_src,nguard+1,:,:) &
           - 2.*solnData(img_soln,nguard,:,:)/delx**2
        endif

      endif

      if ((mg_geometry /= MG_GEOM_2DCYLAXISYM) .or. &
          (.not. quadrant)) then

        if ((ndim >= 2) .and. (nbrs(3) <= -20)) then
          solnData(img_src,:,nguard+1,:) = &
            solnData(img_src,:,nguard+1,:) &
           - 2.*solnData(img_soln,:,nguard,:)/dely**2
        endif

      endif

      if ((ndim == 3) .and. (nbrs(5) <= -20)) then
        solnData(img_src,:,:,nguard+1) = &
          solnData(img_src,:,:,nguard+1) &
          - 2.*solnData(img_soln,:,:,nguard)/delz**2
      endif

      if (nbrs(2) <= -20) then
        solnData(img_src,nguard+nxb,:,:) = &
          solnData(img_src,nguard+nxb,:,:) &
          - 2.*solnData(img_soln,nguard+nxb+1,:,:)/delx**2
      endif

      if ((ndim >= 2) .and. (nbrs(4) <= -20)) then
        solnData(img_src,:,nguard+nyb,:) = &
          solnData(img_src,:,nguard+nyb,:) &
          - 2.*solnData(img_soln,:,nguard+nyb+1,:)/dely**2
      endif

      if ((ndim == 3) .and. (nbrs(6) <= -20)) then
        solnData(img_src,:,:,nguard+nzb) = &
          solnData(img_src,:,:,nguard+nzb) &
          - 2.*solnData(img_soln,:,:,nguard+nzb+1)/delz**2
      endif

      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,solnData,CENTER)    

  enddo

endif

!===============================================================================

return
end
