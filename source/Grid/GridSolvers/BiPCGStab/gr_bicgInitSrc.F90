!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bicgInitSrc
!!
!! NAME
!!
!!  gr_bicgInitSrc
!!
!! SYNOPSIS
!!
!!  call gr_bicgInitSrc(integer(in) :: isrc,
!!                      integer(in) :: isoln,
!!                      real(in) :: poisfact)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isrc : 
!!
!!   isoln : 
!!
!!   poisfact : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     bicg_init_src()

!  Description: Initializes the right-hand side of the equation to be solved
!               with BiPCGStab.  This routine takes the density and sets up
!               the source term for the Poisson equation.  If periodic
!               boundaries are to be used, then the average is subtracted, so
!               that the source term averages to zero.


subroutine gr_bicgInitSrc (isrc, isoln, poisfact)

!===============================================================================
#include "Flash.h"

  use gr_bicgData

  use tree, only : nodetype,neigh

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getDeltas,       &
                                Grid_getBlkPhysicalSize

  use gr_interface, ONLY : gr_findMean

implicit none
#include "constants.h"
include 'mpif.h'

integer, intent(in) :: isrc, isoln
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

! call Grid_getLocalNumBlks(lnblocks)
call Grid_getListOfBlocks(LEAF,blockList,blockCount)

do lb = 1, blockCount
   blockID = blockList(lb)
   ! Point to blocks center vars:
   call Grid_getBlkPtr(blockID,solnData,CENTER)
   do k=kli,kui
      do j=jli,jui
         do i=ili,iui
            solnData(isrc,i,j,k) = solnData(isrc,i,j,k) * poisfact
         enddo
      enddo
   enddo
   ! Release pointers:
   call Grid_releaseBlkPtr(blockID,solnData,CENTER)
enddo


! For periodic or Neumann boundary conditions, we must subtract off the
! average value of the source.

!write(*,*) 'bicg_bnd_cond=',bicg_bnd_cond

if (bicg_bnd_cond == 0) then ! Substract source mean.

  call gr_findMean(isrc,2,.false.,sum)

  do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,solnData,CENTER)
      do k=kli,kui
         do j=jli,jui
            do i=ili,iui 
               solnData(isrc,i,j,k) = solnData(isrc,i,j,k) - sum
            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  enddo

! For given-value boundary conditions, we must subtract exterior boundary
! values (assumed to be stored in the first layer of boundary zones for
! the solution variable).  For cylindrical or spherical coordinates, r=0
! is always a Neumann boundary, and for 2D axisymmetric coordinates with
! one quadrant suppressed, z=0 is also a Neumann boundary.

else if (bicg_bnd_cond == 1) then

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


      if (nbrs(1) <= -20) then
         solnData(isrc,nguard+1,:,:) = &
              solnData(isrc,nguard+1,:,:) &
              - 2.*solnData(isoln,nguard,:,:)/delx**2
      endif


      if ((ndim >= 2) .and. (nbrs(3) <= -20)) then
         solnData(isrc,:,nguard+1,:) = &
              solnData(isrc,:,nguard+1,:) &
              - 2.*solnData(isoln,:,nguard,:)/dely**2
      endif


      if ((ndim == 3) .and. (nbrs(5) <= -20)) then
        solnData(isrc,:,:,nguard+1) = &
          solnData(isrc,:,:,nguard+1) &
          - 2.*solnData(isoln,:,:,nguard)/delz**2
      endif

      if (nbrs(2) <= -20) then
        solnData(isrc,nguard+nxb,:,:) = &
          solnData(isrc,nguard+nxb,:,:) &
          - 2.*solnData(isoln,nguard+nxb+1,:,:)/delx**2
      endif

      if ((ndim >= 2) .and. (nbrs(4) <= -20)) then
        solnData(isrc,:,nguard+nyb,:) = &
          solnData(isrc,:,nguard+nyb,:) &
          - 2.*solnData(isoln,:,nguard+nyb+1,:)/dely**2
      endif

      if ((ndim == 3) .and. (nbrs(6) <= -20)) then
        solnData(isrc,:,:,nguard+nzb) = &
          solnData(isrc,:,:,nguard+nzb) &
          - 2.*solnData(isoln,:,:,nguard+nzb+1)/delz**2
      endif

      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,solnData,CENTER)    

  enddo

endif

!===============================================================================

return
end subroutine gr_bicgInitSrc
