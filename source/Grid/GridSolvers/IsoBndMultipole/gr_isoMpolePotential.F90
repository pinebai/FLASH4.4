!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoMpolePotential
!!
!! NAME
!!
!!  gr_isoMpolePotential
!!
!!
!! SYNOPSIS
!!
!!  call gr_isoMpolePotential(integer(IN) :: ipotvar,
!!                               real(IN) :: poisfact)
!!
!! DESCRIPTION
!!
!! Computes the potential field due to a density field associated
!!               with a given set of multipole moments.  The moments are taken
!!               from the mpole_common variable Moment().  On output, the
!!               variable indexed by ipotvar contains the potential.  This
!!               calculation is entirely local to each processor, as each
!!               processor now has a separate copy of the moments.
!!
!!               Special version for isolated boundary image mass.
!!
!! ARGUMENTS
!!
!!  ipotvar -- output grid variable index (from Flash.h) to contain the potential
!!  poisfact -- scaling factor for the eventual solution (4*pi*newton's constant)
!!
!! NOTES
!!
!!***

!!REORDER(4):solnData

subroutine gr_isoMpolePotential (ipotvar, poisfact)

  use gr_isoMpoleData, ONLY: fourpi_inv, mpole_geometry, Xcm, Ycm, Zcm
  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getDeltas, &
    Grid_getBlkBoundBox, Grid_getBlkPtr, Grid_getBlkBC, &
    Grid_getBlkPhysicalSize, Grid_releaseBlkPtr, Grid_getListOfBlocks, Grid_getBlkIndexLimits

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)   :: ipotvar
  real, intent(IN)      :: poisfact

  integer   :: i, j, k, lb, blockCount, blockList(MAXBLOCKS)
  real      :: potential, delx, dely, delz
  real      :: mpfactor

  integer   :: imin, imax, jmin, jmax, kmin, kmax
  real      :: xx, yy, zz

  real      :: blksize(MDIM), blkCoords(2,MDIM)
  integer, dimension(LOW:HIGH,MDIM)::nbrs,blkLimits,blkLimitsGC

  real, pointer      :: solnData(:,:,:,:)

  real,dimension(MDIM) :: delta
  !======================================================================

  ! Note In general purpose gravity version, this section only done for geometry/=2dspherical
  !   There is a different solution for the 2d spherical case.
   
  !               Scale the Poisson source term factor appropriately.
  mpfactor = poisfact * fourpi_inv

  !               Compute potential on all locally held blocks.


  !!  call Grid_getLocalNumBlks(lnBlocks) !-- ummm... they might not be in order!
  !!     and what about getting the appropriate deltas?
  !!  this would have worked if all blocks needed to be addressed.  If only LEAF blocks, change
  !!  the next line to LEAF from ALL_BLKS.  Note that PRicker's Flash2 version uses all blocks
  call Grid_getListOfBlocks(ALL_BLKS, blockList, blockCount)

  do lb = 1, blockCount

     !               Compute dimensions of each zone and subzone.

     call Grid_getDeltas(blockList(lb),delta)
     call Grid_getBlkBoundBox(blockList(lb),blkCoords)
     delx = delta(IAXIS)
     dely = delta(JAXIS)
     delz = delta(KAXIS)

     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     imin = blkLimits(LOW, IAXIS)
     imax = blkLimits(HIGH, IAXIS)
     jmin = blkLimits(LOW, JAXIS)
     jmax = blkLimits(HIGH, JAXIS)
     kmin = blkLimits(LOW, KAXIS)
     kmax = blkLimits(HIGH, KAXIS)

     !more info about block 
     call Grid_getBlkPtr(blockList(lb),solnData)
     call Grid_getBlkBC(blockList(lb),nbrs)
     call Grid_getBlkPhysicalSize(blockList(lb),blksize)  ! used for spherical case

     select case (mpole_geometry)

     case (CARTESIAN)
        if(NDIM==3) then
           if (nbrs(LOW,IAXIS) /= NOT_BOUNDARY) then
          ! Flash2 version has 
          ! blksize = dBaseBlockSize(lb)
          ! blkcoord = dBaseBlockCoord(lb) - 0.5*blksize
          ! xx = blkcoord(lb) - Xcm ! so this is the left boundary of the block
              xx = blkCoords(LOW,IAXIS) - Xcm
              do k = kmin, kmax
                 ! Flash2 zz= blkcoord(3) + (k-kmin+0.5)*delz - Zcm
                 zz = blkCoords(LOW,KAXIS) + (k-kmin+0.5)*delz - Zcm
                 do j = jmin, jmax
                    yy = blkCoords(LOW,JAXIS) + (j-jmin+0.5)*dely - Ycm
                    call gr_isoZonepotential (xx, yy, zz, potential)
                    solnData(ipotvar,blkLimits(LOW,IAXIS)-1,j,k) = -mpfactor*potential
! the general purpose case here computes by looping over subregions
!  Should be equivalent for constant volume 
!  do k
!     do j
!         ddvol = delx*dely*delz  (or Grid_getVolume should work)
!         dvol = dvol + ddvol
!         potsum = potsum + potential*ddvol
!     end do
!  end do
!  solnData(ipotvar,i,j,k) = -mpfactor*potsum/dvol
                 enddo
              enddo
           endif

           if (nbrs(HIGH,IAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(HIGH,IAXIS) - Xcm
              do k = kmin, kmax
                 zz = blkCoords(LOW,KAXIS) + (k-kmin+0.5)*delz - Zcm
                 do j = jmin, jmax
                    yy = blkCoords(LOW,JAXIS) + (j-jmin+0.5)*dely - Ycm
                    call gr_isoZonepotential (xx, yy, zz, potential)
                    solnData(ipotvar,blkLimits(HIGH,IAXIS)+1,j,k) = -mpfactor*potential
                 enddo
              enddo
           endif

           if (nbrs(LOW,JAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(LOW,JAXIS) - Ycm
              do k = kmin, kmax
                 zz = blkCoords(LOW,KAXIS) + (k-kmin+0.5)*delz - Zcm
                 do i = imin, imax
                    xx = blkCoords(LOW,IAXIS) + (i-imin+0.5)*delx - Xcm
                    call gr_isoZonepotential (xx, yy, zz, potential)
                    solnData(ipotvar,i,blkLimits(LOW,JAXIS)-1,k) = -mpfactor*potential
                 enddo
              enddo
           endif

           if (nbrs(HIGH,JAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(HIGH,JAXIS) - Ycm
              do k = kmin, kmax
                 zz = blkCoords(LOW,KAXIS) + (k-kmin+0.5)*delz - Zcm
                 do i = imin, imax
                    xx = blkCoords(LOW,IAXIS) + (i-imin+0.5)*delx - Xcm
                    call gr_isoZonepotential (xx, yy, zz, potential)
                    solnData(ipotvar,i,blkLimits(HIGH,JAXIS)+1,k) = -mpfactor*potential
                 enddo
              enddo
           endif

           if (nbrs(LOW,KAXIS) /= NOT_BOUNDARY) then
              zz = blkCoords(LOW,KAXIS) - Zcm
              do j = jmin, jmax
                 yy = blkCoords(LOW,JAXIS) + (j-jmin+0.5)*dely - Ycm
                 do i = imin, imax
                    xx = blkCoords(LOW,IAXIS) + (i-imin+0.5)*delx - Xcm
                    call gr_isoZonepotential (xx, yy, zz, potential)
                    solnData(ipotvar,i,j,blkLimits(LOW,KAXIS)-1) = -mpfactor*potential
                 enddo
              enddo
           endif

           if (nbrs(HIGH,KAXIS) /= NOT_BOUNDARY) then
              zz = blkCoords(HIGH,KAXIS) - Zcm
              do j = jmin, jmax
                 yy = blkCoords(LOW,JAXIS) + (j-jmin+0.5)*dely - Ycm
                 do i = imin, imax
                    xx = blkCoords(LOW,IAXIS) + (i-imin+0.5)*delx - Xcm
                    call gr_isoZonepotential (xx, yy, zz, potential)
                    solnData(ipotvar,i,j,blkLimits(HIGH,KAXIS)+1) = -mpfactor*potential
                 enddo
              enddo
           endif
        else
           call Driver_abortFlash("Isobnd_mpole: cartesian works 3D only")
        end if

     case (CYLINDRICAL)
        if(NDIM==2) then
           if (nbrs(HIGH,IAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(HIGH,IAXIS) - Xcm
              zz = 0.
              do j = jmin, jmax
                 yy = blkCoords(LOW,JAXIS) + (j-jmin+0.5)*dely - Ycm
                 call gr_isoZonepotential (xx, 0., yy, potential)
                 solnData(ipotvar,blkLimits(HIGH,IAXIS)+1,j,1) = -mpfactor*potential
              enddo
           endif

           if (nbrs(LOW,JAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(LOW,JAXIS) - Ycm
              zz = 0.
              do i = imin, imax
                 xx = blkCoords(LOW,IAXIS) + (i-imin+0.5)*delx - Xcm
                 call gr_isoZonepotential (xx, 0., yy, potential)
                 solnData(ipotvar,i,blkLimits(LOW,JAXIS)-1,1) = -mpfactor*potential
              enddo
           endif

           if (nbrs(HIGH,JAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(HIGH,JAXIS) - Ycm
              zz = 0.
              do i = imin, imax
                 xx = blkCoords(LOW,IAXIS) + (i-imin+0.5)*delx - Xcm
                 call gr_isoZonepotential (xx, 0., yy, potential)
                 solnData(ipotvar,i,blkLimits(HIGH,JAXIS)+1,1) = -mpfactor*potential
              enddo
           endif
        else
           call Driver_abortFlash("Isobnd_mpole: cylindrical works 2D only")
        end if

     case (SPHERICAL)
        if(NDIM==1) then
           if (nbrs(HIGH,IAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(LOW,IAXIS) + blksize(IAXIS) - Xcm
              yy = 0.
              zz = 0.
              call gr_isoZonepotential (xx, 0., 0., potential)
              solnData(ipotvar,blkLimits(HIGH,IAXIS)+1,1,1) = -mpfactor*potential
           endif
        else
           call Driver_abortFlash("Isobnd_mpole: spherical works 1D only")
        end if

     end select

     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo

  return
end subroutine gr_isoMpolePotential
