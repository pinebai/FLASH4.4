!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoMpoleMoments
!!
!! NAME
!!
!!  gr_isoMpoleMoments
!!
!! 
!! SYNOPSIS
!!
!!  gr_isoMpoleMoments(integer, intent(IN) : idensvar)
!!
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution,
!!  assuming the center of mass and the total mass have first been
!!  computed by gr_isoFindMassCenter.  On output, the Moment()
!!  array declared in gr_isoMpoleData contains the moments.
!!  
!! NOTES
!!  Special version for isolated boundary image mass.
!!
!!***

!!REORDER(4): solnData

subroutine gr_isoMpoleMoments (idensvar)

  !==============================================================

  use gr_isompoleData, ONLY:  mpole_geometry, mpole_quadrant,&
                         & qmax,mpole_INNER, mpole_OUTER, mpole_EVEN, mpole_ODD,&
                         & mpole_lmax,mpole_mmax,&
                         & fourpi,Xcm,Ycm,Zcm, Moment, Momtmp, Leg_fact
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords, Grid_getBlkBC, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkBoundBox, Grid_getDeltas
  use Grid_data, ONLY : gr_meshComm

  implicit none
  
#include "constants.h"
#include "Flash.h"
  include "Flash_mpi.h"
  
  
  integer, intent(IN)   :: idensvar
  
  integer   :: q, i, j, k, l, m, error
  real      :: zonemass, dvol, delx, dely, delz, x, y, z
  
  integer   :: lb, ii, jj, kk, Nint = 2
  real      :: xx, yy, zz
  real      :: delxx, delyy, delzz
  real      :: zonedens
  
  real               :: blockSize(MDIM)
  integer            :: nbrs(LOW:HIGH,MDIM)

  integer                       :: blockCount
  integer, dimension(MAXBLOCKS) :: blockList
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLImitsGC
  real, dimension(2, MDIM)          :: blkCoords
  real, dimension(MDIM)             :: delta
  real, pointer, dimension(:,:,:,:) :: solnData
  

  real, save         :: cylfactor

  !               Sum quantities over all locally held leaf blocks.
  
  Moment(:,:,:,:,:) = 0.
  
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  


  do lb = 1, blockCount
     !               Compute dimensions of each zone and subzone.
     
     call Grid_getDeltas(blockList(lb), delta)
     delx = delta(IAXIS)
     dely = delta(JAXIS)
     delz = delta(KAXIS)
     
     delxx = delx / real(Nint)  ! Nint is parameter of 2
     delyy = dely / real(Nint)
     delzz = delz / real(Nint)
     
     !               Get size information for this block.
     
     call Grid_getBlkPhysicalSize(blockList(lb), blockSize)

     !NO. NOT WHAT TO DO. BAD.  VERY BAD *SMACK*
     !call Grid_getBlkCenterCoords(blockList(lb), blockCenter)


     call Grid_getBlkBoundBox(blockList(lb), blkCoords)

     !               Get pointer to solution data.
     call Grid_getBlkPtr(blockList(lb),solnData, CENTER)
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     ! Find neighbors
     call Grid_getBlkBC(blockList(lb), nbrs)
     
     !               Compute the moment contributions for this block.
     
     select case (mpole_geometry)

     case (CARTESIAN)
        if(NDIM==3) then
           if (nbrs(LOW,IAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(LOW,IAXIS) - Xcm
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 zz = blkCoords(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz - Zcm
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    yy = blkCoords(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely - Ycm
                    zonemass = solnData(idensvar,blkLimits(LOW,IAXIS)-1,j,k) * dely * delz
                    call gr_isoZoneMoments (xx, yy, zz, zonemass)
                 enddo
              enddo
           endif
           
           if (nbrs(HIGH,IAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(HIGH,IAXIS) - Xcm
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 zz = blkCoords(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz - Zcm
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    yy = blkCoords(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely - Ycm
                    zonemass = solnData(idensvar,blkLimits(HIGH,IAXIS)+1,j,k) * dely * delz
                    call gr_isoZoneMoments (xx, yy, zz, zonemass)
                 enddo
              enddo
           endif
           
           if (nbrs(LOW,JAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(LOW,JAXIS) - Ycm
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 zz = blkCoords(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz - Zcm
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    xx = blkCoords(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx - Xcm
                    zonemass = solnData(idensvar,i,blkLimits(LOW,JAXIS)-1,k) * delx * delz
                    call gr_isoZoneMoments (xx, yy, zz, zonemass)
                 enddo
              enddo
           endif
           
           if (nbrs(HIGH,JAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(HIGH,JAXIS) - Ycm
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 zz = blkCoords(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz - Zcm
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    xx = blkCoords(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx - Xcm
                    zonemass = solnData(idensvar,i,blkLimits(HIGH,JAXIS)+1,k) * delx * delz
                    call gr_isoZoneMoments (xx, yy, zz, zonemass)
                 enddo
              enddo
           endif
           
           if (nbrs(LOW,KAXIS) /= NOT_BOUNDARY) then
              zz = blkCoords(LOW,KAXIS) - Zcm
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 yy = blkCoords(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely - Ycm
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    xx = blkCoords(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx - Xcm
                    zonemass = solnData(idensvar,i,j,blkLimits(LOW,KAXIS)-1) * delx * dely
                    call gr_isoZoneMoments (xx, yy, zz, zonemass)
                 enddo
              enddo
           endif
           
           if (nbrs(HIGH,KAXIS) /= NOT_BOUNDARY) then
              zz = blkCoords(HIGH,KAXIS) - Zcm
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 yy = blkCoords(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely - Ycm
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    xx = blkCoords(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx - Xcm
                    zonemass = solnData(idensvar,i,j,blkLimits(HIGH,KAXIS)+1) * delx * dely
                    call gr_isoZoneMoments (xx, yy, zz, zonemass)
                 enddo
              enddo
           endif
        else
           call Driver_abortFlash("Isobnd_mpole: cartesian work with 3D only")
        end if
     case (CYLINDRICAL)
        if(NDIM==2) then
           if (nbrs(LOW,JAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(LOW,IAXIS) + blockSize(IAXIS) - Xcm
              zz = 0.
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 yy = blkCoords(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely - Ycm
                 zonemass = solnData(idensvar,blkLimits(HIGH,IAXIS)+1,j,1) * &
                      cylfactor * xx * dely
                 call gr_isoZoneMoments (xx, 0., yy, zonemass)
              enddo
           endif
           
           if (nbrs(LOW,IAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(LOW,JAXIS) - Ycm
              zz = 0.
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xx = blkCoords(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx - Xcm
                 zonemass = solnData(idensvar,i,blkLimits(LOW,JAXIS)-1,1) * &
                      cylfactor * xx * delx
                 call gr_isoZoneMoments (xx, 0., yy, zonemass)
              enddo
           endif
           
           if (nbrs(HIGH,IAXIS) /= NOT_BOUNDARY) then
              yy = blkCoords(LOW,JAXIS) + blockSize(JAXIS) - Ycm
              zz = 0.
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xx = blkCoords(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx - Xcm
                 zonemass = solnData(idensvar,i,blkLimits(HIGH,JAXIS)+1,1) * &
                      cylfactor * xx * delx
                 call gr_isoZoneMoments (xx, 0., yy, zonemass)
              enddo
           endif
        else
           call Driver_abortFlash("Isobnd_mpole: cylindrical works 2D only")
        end if
           

     case (SPHERICAL)
        if(NDIM==1) then
           if (nbrs(HIGH,IAXIS) /= NOT_BOUNDARY) then
              xx = blkCoords(LOW,IAXIS) + blockSize(KAXIS) - Xcm
              yy = 0.
              zz = 0.
              zonemass = solnData(idensvar,blkLimits(HIGH,IAXIS)+1,1,1) * fourpi * xx**2
              call gr_isoZoneMoments (xx, 0., 0., zonemass)
           endif
        else
           call Driver_abortFlash("Isobnd_mpole: spherical works 1D only")
        end if

     end select

     call Grid_releaseBlkPtr(blockList(lb), solnData)
     
  enddo

  ! In zoneMoments, we only added the contribution of each zone to its
  ! particular radial bin.  Each zone should contribute its inner moment to
  ! all zones with radius greater than it, and its outer moment to all zones
  ! with radius less than it.

  do m = 0, mpole_lmax
     do l = m, mpole_lmax
        
        do i = MPOLE_EVEN, MPOLE_ODD
           
           do q = 2, qmax
              Moment(q,i,MPOLE_INNER,l,m) = Moment(q,i,MPOLE_INNER,l,m) + &
                   Moment(q-1,i,MPOLE_INNER,l,m)
           enddo
           
           do q = qmax-1, 1, -1
              Moment(q,i,MPOLE_OUTER,l,m) = Moment(q,i,MPOLE_OUTER,l,m) + &
                   Moment(q+1,i,MPOLE_OUTER,l,m)
           enddo
           
        enddo
        
     enddo
  enddo

  !=======================================================================

  !               Give all processors a copy of the global sums.

  do m = 0, mpole_lmax
     do l = m, mpole_lmax
        do j = MPOLE_INNER, MPOLE_OUTER
           do i = MPOLE_EVEN, MPOLE_ODD
              call MPI_AllReduce (Moment(0,i,j,l,m), Momtmp(0), qmax+1, & 
                   &                            FLASH_REAL, MPI_Sum,  & 
                   &                            gr_meshComm, error)
              Moment(:,i,j,l,m) = Momtmp(:)
           enddo
        enddo
     enddo
  enddo

  !               Normalize the moments properly.

  do l = 1, mpole_lmax
     do m = 1, l
        Moment(:,:,:,l,m) = Moment(:,:,:,l,m) * Leg_fact(l,m)
     enddo
  enddo

  !==========================================================================

  return
end subroutine gr_isoMpoleMoments
