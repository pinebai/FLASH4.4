!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleRad2Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleRad2Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleRad2Dcylindrical ()
!!
!! DESCRIPTION
!!
!!  This routine determines the radial sampling for accumulating the moments
!!  for the two-dimensional (2D) cylindrical case.
!!
!!***

subroutine gr_mpoleRad2Dcylindrical ()

  use Driver_interface,  ONLY : Driver_abortFlash

  use Grid_data,         ONLY : gr_meshMe,   &
                                gr_meshComm

  use Grid_interface,    ONLY : Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkBoundBox,    &
                                Grid_getDeltas,         &
                                Grid_getBlkIndexLimits, &
                                Grid_getMinCellSizes

  use gr_mpoleInterface, ONLY : gr_mpoleSetInnerZoneGrid, &
                                gr_mpoleSetOuterZoneGrid

  use gr_mpoleData,      ONLY : gr_mpoleDr,              &
                                gr_mpoleDrInv,           &
                                gr_mpoleDrInnerZone,     &
                                gr_mpoleDrInnerZoneInv,  &
                                gr_mpoleMaxR,            &
                                gr_mpoleIgnoreInnerZone, &
                                gr_mpoleInnerZoneExists, &
                                gr_mpoleInnerZoneMaxR,   &
                                gr_mpoleInnerZoneQmax,   &
                                gr_mpoleInnerZoneSize,   &
                                gr_mpoleOuterZoneExists, &
                                gr_mpoleDomainRmin,      &
                                gr_mpoleDomainRmax,      &
                                gr_mpoleDomainZmin,      &
                                gr_mpoleDomainZmax,      &
                                gr_mpoleRcenter,         &
                                gr_mpoleZcenter,         &
                                gr_mpoleBlockCount,      &
                                gr_mpoleBlockList

  implicit none

#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  include "Flash_mpi.h"

  integer :: blockID
  integer :: blockNr
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: nPinnerZone
  integer :: nRinnerZone
  integer :: nBlocal
  integer :: nPlocal
  integer :: nRlocal
  integer :: nRlocalPrev

  integer :: localData   (1:2)
  integer :: globalData  (1:2)
  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: DeltaI
  real    :: DeltaJ
  real    :: DeltaIHalf
  real    :: DeltaJHalf
  real    :: Dxmin,Dymin
  real    :: maxRsqr
  real    :: r1,r2,r3,r4
  real    :: Rcyl
  real    :: rsqr
  real    :: z,r

  real    :: delta        (1:MDIM)
  real    :: minCellSizes (1:MDIM)
  real    :: bndBox       (LOW:HIGH,1:MDIM)

  integer, allocatable :: blockListInnerZone (:)
  real,    allocatable :: RinnerZone         (:)
!
!
!       ...Get the minimum cell sizes for the entire domain.
!
!
  call Grid_getMinCellSizes (minCellSizes)

  Dxmin = minCellSizes (IAXIS)
  Dymin = minCellSizes (JAXIS)
!
!
!       ...Determine the maximum distance from the center of multipole expansion
!          to the computational domain boundaries.
!
!          Determine also the 'atomic' radial spacing. This is the smallest
!          possible size of one radial bin. Half the Geometric mean (n-th root
!          of product of n samples) is used to determine the atomic spacing from
!          the minimum cell spacings in the radial and z direction. The inverse is
!          also calculated for further reference.
!
!
  r1 =   (gr_mpoleDomainRmin - gr_mpoleRcenter) * (gr_mpoleDomainRmin - gr_mpoleRcenter) &
       + (gr_mpoleDomainZmin - gr_mpoleZcenter) * (gr_mpoleDomainZmin - gr_mpoleZcenter)

  r2 =   (gr_mpoleDomainRmax - gr_mpoleRcenter) * (gr_mpoleDomainRmax - gr_mpoleRcenter) &
       + (gr_mpoleDomainZmin - gr_mpoleZcenter) * (gr_mpoleDomainZmin - gr_mpoleZcenter)

  r3 =   (gr_mpoleDomainRmin - gr_mpoleRcenter) * (gr_mpoleDomainRmin - gr_mpoleRcenter) &
       + (gr_mpoleDomainZmax - gr_mpoleZcenter) * (gr_mpoleDomainZmax - gr_mpoleZcenter)

  r4 =   (gr_mpoleDomainRmax - gr_mpoleRcenter) * (gr_mpoleDomainRmax - gr_mpoleRcenter) &
       + (gr_mpoleDomainZmax - gr_mpoleZcenter) * (gr_mpoleDomainZmax - gr_mpoleZcenter)

  maxRsqr       = max  (r1,r2,r3,r4)
  gr_mpoleMaxR  = sqrt (maxRsqr)
  gr_mpoleDr    = HALF * (Dxmin * Dymin) ** HALF
  gr_mpoleDrInv = ONE / gr_mpoleDr
!
!
!       ...Set initial indicators for inner and outer zone.
!
!
  gr_mpoleInnerZoneExists = .not. gr_mpoleIgnoreInnerZone
  gr_mpoleOuterZoneExists = .not. gr_mpoleInnerZoneExists

  if (gr_mpoleInnerZoneExists) then
!
!
!     ...Proceed with establishing the inner zone (if any). From the determined
!        inner zone atomic length and the previously found maximal radial domain
!        distance, readjust the inner zone size variable. Two cases can happen:
!        1) the size of the inner zone fits into the complete radial doamin (no
!        adjustment needed) or 2) the size of the inner zone exceeds the complete
!        radial domain (adjustment needed). Also override existence criterion for
!        the outer zone, if the largest domain radius exceeds the inner zone radius.
!
!
      gr_mpoleOuterZoneExists = (gr_mpoleInnerZoneSize * gr_mpoleDrInnerZone < gr_mpoleMaxR)

      if ( gr_mpoleInnerZoneSize * gr_mpoleDrInnerZone > gr_mpoleMaxR ) then
           gr_mpoleInnerZoneSize = int (ceiling (gr_mpoleMaxR * gr_mpoleDrInnerZoneInv))
      end if
!
!
!     ...Determine the number of radii to be expected in the inner zone.
!        For each processor, store those local blockID's that actually
!        have radii in the inner zone.
!
!
      allocate (blockListInnerZone (1:MAXBLOCKS))

      gr_mpoleInnerZoneMaxR = real (gr_mpoleInnerZoneSize) * gr_mpoleDrInnerZone
      maxRsqr               = gr_mpoleInnerZoneMaxR * gr_mpoleInnerZoneMaxR

      nBlocal = 0
      nRlocal = 0
      nRlocalPrev = 0

      do blockNr = 1,gr_mpoleBlockCount

         blockID = gr_mpoleBlockList (blockNr)

         call Grid_getBlkBoundBox     (blockID, bndBox)
         call Grid_getDeltas          (blockID, delta)
         call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

         imin       = blkLimits (LOW, IAXIS)
         jmin       = blkLimits (LOW, JAXIS)
         imax       = blkLimits (HIGH,IAXIS)
         jmax       = blkLimits (HIGH,JAXIS)

         DeltaI     = delta (IAXIS)
         DeltaJ     = delta (JAXIS)
         DeltaIHalf = DeltaI * HALF
         DeltaJHalf = DeltaJ * HALF

         bndBoxILow = bndBox (LOW,IAXIS)
         bndBoxJLow = bndBox (LOW,JAXIS)

         z = bndBoxJLow + DeltaJHalf - gr_mpoleZcenter
         do j = jmin,jmax
            Rcyl = bndBoxILow + DeltaIHalf - gr_mpoleRcenter
            do i = imin,imax

               rsqr = z * z + Rcyl * Rcyl

               if (rsqr <= maxRsqr) then
                   nRlocal = nRlocal + 1
               end if

               Rcyl = Rcyl + DeltaI
            end do
            z = z + DeltaJ
         end do

         if (nRlocal > nRlocalPrev) then
             nBlocal = nBlocal + 1
             blockListInnerZone (nBlocal) = blockID
         end if

         nRlocalPrev = nRlocal

      end do
!
!
!     ...Calculate the total number of processors contributing to the inner
!        zone radii and the overall total number of inner zone radii to be
!        expected. Allocate the array that will contain all inner zone radii
!        on all processors. If no inner zone radii are found globally, there
!        is something wrong and the program has to stop.
!
!
      nPlocal = min (nRlocal,1)  ! current processor adds +1, if inner zone radii found

      localData (1) = nPlocal
      localData (2) = nRlocal

      call MPI_AllReduce (localData,     &
                          globalData,    &
                          2,             &
                          FLASH_INTEGER, &
                          MPI_SUM,       &
                          gr_meshComm,   &
                          error          )

      nPinnerZone = globalData (1)
      nRinnerZone = globalData (2)

      if (nRinnerZone == 0) then
          call Driver_abortFlash ('[gr_mpoleRad2Dcylindrical] ERROR: no inner zone radii found')
      end if

      allocate (RinnerZone (1:nRinnerZone))
!
!
!     ...Calculate and store now all inner zone radii on each processor.
!        Loop only over those local blocks which actually contribute to the
!        inner zone (skip, if no blocks).
!
!
      nRlocal = 0

      if (nBlocal > 0) then

          do blockNr = 1,nBlocal

             blockID = blockListInnerZone (blockNr)

             call Grid_getBlkBoundBox     (blockID, bndBox)
             call Grid_getDeltas          (blockID, delta)
             call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

             imin       = blkLimits (LOW, IAXIS)
             jmin       = blkLimits (LOW, JAXIS)
             imax       = blkLimits (HIGH,IAXIS)
             jmax       = blkLimits (HIGH,JAXIS)

             DeltaI     = delta (IAXIS)
             DeltaJ     = delta (JAXIS)
             DeltaIHalf = DeltaI * HALF
             DeltaJHalf = DeltaJ * HALF

             bndBoxILow = bndBox (LOW,IAXIS)
             bndBoxJLow = bndBox (LOW,JAXIS)

             z = bndBoxJLow + DeltaJHalf - gr_mpoleZcenter
             do j = jmin,jmax
                Rcyl = bndBoxILow + DeltaIHalf - gr_mpoleRcenter
                do i = imin,imax

                   r = sqrt (z * z + Rcyl * Rcyl)

                   if (r <= gr_mpoleInnerZoneMaxR) then
                       nRlocal = nRlocal + 1
                       RinnerZone (nRlocal) = r
                   end if

                   Rcyl = Rcyl + DeltaI
                end do
                z = z + DeltaJ
             end do

          end do
      end if

      deallocate (blockListInnerZone)
!
!
!       ...Set up the inner zone radial grid.
!
!
      call gr_mpoleSetInnerZoneGrid (nRlocal,     &
                                     nRinnerZone, &
                                     nPinnerZone, &
                                     RinnerZone   )

      deallocate (RinnerZone)

  else
!
!
!       ...No inner zone! Set the inner zone variables to nonexistent.
!
!
      gr_mpoleDrInnerZone   = ZERO
      gr_mpoleInnerZoneMaxR = ZERO
      gr_mpoleInnerZoneQmax = 0

  end if  ! inner zone condition
!
!
!       ...Complete the radial grid picture by setting up the outer (statistical)
!          zone radial grid.
!
!
  call gr_mpoleSetOuterZoneGrid ()
!
!
!       ...Ready!
!
!
  return
end subroutine gr_mpoleRad2Dcylindrical
