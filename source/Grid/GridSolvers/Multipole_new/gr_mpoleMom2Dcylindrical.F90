!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMom2Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleMom2Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom2Dcylindrical (integer (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Prepares for evaluation of the moments in 2D cylindrical geometry. In this
!!  routine, all the necessary arrays are prepared to enable evaluation of
!!  the moments in radial bin order. Each of the moments are grouped together
!!  according to their radial bins. This will ensure optimum unit stride values
!!  when accessing the big moment arrays and makes threading trivial.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleMom2Dcylindrical (idensvar)

  use Grid_interface,    ONLY : Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkBoundBox,    &
                                Grid_getDeltas,         &
                                Grid_getBlkIndexLimits

  use gr_mpoleInterface, ONLY : gr_mpoleMomBins2Dcylindrical

  use gr_mpoleData,      ONLY : gr_mpoleTwoPi,                  &
                                gr_mpoleDrInv,                  &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxQ,                   &
                                gr_mpoleMaxRadialZones,         &
                                gr_mpoleMinRadialZone,          &
                                gr_mpoleZoneRmax,               &
                                gr_mpoleZoneQmax,               &
                                gr_mpoleZoneType,               &
                                gr_mpoleZoneScalarInv,          &
                                gr_mpoleZoneLogNormInv,         &
                                gr_mpoleZoneExponentInv,        &
                                gr_mpoleInnerZoneMaxR,          &
                                gr_mpoleInnerZoneDrRadii,       &
                                gr_mpoleInnerZoneQlower,        &
                                gr_mpoleInnerZoneQupper,        &
                                gr_mpoleInnerZoneResolution,    &
                                gr_mpoleInnerZoneResolutionInv, &
                                gr_mpoleOuterZoneQshift,        &
                                gr_mpoleQ,                      &
                                gr_mpoleQused,                  &
                                gr_mpoleQnumberOfCells,         &
                                gr_mpoleQdataCells2D

  use gr_mpoleData,      ONLY : gr_mpoleRcenter,                &
                                gr_mpoleZcenter,                &
                                gr_mpoleBlockCount,             &
                                gr_mpoleBlockList

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: idensvar

  logical :: innerZonePotential

  integer :: blockNr, blockID
  integer :: DrUnit
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: maxCells
  integer :: nC, nQ
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: used
  integer :: zone

  integer, save :: maxQtype                ! for multithreading needs to be on stack (save)

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: bndBoxILow, bndBoxJLow
  real    :: cellDensity, cellMass, cellPlane, cellVolume
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf
  real    :: r, rlocal, rinDrs
  real    :: Rcyl
  real    :: sclInv, lgnInv, expInv
  real    :: z

  real    :: delta           (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real, pointer :: solnData (:,:,:,:)
!
!
!     ...The first pass over all blocks on the current processor will get us information
!        about how many different radial bin indices will be addressed and for each such
!        radial bin index, how many cells it will contain.
!
!
!$omp single
  gr_mpoleQused (:) = 0 

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
!
!
!          ...The 2D cylindrical case:
!
!
!                               ------
!                             /        \
!                            /     z    \
!                           |\     |    /|
!                           | \    |   / |
!                           |   ------   |
!                           |      |     |
!                           |      |     |         Rcyl --> stored in i-index (FLASH x)
!                           |      ----->|
!                           |       Rcyl |            z --> stored in j-index (FLASH y)
!                           |            |
!                           |            |
!                           |   ------   |
!                           | /        \ |
!                           |/          \|
!                            \          /
!                             \        /
!                               ------
!
!             Here the convention used in FLASH is to store the Rcyl values into the
!             i-index and the z values into the j-index.
!
!
     z = bndBoxJLow + DeltaJHalf - gr_mpoleZcenter
     do j = jmin,jmax
      Rcyl = bndBoxILow + DeltaIHalf - gr_mpoleRcenter
      do i = imin,imax

       r = sqrt (z * z + Rcyl * Rcyl)
!
!
!        ...Find the radial bin and increment radial bin counter.
!
!
       innerZonePotential = r <= gr_mpoleInnerZoneMaxR

       if (innerZonePotential) then

           rinDrs = r * gr_mpoleDrInnerZoneInv
           DrUnit = int (ceiling (rinDrs))
           Qlower = gr_mpoleInnerZoneQlower (DrUnit)
           Qupper = gr_mpoleInnerZoneQupper (DrUnit)

           do Q = Qlower,Qupper
              if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
           end do

       else

           do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones
              if (r - gr_mpoleZoneRmax (zone) <= ZERO) exit
           end do

           rlocal = r - gr_mpoleZoneRmax    (zone - 1)
           type   = gr_mpoleZoneType        (zone)
           sclInv = gr_mpoleZoneScalarInv   (zone)
           expInv = gr_mpoleZoneExponentInv (zone)

           if (type == ZONE_EXPONENTIAL) then
               Qlocal = ceiling ( (rlocal * sclInv * gr_mpoleDrInv) ** expInv )
           else if (type == ZONE_LOGARITHMIC) then
               lgnInv = gr_mpoleZoneLogNormInv (zone)
               Qlocal = ceiling ( expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE) )
           end if

           Q = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

       end if

       gr_mpoleQused (Q) = gr_mpoleQused (Q) + 1

       Rcyl = Rcyl + DeltaI
      end do
      z = z + DeltaJ
     end do

  end do
!
!
!     ...Create the arrays that will contain the radial info.
!
!
  maxQtype = count  (gr_mpoleQused /= 0)
  maxCells = maxval (gr_mpoleQused     )

  allocate (gr_mpoleQ              (1:maxQtype))
  allocate (gr_mpoleQnumberOfCells (1:maxQtype))
  allocate (gr_mpoleQdataCells2D   (1:maxCells , 1:maxQtype))
!
!
!     ...The second pass over all blocks on the current processor will scatter all
!        the radial bin information into the radial bin info array.
!
!
  gr_mpoleQused (:) = 0 

  nQ = 0

  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList (blockNr)

     call Grid_getBlkBoundBox     (blockID, bndBox)
     call Grid_getDeltas          (blockID, delta)
     call Grid_getBlkPtr          (blockID, solnData)
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
!
!
!          ...Create all the cell info needed and place into proper radial bin array places.
!             The cell volume is:
!
!                               pi * (R^2 - r^2) * Dz
!
!             where r is the left-most (smaller) and R is the right-most (larger)
!             radial cell distance and Dz is the cell's z-axis delta value. Since our
!             radial measure is based on the cell's center, we have: r = Rcyl - Dr/2
!             and R = Rcyl + Dr/2 with Dr being the cell's radial delta value.
!             Hence the cell volume becomes:
!
!                             2 * pi * Rcyl * Dr * Dz
!
!
     cellPlane = DeltaI * DeltaJ

     z = bndBoxJLow + DeltaJHalf - gr_mpoleZcenter
     do j = jmin,jmax
      Rcyl = bndBoxILow + DeltaIHalf - gr_mpoleRcenter
      do i = imin,imax

       cellVolume  = gr_mpoleTwoPi * Rcyl * cellPlane
       cellDensity = solnData (idensvar,i,j,1)
       cellMass    = cellDensity * cellVolume

       r = sqrt (z * z + Rcyl * Rcyl)
!
!
!        ...Find the radial bin.
!
!
       innerZonePotential = r <= gr_mpoleInnerZoneMaxR

       if (innerZonePotential) then

           rinDrs = r * gr_mpoleDrInnerZoneInv
           DrUnit = int (ceiling (rinDrs))
           Qlower = gr_mpoleInnerZoneQlower (DrUnit)
           Qupper = gr_mpoleInnerZoneQupper (DrUnit)

           do Q = Qlower,Qupper
              if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
           end do

       else

           do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones
              if (r - gr_mpoleZoneRmax (zone) <= ZERO) exit
           end do

           rlocal = r - gr_mpoleZoneRmax    (zone - 1)
           type   = gr_mpoleZoneType        (zone)
           sclInv = gr_mpoleZoneScalarInv   (zone)
           expInv = gr_mpoleZoneExponentInv (zone)

           if (type == ZONE_EXPONENTIAL) then
               Qlocal = ceiling ( (rlocal * sclInv * gr_mpoleDrInv) ** expInv )
           else if (type == ZONE_LOGARITHMIC) then
               lgnInv = gr_mpoleZoneLogNormInv (zone)
               Qlocal = ceiling ( expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE) )
           end if

           Q = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

       end if

       used = gr_mpoleQused (Q)

       if (used == 0) then

           nQ = nQ + 1

           gr_mpoleQused                (Q)             = nQ
           gr_mpoleQ                   (nQ)             = Q
           gr_mpoleQnumberOfCells      (nQ)             = 1
           gr_mpoleQdataCells2D      (1,nQ) % coord1    = z
           gr_mpoleQdataCells2D      (1,nQ) % cellMass  = cellMass
           gr_mpoleQdataCells2D      (1,nQ) % radius    = r

       else

           nC = gr_mpoleQnumberOfCells (used) + 1

           gr_mpoleQnumberOfCells    (used)             = nC
           gr_mpoleQdataCells2D   (nC,used) % coord1    = z
           gr_mpoleQdataCells2D   (nC,used) % cellMass  = cellMass
           gr_mpoleQdataCells2D   (nC,used) % radius    = r

       end if

       Rcyl = Rcyl + DeltaI
      end do
      z = z + DeltaJ
     end do

     call Grid_releaseBlkPtr (blockID, solnData)

  end do
!$omp end single
!
!
!    ...Call the radial bin clustered moment evaluation routine (all threads).
!
!
  call gr_mpoleMomBins2Dcylindrical (maxQtype)
!
!
!    ...Deallocate used arrays.
!
!
!$omp single
  deallocate (gr_mpoleQ             )
  deallocate (gr_mpoleQnumberOfCells)
  deallocate (gr_mpoleQdataCells2D  )
!$omp end single
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleMom2Dcylindrical
