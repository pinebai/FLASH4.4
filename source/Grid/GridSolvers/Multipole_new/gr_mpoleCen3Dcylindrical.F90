!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCen3Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleCen3Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen3Dcylindrical (integer, intent (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 3D cylindrical
!!  geometry. It computes the center of expansion for the multipoles for 3D cylindrical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                            integral (r * rho * rho  dr)
!!              Cen (x,y,z) = ----------------------------
!!                              integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                     sum cells (cell center r * cell mass * cell rho)
!!       Cen (x,y,z) = ------------------------------------------------
!!                             sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (x,y,z) has been determined, it is placed on the
!!  the nearest cell corner. The following is computed here:
!!
!!                  1) multipole expansion center (placed on nearest cell corner)
!!                  2) total mass (aborts, if <= 0)
!!                  3) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!! NOTES
!!
!!  gr_mpoleXcenter, gr_mpoleYcenter and gr_mpoleZcenter denote the location of
!!  the center of multipole expansion in the 3D cartesian framework (going from
!!  R,phi,z -> x,y,z).
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen3Dcylindrical (idensvar)

  use Grid_data,         ONLY : gr_meshMe,   &
                                gr_meshComm

  use Driver_interface,  ONLY : Driver_abortFlash

  use Grid_interface,    ONLY : Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkBoundBox,    &
                                Grid_getDeltas,         &
                                Grid_getBlkIndexLimits, &
                                Grid_getCellCoords

  use gr_mpoleData,      ONLY : gr_mpoleDomainRmax,     &
                                gr_mpoleDomainZmax,     &
                                gr_mpoleDomainPhiMax,   &
                                gr_mpoleDrInnerZone,    &
                                gr_mpoleDrInnerZoneInv, &
                                gr_mpolePi,             & 
                                gr_mpoleHalfPi,         & 
                                gr_mpoleTwoPi,          & 
                                gr_mpoleXcenter,        &
                                gr_mpoleYcenter,        &
                                gr_mpoleZcenter,        &
                                gr_mpoleRcenter,        &
                                gr_mpolePhiCenter,      &
                                gr_mpoleTotalMass,      &
                                gr_mpoleBlockCount,     &
                                gr_mpoleBlockList

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  include "Flash_mpi.h"
  
  integer, intent (in) :: idensvar

  logical :: domainRmax, domainZmax, domainPhiMax
  logical :: guardCells
  logical :: insideBlock
  logical :: invokeRecv
  logical :: positiveYaxis, negativeYaxis
  logical :: quadrantI, quadrantII, quadrantIII, quadrantIV

  integer :: blockNr
  integer :: blockID
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: k,kmin,kmax
  integer :: maxEdges
  integer :: messageTag
  integer :: nCellsI, nCellsJ, nCellsK
  integer :: nEdgesI, nEdgesJ, nEdgesK

  integer :: locate      (1:3)
  integer :: status      (MPI_STATUS_SIZE)
  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: alpha, beta
  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: bndBoxKLow
  real    :: cellDelta, cellDensity, cellMass, cellMassDensity, cellVolume
  real    :: cmRcyl, cmZ, cmPhi
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf
  real    :: DeltaK, DeltaKHalf, DeltaKSine, DeltaKHalfSine
  real    :: localMsum, localMDsum, localMDXsum, localMDYsum, localMDZsum
  real    :: maxRcyl, maxZ, maxPhi
  real    :: minRcyl, minZ, minPhi
  real    :: phi, phiCosine, phiSine, phiSineSave
  real    :: Rcyl
  real    :: totalMassDensityInv
  real    :: x,y,z

  real    :: delta     (1:MDIM)
  real    :: localData (1:6)
  real    :: totalData (1:5)
  real    :: bndBox    (LOW:HIGH,1:MDIM)

  real, allocatable :: shifts   (:,:)
  real, pointer     :: solnData (:,:,:,:)
!
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
  localMsum   = ZERO
  localMDsum  = ZERO
  localMDXsum = ZERO
  localMDYsum = ZERO
  localMDZsum = ZERO

  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList  (blockNr)

     call Grid_getBlkBoundBox     (blockID, bndBox)
     call Grid_getDeltas          (blockID, delta)
     call Grid_getBlkPtr          (blockID, solnData)
     call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

     imin = blkLimits (LOW, IAXIS)
     jmin = blkLimits (LOW, JAXIS)
     kmin = blkLimits (LOW, KAXIS)  
     imax = blkLimits (HIGH,IAXIS)
     jmax = blkLimits (HIGH,JAXIS)
     kmax = blkLimits (HIGH,KAXIS)

     DeltaI         = delta (IAXIS)
     DeltaJ         = delta (JAXIS)
     DeltaK         = delta (KAXIS)
     DeltaIHalf     = DeltaI * HALF
     DeltaJHalf     = DeltaJ * HALF
     DeltaKHalf     = DeltaK * HALF
     DeltaKSine     = sin (DeltaK)
     DeltaKHalfSine = sin (DeltaKHalf)

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
     bndBoxKLow = bndBox (LOW,KAXIS)
!
!
!          ...The 3D cylindrical case:
!
!
!                               ------
!                             /        \
!                            /     z    \
!                           |\     |    /|
!                           | \    |   / |
!                           |   ------   |
!                           |      | /   |
!                           |      |/phi |         Rcyl --> stored in i-index (FLASH x)
!                           |      ----->|            z --> stored in j-index (FLASH y)
!                           |       Rcyl |          phi --> stored in k-index (FLASH z)
!                           |            |
!                           |            |
!                           |   ------   |
!                           | /        \ |
!                           |/          \|
!                            \          /
!                             \        /
!                               ------
!
!             The (Rcyl,phi) pair needs to be converted to the corresponding cartesian
!             (x,y) pair, because the location of the center of multipole expansion
!             will be determined in 3D cartesian (x,y,z) coordinates.
!
!             The cell volume is:
!
!                              (1/2) * Dphi * Dz * (R^2 - r^2)
!
!             where r is the left-most (smaller) and R is the right-most (larger)
!             radial cell distance, Dz is the cell's z-axis delta value and Dphi is the
!             cell's angular delta value in radians. Since our radial measure is based
!             on the cell's center, we have: r = Rcyl - Dr/2 and R = Rcyl + Dr/2 with
!             Dr being the cell's radial delta value. Hence the cell volume becomes:
!
!                                  Rcyl * Dr * Dphi * Dz
!
!
     cellDelta = DeltaI * DeltaJ * DeltaK

     alpha     = TWO * DeltaKHalfSine * DeltaKHalfSine
     beta      = DeltaKSine

     phi       = bndBoxKLow + DeltaKHalf
     phiSine   = sin (phi)
     phiCosine = cos (phi)

     do k = kmin, kmax
        z = bndBoxJLow + DeltaJHalf
        do j = jmin, jmax
           Rcyl = bndBoxILow + DeltaIHalf
           do i = imin,imax

              cellVolume      = Rcyl * cellDelta
              cellDensity     = solnData (idensvar,i,j,k)
              cellMass        = cellDensity * cellVolume
              cellMassDensity = cellMass * cellDensity

              x = Rcyl * phiCosine
              y = Rcyl * phiSine

              localMsum   = localMsum   + cellMass
              localMDsum  = localMDsum  + cellMassDensity
              localMDXsum = localMDXsum + cellMassDensity * x
              localMDYsum = localMDYsum + cellMassDensity * y
              localMDZsum = localMDZsum + cellMassDensity * z

              Rcyl = Rcyl + DeltaI
           end do
           z = z + DeltaJ
        end do

        phiSineSave = phiSine
        phiSine     = phiSine   - (alpha * phiSine   - beta * phiCosine  )
        phiCosine   = phiCosine - (alpha * phiCosine + beta * phiSineSave)

     end do

     call Grid_releaseBlkPtr (blockID, solnData)

  end do
!
!
!     ...Prepare for a one-time all reduce call.
!
!
  localData (1) = localMsum
  localData (2) = localMDsum
  localData (3) = localMDXsum
  localData (4) = localMDYsum
  localData (5) = localMDZsum
!
!
!     ...Calculate the total sums and give a copy to each processor.
!
!
  call  MPI_AllReduce (localData,   &
                       totalData,   &
                       5,           &
                       FLASH_REAL,  & 
                       MPI_Sum,     &
                       gr_meshComm, &
                       error        )
!
!
!     ...Analyze total mass obtained. If nonsense, abort.
!
!
  gr_mpoleTotalMass = totalData (1)
     
  if (abs (gr_mpoleTotalMass) < tiny (gr_mpoleTotalMass)) then
      call Driver_abortFlash ('[gr_mpoleCen3Dcylindrical] ERROR:  gr_mpoleTotalMass <= 0')
  end if
!
!
!     ...Calculate center of multipole expansion cartesian coordinates.
!
!
  totalMassDensityInv = ONE / totalData (2)

  gr_mpoleXcenter = totalData (3) * totalMassDensityInv
  gr_mpoleYcenter = totalData (4) * totalMassDensityInv
  gr_mpoleZcenter = totalData (5) * totalMassDensityInv
!
!
!     ...We need to calculate the 3D cylindrical coordinates of the center
!        of multipole expansion in order to determine to which block it belongs.
!        FLASH convention is such that the phi angle of the 3D cylindrical
!        geometry ranges from 0 to 2pi. For accurate angle determination
!        we use the atan function, which returns angles in the range of
!        -pi/2 < phi < pi/2. The following inverse angle calculations are
!        thus performed, depending on which quadrant in the xy-plane the center
!        of multipole expansion lays:
!
!
!                                          y
!                                          |
!                         Quadrant II      |     Quadrant I
!                                          |
!                   phi = pi + atan (y/x)  |  phi = atan (y/x)
!                                          |  
!                                          | /
!                                          |/phi
!                       ----------------------------------------> x
!                                          |
!                                          |
!                         Quadrant III     |     Quadrant I
!                                          |
!                   phi = pi + atan (y/x)  |  phi = 2pi + atan (y/x)
!                                          |
!                                          |
!                                          |
!
!
!        If the center of multipole expansion lays exactly on the y-axis (x = 0),
!        we handle those cases separately.
!
!
  cmZ    = gr_mpoleZcenter
  cmRcyl = sqrt (gr_mpoleXcenter * gr_mpoleXcenter + gr_mpoleYcenter * gr_mpoleYcenter)

  quadrantI     = (gr_mpoleXcenter >  ZERO) .and. (gr_mpoleYcenter >= ZERO)
  quadrantII    = (gr_mpoleXcenter <  ZERO) .and. (gr_mpoleYcenter >= ZERO)
  quadrantIII   = (gr_mpoleXcenter <  ZERO) .and. (gr_mpoleYcenter <  ZERO)
  quadrantIV    = (gr_mpoleXcenter >  ZERO) .and. (gr_mpoleYcenter <  ZERO)
  positiveYaxis = (gr_mpoleXcenter == ZERO) .and. (gr_mpoleYcenter >  ZERO)
  negativeYaxis = (gr_mpoleXcenter == ZERO) .and. (gr_mpoleYcenter <  ZERO)

  if (quadrantI) then
      cmPhi = atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (quadrantII) then
      cmPhi = gr_mpolePi + atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (quadrantIII) then
      cmPhi = gr_mpolePi + atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (quadrantIV) then
      cmPhi = gr_mpoleTwoPi + atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (positiveYaxis) then
      cmPhi = gr_mpoleHalfPi
  else if (negativeYaxis) then
      cmPhi = gr_mpolePi + gr_mpoleHalfPi
  else
      cmPhi = ZERO
  end if
!
!
!     ...Find the local blockID to which the center of multipole expansion
!        belongs and place the center on the nearest cell corner. Also at
!        this point we determine the inner zone atomic length, since the
!        inner zone is defined around the center of multipole expansion.
!        Whatever processor is doing the relevant calculation sends its
!        final data (updated center of multipole expansion and inner zone
!        atomic length) to the master, which then broadcasts the info.
!
!
  messageTag = 1
  invokeRecv = .true.

  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList (blockNr)

     call Grid_getBlkBoundBox (blockID, bndBox)

     minRcyl = bndBox (LOW ,IAXIS)
     maxRcyl = bndBox (HIGH,IAXIS)
     minZ    = bndBox (LOW ,JAXIS)
     maxZ    = bndBox (HIGH,JAXIS)
     minPhi  = bndBox (LOW ,KAXIS)
     maxPhi  = bndBox (HIGH,KAXIS)

     insideBlock  =       (cmRcyl >= minRcyl) &
                    .and. (cmZ    >= minZ   ) &
                    .and. (cmPhi  >= minPhi ) &
                    .and. (cmRcyl <  maxRcyl) &     ! the < instead of <= is
                    .and. (cmZ    <  maxZ   ) &     ! necessary for finding
                    .and. (cmPhi  <  maxPhi )       ! the unique block

     domainRmax   =       (cmRcyl == maxRcyl             ) &    ! include (however unlikely) the
                    .and. (cmRcyl == gr_mpoleDomainRmax  )      ! missing R part of the domain
     domainZmax   =       (cmZ    == maxZ                ) &    ! include (however unlikely) the
                    .and. (cmZ    == gr_mpoleDomainZmax  )      ! missing Z part of the domain
     domainPhiMax =       (cmPhi  == maxPhi              ) &    ! include (however unlikely)) the
                    .and. (cmPhi  == gr_mpoleDomainPhiMax)      ! missing phi part of the domain

     insideBlock  = insideBlock .or. domainRmax .or. domainZmax .or. domainPhiMax

     if (insideBlock) then

         call Grid_getDeltas          (blockID, delta)
         call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

         DeltaI = delta (IAXIS)
         DeltaJ = delta (JAXIS)
         DeltaK = delta (KAXIS)

         gr_mpoleDrInnerZone = HALF * sqrt (DeltaI * DeltaJ)   ! based on Rcyl and Z only

         imin = blkLimits (LOW, IAXIS)
         jmin = blkLimits (LOW, JAXIS)
         kmin = blkLimits (LOW, KAXIS)  
         imax = blkLimits (HIGH,IAXIS)
         jmax = blkLimits (HIGH,JAXIS)
         kmax = blkLimits (HIGH,KAXIS)

         nCellsI = imax - imin + 1
         nCellsJ = jmax - jmin + 1
         nCellsK = kmax - kmin + 1

         nEdgesI = nCellsI + 1
         nEdgesJ = nCellsJ + 1
         nEdgesK = nCellsK + 1

         maxEdges = max (nEdgesI, nEdgesJ, nEdgesK)

         allocate (shifts (1:maxEdges,3))

         guardCells = .false.

         call Grid_getCellCoords (IAXIS, blockID, FACES, guardCells, shifts (1:nEdgesI,1), nEdgesI)
         call Grid_getCellCoords (JAXIS, blockID, FACES, guardCells, shifts (1:nEdgesJ,2), nEdgesJ)
         call Grid_getCellCoords (KAXIS, blockID, FACES, guardCells, shifts (1:nEdgesK,3), nEdgesK)

         shifts (1:nEdgesI,1) = shifts (1:nEdgesI,1) - cmRcyl
         shifts (1:nEdgesJ,2) = shifts (1:nEdgesJ,2) - cmZ
         shifts (1:nEdgesK,3) = shifts (1:nEdgesK,3) - cmPhi

         locate (1) = minloc (abs (shifts (1:nEdgesI,1)), dim = 1)
         locate (2) = minloc (abs (shifts (1:nEdgesJ,2)), dim = 1)
         locate (3) = minloc (abs (shifts (1:nEdgesK,3)), dim = 1)

         cmRcyl = cmRcyl + shifts (locate (1),1)  ! move multipole center to nearest R edge
         cmZ    = cmZ    + shifts (locate (2),2)  ! move multipole center to nearest Z edge
         cmPhi  = cmPhi  + shifts (locate (3),3)  ! move multipole center to nearest Phi edge

         deallocate (shifts)

         gr_mpoleXcenter   = cmRcyl * cos (cmPhi)    !
         gr_mpoleYcenter   = cmRcyl * sin (cmPhi)    ! convert back to 3D cartesian
         gr_mpoleZcenter   = cmZ                     !
         gr_mpoleRcenter   = cmRcyl
         gr_mpolePhiCenter = cmPhi

         localData (1) = gr_mpoleDrInnerZone
         localData (2) = gr_mpoleXcenter
         localData (3) = gr_mpoleYcenter
         localData (4) = gr_mpoleZcenter
         localData (5) = gr_mpoleRcenter
         localData (6) = gr_mpolePhiCenter

         if (gr_meshMe /= MASTER_PE) then

             call MPI_Send (localData,    &
                            6,            &
                            FLASH_REAL,   &
                            MASTER_PE,    &
                            messageTag,   &
                            gr_meshComm,  &
                            error         )
         else
             invokeRecv = .false.
         end if

         exit

     end if
  end do

  if ((gr_meshMe == MASTER_PE) .and. invokeRecv) then

       call MPI_Recv (localData,      &
                      6,              &
                      FLASH_REAL,     &
                      MPI_ANY_SOURCE, &
                      messageTag,     &
                      gr_meshComm,    &
                      status,         &
                      error           )
  end if
!
!
!     ...At this point, the master has all the info. Broadcast and update all
!        other processors.
!
!
  call MPI_Bcast (localData,   &
                  6,           &
                  FLASH_REAL,  &
                  MASTER_PE,   &
                  gr_meshComm, &
                  error        )

  gr_mpoleDrInnerZone    = localData (1)
  gr_mpoleDrInnerZoneInv = ONE / gr_mpoleDrInnerZone
  gr_mpoleXcenter        = localData (2)
  gr_mpoleYcenter        = localData (3)
  gr_mpoleZcenter        = localData (4)
  gr_mpoleRcenter        = localData (5)
  gr_mpolePhiCenter      = localData (6)
!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleCen3Dcylindrical
