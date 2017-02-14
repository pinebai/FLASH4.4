!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCen2Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleCen2Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen2Dcylindrical (integer, intent (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 2D cylindrical
!!  geometry. It computes the center of expansion for the multipoles for 2D cylindrical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                          integral (r * rho * rho  dr)
!!              Cen (R,z) = ----------------------------
!!                            integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                   sum cells (cell center r * cell mass * cell rho)
!!       Cen (R,z) = ------------------------------------------------
!!                           sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (R,z) has been determined, it is placed on the
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
!!  gr_mpoleRcenter and gr_mpoleZcenter denote the location of the center
!!  of multipole expansion in the 2D cylindrical framework.
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen2Dcylindrical (idensvar)

  use Grid_data,         ONLY : gr_meshMe,   &
                                gr_meshComm

  use Driver_interface,  ONLY : Driver_abortFlash

  use Grid_interface,    ONLY : Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkBoundBox,    &
                                Grid_getDeltas,         &
                                Grid_getBlkIndexLimits, &
                                Grid_getCellCoords

  use gr_mpoleData,      ONLY : gr_mpoleSymmetryPlane2D, &
                                gr_mpoleDomainRmax,      &
                                gr_mpoleDomainZmax,      &
                                gr_mpoleDrInnerZone,     &
                                gr_mpoleDrInnerZoneInv,  &
                                gr_mpoleTwoPi,           & 
                                gr_mpoleRcenter,         &
                                gr_mpoleZcenter,         &
                                gr_mpoleTotalMass,       &
                                gr_mpoleBlockCount,      &
                                gr_mpoleBlockList,       &
                                gr_mpoleXdens2CoM,       &
                                gr_mpoleYdens2CoM,       &
                                gr_mpoleXcenterOfMass,   &
                                gr_mpoleYcenterOfMass

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  include "Flash_mpi.h"
  
  integer, intent (in) :: idensvar

  logical :: domainRmax, domainZmax
  logical :: guardCells
  logical :: insideBlock
  logical :: invokeRecv

  integer :: blockNr
  integer :: blockID
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: maxEdges
  integer :: messageTag
  integer :: nCellsI, nCellsJ
  integer :: nEdgesI, nEdgesJ

  integer :: locate      (1:2)
  integer :: status      (MPI_STATUS_SIZE)
  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: cellPlane, cellDensity, cellMass, cellMassDensity, cellVolume
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf
  real    :: localMsum, localMDsum, localMDYsum, localMYsum
  real    :: maxRcyl, maxZ
  real    :: minRcyl, minZ
  real    :: Rcyl
  real    :: z

  real    :: delta     (1:MDIM)
  real    :: localData (1:4)
  real    :: totalData (1:4)
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
  localMDYsum = ZERO
  localMYsum  = ZERO

  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList  (blockNr)

     call Grid_getBlkBoundBox     (blockID, bndBox)
     call Grid_getDeltas          (blockID, delta)
     call Grid_getBlkPtr          (blockID, solnData)
     call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

     imin = blkLimits (LOW, IAXIS)
     jmin = blkLimits (LOW, JAXIS)
     imax = blkLimits (HIGH,IAXIS)
     jmax = blkLimits (HIGH,JAXIS)

     DeltaI = delta (IAXIS)
     DeltaJ = delta (JAXIS)

     DeltaIHalf = DeltaI * HALF
     DeltaJHalf = DeltaJ * HALF

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
!
!
!          ...The 2D cylindrical case:
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
!             Due to the implied rotation around the z-axis, the center of multipole
!             expansion must have a  FLASH x-coordinate of zero. The only nonzero
!             coordinate of the center of multipole expansion must be along the FLASH
!             y-coordinate.
!
!             The cell volume is:
!
!                               pi * (R^2 - r^2) * Dz
!
!             where r is the left-most (smaller) and R is the right-most (larger)
!             radial cell distance and Dz is the cell's z-axis delta value. Since
!             our radial measure is based on the cell's center, we have: r = Rcyl - Dr/2
!             and R = Rcyl + Dr/2 with Dr being the cell's radial delta value.
!             Hence the cell volume becomes:
!
!                             2 * pi * Rcyl * Dr * Dz
!
!
     cellPlane = DeltaI * DeltaJ

     z = bndBoxJLow + DeltaJHalf
     do j = jmin,jmax
        Rcyl = bndBoxILow + DeltaIHalf
        do i = imin,imax

           cellVolume      = gr_mpoleTwoPi * Rcyl * cellPlane
           cellDensity     = solnData (idensvar,i,j,1)
           cellMass        = cellDensity * cellVolume
           cellMassDensity = cellMass * cellDensity

           localMsum   = localMsum   + cellMass
           localMDsum  = localMDsum  + cellMassDensity
           localMDYsum = localMDYsum + cellMassDensity * z
           localMYsum  = localMYsum  + cellMass * z

           Rcyl = Rcyl + DeltaI
        end do
        z = z + DeltaJ
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
  localData (3) = localMDYsum
  localData (4) = localMYsum
!
!
!     ...Calculate the total sums and give a copy to each processor.
!
!
  call  MPI_AllReduce (localData,   &
                       totalData,   &
                       4,           &
                       FLASH_REAL,  & 
                       MPI_Sum,     &
                       gr_meshComm, &
                       error        )
!
!
!     ...Proceed according to symmetry situation. Take the total mass x 2 if
!        a symmetry plane has been specified along the radial (FLASH x-axis)
!        2D cylindrical boundary. Analyze total mass obtained. If nonsense,
!        abort. Calculate center of multipole expansion in FLASH coordinates.
!        Due to implicit rotational cylindrical symmetry the center of multipole
!        expansion must have a R-coordinate of zero. For the z-coordinate,
!        override the center of multipole expansion if a symmetry plane is
!        present. In this case the z-coordinate is also equal to zero.
!
!
  if (gr_mpoleSymmetryPlane2D) then
      gr_mpoleTotalMass = TWO * totalData (1)
      gr_mpoleRcenter   = ZERO
      gr_mpoleZcenter   = ZERO
      gr_mpoleXcenterOfMass = ZERO
      gr_mpoleYcenterOfMass = ZERO
      gr_mpoleXdens2CoM = ZERO
      gr_mpoleYdens2CoM = ZERO
  else
      gr_mpoleTotalMass = totalData (1)
      gr_mpoleRcenter   = ZERO
      gr_mpoleZcenter   = totalData (3) / totalData (2)
      gr_mpoleXcenterOfMass = ZERO
      gr_mpoleYcenterOfMass = totalData(4) / totalData(1)
      gr_mpoleXdens2CoM = ZERO
      gr_mpoleYdens2CoM = gr_mpoleZcenter
  end if

  if (abs (gr_mpoleTotalMass) < tiny (gr_mpoleTotalMass)) then
      call Driver_abortFlash ('[gr_mpoleCen2Dcylindrical] ERROR:  gr_mpoleTotalMass <= 0')
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

     insideBlock =       (gr_mpoleRcenter >= minRcyl) &
                   .and. (gr_mpoleZcenter >= minZ   ) &
                   .and. (gr_mpoleRcenter <  maxRcyl) &               ! the < instead of <= is necessary
                   .and. (gr_mpoleZcenter <  maxZ   )                 ! for finding the unique block

     domainRmax  =       (gr_mpoleRcenter == maxRcyl           ) &    ! include (however unlikely) the
                   .and. (gr_mpoleRcenter == gr_mpoleDomainRmax)      ! missing R part of the domain
     domainZmax  =       (gr_mpoleZcenter == maxZ              ) &    ! include (however unlikely) the
                   .and. (gr_mpoleZcenter == gr_mpoleDomainZmax)      ! missing Z part of the domain

     insideBlock = insideBlock .or. domainRmax .or. domainZmax

     if (insideBlock) then

         call Grid_getDeltas          (blockID, delta)
         call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

         DeltaI = delta (IAXIS)
         DeltaJ = delta (JAXIS)

         gr_mpoleDrInnerZone = HALF * sqrt (DeltaI * DeltaJ)

         imin = blkLimits (LOW, IAXIS)
         jmin = blkLimits (LOW, JAXIS)
         imax = blkLimits (HIGH,IAXIS)
         jmax = blkLimits (HIGH,JAXIS)

         nCellsI = imax - imin + 1
         nCellsJ = jmax - jmin + 1

         nEdgesI = nCellsI + 1
         nEdgesJ = nCellsJ + 1

         maxEdges = max (nEdgesI, nEdgesJ)

         allocate (shifts (1:maxEdges,2))

         guardCells = .false.

         call Grid_getCellCoords (IAXIS, blockID, FACES, guardCells, shifts (1:nEdgesI,1), nEdgesI)
         call Grid_getCellCoords (JAXIS, blockID, FACES, guardCells, shifts (1:nEdgesJ,2), nEdgesJ)

         shifts (1:nEdgesI,1) = shifts (1:nEdgesI,1) - gr_mpoleRcenter
         shifts (1:nEdgesJ,2) = shifts (1:nEdgesJ,2) - gr_mpoleZcenter

         locate (1) = minloc (abs (shifts (1:nEdgesI,1)), dim = 1)
         locate (2) = minloc (abs (shifts (1:nEdgesJ,2)), dim = 1)

         gr_mpoleRcenter = gr_mpoleRcenter + shifts (locate (1),1)  ! move multipole center to nearest R edge
         gr_mpoleZcenter = gr_mpoleZcenter + shifts (locate (2),2)  ! move multipole center to nearest Z edge

         deallocate (shifts)

         localData (1) = gr_mpoleDrInnerZone
         localData (2) = gr_mpoleRcenter
         localData (3) = gr_mpoleZcenter

         if (gr_meshMe /= MASTER_PE) then

             call MPI_Send (localData,    &
                            3,            &
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
                      3,              &
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
                  3,           &
                  FLASH_REAL,  &
                  MASTER_PE,   &
                  gr_meshComm, &
                  error        )

  gr_mpoleDrInnerZone    = localData (1)
  gr_mpoleDrInnerZoneInv = ONE / gr_mpoleDrInnerZone
  gr_mpoleRcenter        = localData (2)
  gr_mpoleZcenter        = localData (3)
!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleCen2Dcylindrical
