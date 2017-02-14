!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCen2Dspherical
!!
!! NAME
!!
!!  gr_mpoleCen2Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen2Dspherical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 2D spherical
!!  geometry. It computes the center of expansion for the multipoles for 2D spherical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                          integral (r * rho * rho  dr)
!!                Cen (z) = ----------------------------
!!                            integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                   sum cells (cell center r * cell mass * cell rho)
!!         Cen (z) = ------------------------------------------------
!!                           sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (z) has been determined, it is placed on the nearest
!!  nearest cell corner. The following is computed here:
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
!!  gr_mpoleZcenter denotes the location of the center of multipole expansion
!!  in the 2D cartesian version of the 2D spherical framework (R,theta -> x,z).
!!  Due to symmetry, there is no need to evaluate and store a x-coordinate
!!  component, which is always equal to zero.
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen2Dspherical (idensvar)

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
                                gr_mpoleDomainThetaMax,  &
                                gr_mpoleDrInnerZone,     &
                                gr_mpoleDrInnerZoneInv,  &
                                gr_mpolePi,              &
                                gr_mpoleThirdPi,         &
                                gr_mpoleFourPi,          & 
                                gr_mpoleZcenter,         &
                                gr_mpoleRcenter,         &
                                gr_mpoleThetaCenter,     &
                                gr_mpoleTotalMass,       &
                                gr_mpoleBlockCount,      &
                                gr_mpoleBlockList

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  include "Flash_mpi.h"
  
  integer, intent (in) :: idensvar

  logical :: domainRmax, domainThetaMax
  logical :: guardCells
  logical :: insideBlock
  logical :: invokeRecv

  integer :: blockNr
  integer :: blockID
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: messageTag
  integer :: nCellsI
  integer :: nEdgesI

  integer :: locate      (1:1)
  integer :: status      (MPI_STATUS_SIZE)
  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: alpha, beta
  real    :: angularVolumePart
  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: cellDensity, cellMass, cellMassDensity, cellVolume
  real    :: cmRsph, cmTheta
  real    :: DeltaI, DeltaIHalf, DeltaIcube
  real    :: DeltaJ, DeltaJHalf, DeltaJSine, DeltaJHalfSine
  real    :: localMsum, localMDsum, localMDZsum
  real    :: maxRsph, maxTheta
  real    :: minRsph, minTheta
  real    :: Rsph
  real    :: theta, thetaCosine, thetaSine, thetaSineSave
  real    :: z

  real    :: delta     (1:MDIM)
  real    :: localData (1:3)
  real    :: totalData (1:3)
  real    :: bndBox    (LOW:HIGH,1:MDIM)

  real, allocatable :: shifts   (:)
  real, pointer     :: solnData (:,:,:,:)
!
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
  localMsum   = ZERO
  localMDsum  = ZERO
  localMDZsum = ZERO

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

     DeltaI         = delta (IAXIS)
     DeltaIHalf     = DeltaI * HALF
     DeltaIcube     = DeltaI * DeltaI * DeltaI
     DeltaJ         = delta (JAXIS)
     DeltaJHalf     = DeltaJ * HALF
     DeltaJSine     = sin (DeltaJ)
     DeltaJHalfSine = sin (DeltaJHalf)

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
!
!
!          ...The 2D spherical case:
!
!
!
!                         |*
!                         |     *
!                         |theta/  *
!                         |    /     *
!                         |   /       *
!                         |  /         *              Rsph --> stored in i-index
!                         | / Rsph      *            theta --> stored in j-index
!                         |/            *
!                         |             *
!                         |            *
!                         |           *
!                         |          *
!                         |        *
!                         |     *
!                         |*
!
!
!             Here the convention used in FLASH is to store the Rsph values into the
!             i-index (FLASH x-axis) and the angular theta values (in radians) into
!             the j-index (FLASH y-axis). Due to symmetry, the center of mass only
!             has a radial component and must have theta = 0.
!
!             The cell volume is:
!
!                      (2/3) * pi * (R^3 - r^3) * (cos[t] - cos[T])
!
!             where
!                       r  =  left-most (smaller) radial cell distance
!                       R  =  right-most (larger) radial cell distance
!                       t  =  smaller cell angle
!                       T  =  larger cell angle
!
!             The above formula can be obtained from the simple formula for a spherical
!             cone and noting that the cell can be obtained by adding and subtracting
!             the appropriate spherical cones. The cell looks like a torus with a 2D polar
!             cell cross-sectional shape. Since our angular and radial measures are based
!             on the 2D polar cell's center, we have:
!
!                       r  =  Rsph - Dr/2
!                       R  =  Rsph + Dr/2
!                       t  =  theta - Dt/2
!                       T  =  theta + Dt/2
!
!             with Dt and Dr being the cell's angular and radial delta value.
!             Hence the cell volume becomes, after using some trig identities:
!
!               (4pi * (Rsph)^2 * Dr + (pi/3) * Dr^3) * sin[Dt/2] * sin[theta]
!
!             Note, that the 'cartesian' center of multipole expansion coordinates
!             are obtained. This 'cartesian' framework is solely used to facilitate
!             computations in the spherical framework, especially if the center of
!             multipole expansion shifts away from the spherical domain origin.
!             The two frameworks are related via:
!
!                           x = Rsph * sin (theta)
!                           z = Rsph * cos (theta)
!
!             with the y-cartesian direction being covered by the implicitly assumed
!             symmetry around the z-axis. The center of multipole expansion must be
!             along the z-axis (i.e. have theta = 0 for z >= 0 or theta = pi for z < 0)
!             due to the assumed symmetry around the z-axis.
!
!             The recursions for evaluating the sine and cosine functions for angles
!             which increase by a fixed amount are taken from the Numerical Recipies.
!
!
     alpha       = TWO * DeltaJHalfSine * DeltaJHalfSine
     beta        = DeltaJSine
     theta       = bndBoxJLow + DeltaJHalf
     thetaSine   = sin (theta)
     thetaCosine = cos (theta)

     do j = jmin,jmax

        angularVolumePart = DeltaJHalfSine * thetaSine

        Rsph = bndBoxILow + DeltaIHalf
        do i = imin,imax

           cellVolume      = angularVolumePart * (  gr_mpoleFourPi * Rsph * Rsph * DeltaI  &
                                                  + gr_mpoleThirdPi * DeltaIcube )
           cellDensity     = solnData (idensvar,i,j,1)
           cellMass        = cellDensity * cellVolume
           cellMassDensity = cellMass * cellDensity

           z = Rsph * thetaCosine

           localMsum   = localMsum   + cellMass
           localMDsum  = localMDsum  + cellMassDensity
           localMDZsum = localMDZsum + cellMassDensity * z

           Rsph = Rsph + DeltaI
        end do

        thetaSineSave = thetaSine
        thetaSine     = thetaSine   - (alpha * thetaSine   - beta * thetaCosine  )
        thetaCosine   = thetaCosine - (alpha * thetaCosine + beta * thetaSineSave)
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
  localData (3) = localMDZsum
!
!
!     ...Calculate the total sums and give a copy to each processor.
!
!
  call  MPI_AllReduce (localData,   &
                       totalData,   &
                       3,           &
                       FLASH_REAL,  & 
                       MPI_Sum,     &
                       gr_meshComm, &
                       error        )
!
!
!     ...Proceed according to symmetry situation. Take the total mass x 2 if
!        a symmetry plane has been specified along the theta = 90 degree (pi/2)
!        2D spherical boundary. Analyze total mass obtained. If nonsense, abort.
!        Calculate center of multipole expansion in cartesian z-coordinates. Due
!        to implicit rotational apherical symmetry around the z-axis, the center
!        of multipole expansion must have an x-coordinate of zero and needs not
!        be explicitly stored. For the z-coordinate, override the center of
!        multipole expansion if a symmetry plane is present. In this case the
!        z-coordinate is equal to zero.
!
!
  if (gr_mpoleSymmetryPlane2D) then
      gr_mpoleTotalMass = TWO * totalData (1)
      gr_mpoleZcenter   = ZERO
  else
      gr_mpoleTotalMass = totalData (1)
      gr_mpoleZcenter   = totalData (3) / totalData (2)
  end if

  if (abs (gr_mpoleTotalMass) < tiny (gr_mpoleTotalMass)) then
      call Driver_abortFlash ('[gr_mpoleCen2Dspherical] ERROR:  gr_mpoleTotalMass <= 0')
  end if
!
!
!     ...We need to calculate the 2D spherical coordinates of the center
!        of multipole expansion in order to determine to which block it belongs.
!        Only the radial part is of relevance here.
!
!
  cmRsph = abs (gr_mpoleZcenter)

  if (gr_mpoleZcenter >= ZERO) then
      cmTheta = ZERO
  else
      cmTheta = gr_mpolePi
  end if

  gr_mpoleThetaCenter = cmTheta
!
!
!     ...Find the local blockID to which the center of multipole expansion
!        belongs and place the center on the nearest cell corner. Also at
!        this point we determine the inner zone atomic length, since the
!        inner zone is defined around the center of multipole expansion.
!        Whatever processor is doing the relevant calculation sends its
!        final data (updated center of multipole expansion and inner zone
!        atomic length) to the master, which then broadcasts the info.
!        The inner zone atomic length is based entirely on the radial
!        component.
!
!
  messageTag = 1
  invokeRecv = .true.

  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList (blockNr)

     call Grid_getBlkBoundBox (blockID, bndBox)

     minRsph  = bndBox (LOW ,IAXIS)
     maxRsph  = bndBox (HIGH,IAXIS)
     minTheta = bndBox (LOW ,JAXIS)    ! needed to find the relevant unique block
     maxTheta = bndBox (HIGH,JAXIS)    ! needed to find the relevant unique block

     insideBlock =       (cmRsph  >= minRsph ) &
                   .and. (cmTheta >= minTheta) &
                   .and. (cmRsph  <  maxRsph ) &     ! the < instead of <= is necessary
                   .and. (cmTheta <  maxTheta)       ! for finding the unique block

     domainRmax     =       (cmRsph  == maxRsph               ) &  ! include (however unlikely) the
                      .and. (cmRsph  == gr_mpoleDomainRmax    )    ! missing R part of the domain
     domainThetaMax =       (cmTheta == maxTheta              ) &  ! include (very likely) the
                      .and. (cmTheta == gr_mpoleDomainThetaMax)    ! missing theta part of the domain

     insideBlock = insideBlock .or. domainRmax .or. domainThetaMax

     if (insideBlock) then

         call Grid_getDeltas          (blockID, delta)
         call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

         DeltaI = delta (IAXIS)
         DeltaJ = delta (JAXIS)

         gr_mpoleDrInnerZone = HALF * DeltaI   ! based on Rsph only

         imin = blkLimits (LOW, IAXIS)
         imax = blkLimits (HIGH,IAXIS)

         nCellsI = imax - imin + 1
         nEdgesI = nCellsI + 1

         allocate (shifts (1:nEdgesI))

         guardCells = .false.

         call Grid_getCellCoords (IAXIS, blockID, FACES, guardCells, shifts (1:nEdgesI), nEdgesI)

         shifts (1:nEdgesI) = shifts (1:nEdgesI) - cmRsph
         locate (1) = minloc (abs (shifts (1:nEdgesI)), dim = 1)

         cmRsph = cmRsph + shifts (locate (1))  ! move multipole center to nearest R edge

         deallocate (shifts)

         if (cmTheta == ZERO) then
             gr_mpoleZcenter = cmRsph
         else
             gr_mpoleZcenter = - cmRsph
         end if

         localData (1) = gr_mpoleDrInnerZone
         localData (2) = gr_mpoleZcenter

         if (gr_meshMe /= MASTER_PE) then

             call MPI_Send (localData,    &
                            2,            &
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
                      2,              &
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
                  2,           &
                  FLASH_REAL,  &
                  MASTER_PE,   &
                  gr_meshComm, &
                  error        )

  gr_mpoleDrInnerZone    = localData (1)
  gr_mpoleDrInnerZoneInv = ONE / gr_mpoleDrInnerZone
  gr_mpoleZcenter        = localData (2)
  gr_mpoleRcenter        = abs (gr_mpoleZcenter)
!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleCen2Dspherical
