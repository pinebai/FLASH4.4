!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/pi_createProtons3DRec
!!
!! NAME
!!
!!  pi_createProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_createProtons3DRec (integer, intent (in) :: blockCount,
!!                              integer, intent (in) :: blockList (:),
!!                              real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates protons and places them in their initial blocks for those geometries consisting
!!  formally of 3D rectangular grids (cartesian). On exit, all protons hitting the domain boundary
!!  have been generated for the current processor. Their block ID's are not ordered as the
!!  outer loop is over all beams.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeSimulation : current simulation time
!!
!! NOTES
!!
!!***

subroutine pi_createProtons3DRec (blockCount, blockList, timeSimulation)

  use Driver_interface,    ONLY : Driver_abortFlash
 
  use Grid_interface,      ONLY : Grid_getBlkBC,      &
                                  Grid_getBlkBoundBox

  use pi_interface,        ONLY : pi_statisticalXYcircle,     &
                                  pi_statisticalXYZsphere,    &
                                  pi_capsuleGrainIndices2xyz, &
                                  pi_capsuleNextGrainIndices
 
  use ProtonImaging_data,  ONLY : pi_beams,               &
                                  pi_cellWallThickness,   &
                                  pi_domainErrorMarginX,  &
                                  pi_domainErrorMarginY,  &
                                  pi_domainErrorMarginZ,  &
                                  pi_globalMe,            &
                                  pi_largestPositiveReal, &
                                  pi_maxProtonCount,      &
                                  pi_numberOfBeams,       &
                                  pi_protonCount,         &
                                  pi_protons,             &
                                  pi_xCircle,             &
                                  pi_yCircle,             &
                                  pi_xSphere,             &
                                  pi_ySphere,             &
                                  pi_zSphere,             &
                                  pi_xminDomain,          &
                                  pi_xmaxDomain,          &
                                  pi_yminDomain,          &
                                  pi_ymaxDomain,          &
                                  pi_zminDomain,          &
                                  pi_zmaxDomain
  
  implicit none

#include "ProtonImaging.h"
#include "constants.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeSimulation

  logical :: beamActive
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: blockReflectMinX, blockReflectMaxX
  logical :: blockReflectMinY, blockReflectMaxY
  logical :: blockReflectMinZ, blockReflectMaxZ
  logical :: ignoreBoundary
  logical :: inBlock
  logical :: onDomainFace
  logical :: parallel2XYplane
  logical :: parallel2XZplane
  logical :: parallel2YZplane
  logical :: proceed
  logical :: protonMissesDomain, protonReflects
  logical :: xmaxLimit, ymaxLimit, zmaxLimit
  logical :: xminLimit, yminLimit, zminLimit
  logical :: valid

  integer :: beam
  integer :: blockID
  integer :: block
  integer :: detector
  integer :: grainIndexI, grainIndexJ, grainIndexK
  integer :: grainLevel
  integer :: grainProtons
  integer :: n
  integer :: nCircle, nSphere
  integer :: nProtons, nProtonsBeam, nProtonsCapsule, &
             nProtonsCreate, nProtonsPerGrain, nProtonsProcessed
  integer :: seed

  real    :: capsuleRadius
  real    :: capsuleX, capsuleY, capsuleZ
  real    :: cellWallThicknessHalf
  real    :: compDomainXmin, compDomainXmax
  real    :: compDomainYmin, compDomainYmax
  real    :: compDomainZmin, compDomainZmax
  real    :: Cx, Cy, Cz
  real    :: distToFaceMinX, distToFaceMinY, distToFaceMinZ
  real    :: distToFaceMaxX, distToFaceMaxY, distToFaceMaxZ
  real    :: grainLocalX, grainLocalY, grainLocalZ
  real    :: grainSize
  real    :: initialProtonSpeed
  real    :: Px, Py, Pz
  real    :: posX, posY, posZ
  real    :: RsizeInv
  real    :: Rx, Ry, Rz
  real    :: RxInv, RyInv, RzInv
  real    :: targetRadius
  real    :: targetX, targetY, targetZ
  real    :: time2Launch
  real    :: Tx, Ty, Tz
  real    :: u1x, u1y, u1z
  real    :: u2x, u2y, u2z
  real    :: u3x, u3y, u3z
  real    :: uRx, uRy, uRz
  real    :: userDomainXmin, userDomainXmax
  real    :: userDomainYmin, userDomainYmax
  real    :: userDomainZmin, userDomainZmax
  real    :: velX, velY, velZ
  real    :: w, wtest
  real    :: xminBlock, yminBlock, zminBlock
  real    :: xmaxBlock, ymaxBlock, zmaxBlock

  integer :: faces  (LOW:HIGH,1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  integer, allocatable :: blockBoundary (:,:)
  real,    allocatable :: blockBndBox   (:,:)
!
!
!     ...In order to avoid repeated calls to the grid get bounding box and boundary
!        condition function inside the innermost loop for each proton, we determine the
!        bounding box and boundary consitions for each block beforehand and store
!        this info away for future reference.
!
!
  allocate (blockBndBox   (1:6,1:blockCount))
  allocate (blockBoundary (1:6,1:blockCount))

  do block = 1, blockCount

     blockID = blockList (block)

     call Grid_getBlkBC       (blockID, faces)
     call Grid_getBlkBoundBox (blockID, bndBox)

     blockBndBox   (1,block) = bndBox (LOW ,IAXIS)
     blockBndBox   (2,block) = bndBox (HIGH,IAXIS)
     blockBndBox   (3,block) = bndBox (LOW ,JAXIS)
     blockBndBox   (4,block) = bndBox (HIGH,JAXIS)
     blockBndBox   (5,block) = bndBox (LOW ,KAXIS)
     blockBndBox   (6,block) = bndBox (HIGH,KAXIS)

     blockBoundary (1,block) = faces  (LOW ,IAXIS)
     blockBoundary (2,block) = faces  (HIGH,IAXIS)
     blockBoundary (3,block) = faces  (LOW ,JAXIS)
     blockBoundary (4,block) = faces  (HIGH,JAXIS)
     blockBoundary (5,block) = faces  (LOW ,KAXIS)
     blockBoundary (6,block) = faces  (HIGH,KAXIS)

  end do
!
!
!     ...Set the computational domain boundaries. These will be used to decide, whether
!        protons hit the domain boundaries. The computational domain is slightly larger
!        in each dimensional direction to allow for computational rounding errors. Usage
!        of the user defined domain only results in potential loss of protons due to boundary
!        misses.
!
!
  userDomainXmin = pi_xminDomain
  userDomainXmax = pi_xmaxDomain
  userDomainYmin = pi_yminDomain
  userDomainYmax = pi_ymaxDomain
  userDomainZmin = pi_zminDomain
  userDomainZmax = pi_zmaxDomain

  compDomainXmin = pi_xminDomain - pi_domainErrorMarginX
  compDomainXmax = pi_xmaxDomain + pi_domainErrorMarginX
  compDomainYmin = pi_yminDomain - pi_domainErrorMarginY
  compDomainYmax = pi_ymaxDomain + pi_domainErrorMarginY
  compDomainZmin = pi_zminDomain - pi_domainErrorMarginZ
  compDomainZmax = pi_zmaxDomain + pi_domainErrorMarginZ
!
!
!     ...Set some extra needed data.
!
!
  cellWallThicknessHalf = 0.5 * pi_cellWallThickness
!
!
!     ...Outer loop over all (active) proton beams.
!
!
  nProtonsProcessed = 0

  do beam = 1, pi_numberOfBeams

     time2Launch  = pi_beams (beam) % time2Launch
     nProtonsBeam = pi_beams (beam) % numberOfProtons

     beamActive = (time2Launch <= timeSimulation) .and. (nProtonsBeam > 0)

     if (beamActive) then

         capsuleRadius        = pi_beams (beam) % capsuleRadius
         capsuleX             = pi_beams (beam) % capsuleX
         capsuleY             = pi_beams (beam) % capsuleY
         capsuleZ             = pi_beams (beam) % capsuleZ
         detector             = pi_beams (beam) % detector
         grainIndexI          = pi_beams (beam) % capsuleGrainIndexI
         grainIndexJ          = pi_beams (beam) % capsuleGrainIndexJ
         grainIndexK          = pi_beams (beam) % capsuleGrainIndexK
         grainLevel           = pi_beams (beam) % capsuleGrainLevel
         grainLocalX          = pi_beams (beam) % capsuleGrainLocalX
         grainLocalY          = pi_beams (beam) % capsuleGrainLocalY
         grainLocalZ          = pi_beams (beam) % capsuleGrainLocalZ
         grainProtons         = pi_beams (beam) % capsuleGrainProtons
         grainSize            = pi_beams (beam) % capsuleGrainSize
         ignoreBoundary       = pi_beams (beam) % noBoundaryCondition
         initialProtonSpeed   = pi_beams (beam) % initialProtonSpeed
         nProtonsPerGrain     = pi_beams (beam) % numberOfProtonsPerGrain
         seed                 = pi_beams (beam) % randomNumberSeed
         targetRadius         = pi_beams (beam) % targetRadius
         targetX              = pi_beams (beam) % targetX
         targetY              = pi_beams (beam) % targetY
         targetZ              = pi_beams (beam) % targetZ
         u1x                  = pi_beams (beam) % axisUnit1X
         u1y                  = pi_beams (beam) % axisUnit1Y
         u1z                  = pi_beams (beam) % axisUnit1Z
         u2x                  = pi_beams (beam) % axisUnit2X
         u2y                  = pi_beams (beam) % axisUnit2Y
         u2z                  = pi_beams (beam) % axisUnit2Z
         u3x                  = pi_beams (beam) % axisUnit3X
         u3y                  = pi_beams (beam) % axisUnit3Y
         u3z                  = pi_beams (beam) % axisUnit3Z
!
!
!     ...Decide how many protons of the beam can be safely created to fit into memory.
!        The number of remaining protons in the beam is adjusted for later creation.
!        Once we have determined the number of protons that can be created, we start
!        getting them from the (still active) capsule grain(s). As many active capsule
!        grains are depleted of their protons as can safely be accommodated into the
!        proton array.
!
!
         if (nProtonsProcessed + nProtonsBeam > pi_maxProtonCount) then
             nProtonsCreate    = pi_maxProtonCount - nProtonsProcessed
             nProtonsProcessed = pi_maxProtonCount
             nProtonsBeam      = nProtonsBeam - nProtonsCreate
         else
             nProtonsCreate    = nProtonsBeam
             nProtonsProcessed = nProtonsProcessed + nProtonsCreate
             nProtonsBeam      = 0
         end if

         pi_beams (beam) % numberOfProtons = nProtonsBeam                       ! reset remaining # of protons

         if (nProtonsCreate == 0) then                                          ! exits the beam loop -> no more
             exit                                                               ! protons can be accommodated
         end if                                                                 ! for the beams that follow

         do while (nProtonsCreate > 0)

            if (grainLevel > 0) then                                             ! non-statistical capsule
                if (grainProtons == 0) then                                      ! current capsule grain is empty
                    call pi_capsuleNextGrainIndices (grainLevel,              &  ! get new capsule grain
                                                                 grainIndexI, &  ! new capsule grain index I
                                                                 grainIndexJ, &  ! new capsule grain index J
                                                                 grainIndexK, &  ! new capsule grain index K
                                                                 valid        )  !
                    if (.not.valid) then
                        call Driver_abortFlash ("pi_createProtons3DRec: Beam has protons, but 0 capsule grains!")
                    end if

                    call pi_capsuleGrainIndices2xyz (grainIndexI,             &
                                                     grainIndexJ,             &
                                                     grainIndexK,             &
                                                     grainSize,               &
                                                                 grainLocalX, &  ! capsule grain local X
                                                                 grainLocalY, &  ! capsule grain local Y
                                                                 grainLocalZ  )  ! capsule grain local Z

                    pi_beams (beam) % capsuleGrainIndexI = grainIndexI           !
                    pi_beams (beam) % capsuleGrainIndexJ = grainIndexJ           !
                    pi_beams (beam) % capsuleGrainIndexK = grainIndexK           ! save current capsule data into
                    pi_beams (beam) % capsuleGrainLocalX = grainLocalX           ! beams array for further reference
                    pi_beams (beam) % capsuleGrainLocalY = grainLocalY           !
                    pi_beams (beam) % capsuleGrainLocalZ = grainLocalZ           !

                    grainProtons = nProtonsPerGrain                              ! update grain # of protons
                end if

                nProtonsCapsule = min (grainProtons, nProtonsCreate)    ! # of protons to be created from capsule
                nProtonsCreate  = nProtonsCreate  -  nProtonsCapsule    ! adjust # of total protons to be created
                grainProtons    = grainProtons    -  nProtonsCapsule    ! remaining # of protons in capsule grain

                pi_beams (beam) % capsuleGrainProtons = grainProtons    ! update capsule grain # of protons
            else
                nProtonsCapsule = nProtonsCreate                        ! statistical capsule
                nProtonsCreate  = 0                                     ! create all protons at once
            end if
!
!
!     ...Create all protons for the current capsule grain. This is done by looping statistically
!        over all possible grid points in the circular target area.
!
!        For each proton trajectory we have already a global x,y,z position C for the capsule
!        grain. From the (to be caluclated) global x,y,z position T in the target area, the
!        proton trajectory can be represented by a vector 'R' with tail point 'C' on the capsule
!        grain and head point 'T' at the target. It therefore lays on a line in 3D, given by the
!        following parametric line vector equation:
!
!                                      P = C + w * R
!                                        = C + w * (T - C)
!
!        where 'w' is a real parameter that can attain all possible values. In particular we
!        must have:
!
!                              w >= 0    -->   proton moves from C to T
!
!        Having the proton line equation, we now must see if and where the domain boundaries
!        are crossed. For each domain boundary we have very simple plane equations:
!
!                 domain yz-faces  -->  plane equations: x = xminDomain and x = xmaxDomain
!                 domain xz-faces  -->  plane equations: y = yminDomain and y = ymaxDomain
!                 domain xy-faces  -->  plane equations: z = zminDomain and z = zmaxDomain
!
!        Lets assume we have a plane equation x = a. Then we find the 'w' value, such that Px = a
!        is on the plane. We obtain from the proton line equation:
!
!                                    w = (a - Cx) / (Tx - Cx)
!
!        From this we get the corresponding Py and Pz coordinates:
!
!                                   Py = Cy + w * (Ty - Cy)
!                                   Pz = Cz + w * (Tz - Cz)
!
!        Any w > 0 values are acceptable for a specific plane crossing. w <= 0 must be excluded,
!        since this would mean that the capsule has been placed (partially) inside the domain,
!        a situation which we must exclude. Hence we have:
!
!                     i)     w <= 0  -->  capsule grain (partially) inside the domain (exclude)
!                    ii)     w >  0  -->  capsule grain completely outside the domain (ok)
!
!        If the 'w' value is ok, we have to check, if the point (Py,Pz) is actually contained
!        on the domain yz-face. If yes, add this 'w' to the allowed w-list. To find the
!        relevant 'w' we need to pick the smallest 'w' from the allowed w-list. If the list
!        turns out to be empty, then we know that the proton will not hit the domain.
!
!
!
            do while (nProtonsCapsule > 0)

               if (grainLevel > 0) then
                   call Driver_abortFlash ("pi_createProtons3DRec: No capsule grains > 0 allowed!")
               end if

               call pi_statisticalXYZsphere (capsuleRadius,                  &
                                             BEAM_GRIDARRAYSIZE,             &
                                                                 nSphere,    &
                                                                 pi_xSphere, &   ! capsule grain local X
                                                                 pi_ySphere, &   ! capsule grain local Y
                                                                 pi_zSphere  )   ! capsule grain local Z


               call pi_statisticalXYcircle (targetRadius,                   &
                                            BEAM_GRIDARRAYSIZE,             &
                                                                nCircle,    &
                                                                pi_xCircle, &        ! target local X
                                                                pi_yCircle  )        ! target local Y

               nProtons = min (nProtonsCapsule, nCircle, nSphere)

               nProtonsCapsule = nProtonsCapsule - nProtons

               do n = 1, nProtons    ! Tx,Ty,Tz = target global xyz's , Cx,Cy,Cz = capsule grain global xyz's

                  Tx =  targetX + pi_xCircle (n) * u1x + pi_yCircle (n) * u2x
                  Ty =  targetY + pi_xCircle (n) * u1y + pi_yCircle (n) * u2y
                  Tz =  targetZ + pi_xCircle (n) * u1z + pi_yCircle (n) * u2z

                  Cx = capsuleX + pi_xSphere (n) * u1x + pi_ySphere (n) * u2x + pi_zSphere (n) * u3x
                  Cy = capsuleY + pi_xSphere (n) * u1y + pi_ySphere (n) * u2y + pi_zSphere (n) * u3y
                  Cz = capsuleZ + pi_xSphere (n) * u1z + pi_ySphere (n) * u2z + pi_zSphere (n) * u3z

                  Rx = Tx - Cx           ! proton trajectory vector global X
                  Ry = Ty - Cy           ! proton trajectory vector global Y
                  Rz = Tz - Cz           ! proton trajectory vector global Z

                  parallel2YZplane = (Rx == 0.0)
                  parallel2XZplane = (Ry == 0.0)
                  parallel2XYplane = (Rz == 0.0)
!
!
!     ...Start the search for the smallest 'w' such that w > 0. Check first, if the proton hits
!        any of the two boundary yz-faces of the domain. Note, that while the check is being performed
!        using the computational domain, the protons will be placed exactly on the user defined domain.
!
!
                  protonMissesDomain = .true.

                  w = pi_largestPositiveReal

                  if (.not. parallel2YZplane) then

                       RxInv   = 1.0 / Rx
                       wtest   = (userDomainXmin - Cx) * RxInv
                       proceed = (wtest > 0.0) .and. (wtest <= w)

                       if (proceed) then
                           Py = Cy + wtest * Ry
                           Pz = Cz + wtest * Rz
                           onDomainFace =      (Py >= compDomainYmin) .and. (Py <= compDomainYmax) &
                                         .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                           if (onDomainFace) then
                               w = wtest
                               posX = userDomainXmin
                               posY = min (  Py, userDomainYmax)     !
                               posY = max (posY, userDomainYmin)     ! force protons on user defined domain
                               posZ = min (  Pz, userDomainZmax)     !
                               posZ = max (posZ, userDomainZmin)     !
                               protonMissesDomain = .false.
                           end if
                       end if

                       wtest   = (userDomainXmax - Cx) * RxInv
                       proceed = (wtest > 0.0) .and. (wtest <= w)

                       if (proceed) then
                           Py = Cy + wtest * Ry
                           Pz = Cz + wtest * Rz
                           onDomainFace =      (Py >= compDomainYmin) .and. (Py <= compDomainYmax) &
                                         .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                           if (onDomainFace) then
                               w = wtest
                               posX = userDomainXmax
                               posY = min (  Py, userDomainYmax)
                               posY = max (posY, userDomainYmin)
                               posZ = min (  Pz, userDomainZmax)
                               posZ = max (posZ, userDomainZmin)
                               protonMissesDomain = .false.
                           end if
                       end if

                  end if
!
!
!     ...Next check, if proton hits any of the two boundary xz-faces of the domain.
!
!
                  if (.not. parallel2XZplane) then

                       RyInv   = 1.0 / Ry
                       wtest   = (userDomainYmin - Cy) * RyInv
                       proceed = (wtest > 0.0) .and. (wtest <= w)

                       if (proceed) then
                           Px = Cx + wtest * Rx
                           Pz = Cz + wtest * Rz
                           onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                         .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                           if (onDomainFace) then
                               w = wtest
                               posY = userDomainYmin
                               posX = min (  Px, userDomainXmax)
                               posX = max (posX, userDomainXmin)
                               posZ = min (  Pz, userDomainZmax)
                               posZ = max (posZ, userDomainZmin)
                               protonMissesDomain = .false.
                           end if
                       end if

                       wtest   = (userDomainYmax - Cy) * RyInv
                       proceed = (wtest > 0.0) .and. (wtest <= w)

                       if (proceed) then
                           Px = Cx + wtest * Rx
                           Pz = Cz + wtest * Rz
                           onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                         .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                           if (onDomainFace) then
                               w = wtest
                               posY = userDomainYmax
                               posX = min (  Px, userDomainXmax)
                               posX = max (posX, userDomainXmin)
                               posZ = min (  Pz, userDomainZmax)
                               posZ = max (posZ, userDomainZmin)
                               protonMissesDomain = .false.
                           end if
                       end if

                  end if
!
!
!     ...Finally, check, if proton hits any of the two boundary xy-faces of the domain.
!
!
                  if (.not. parallel2XYplane) then

                       RzInv   = 1.0 / Rz
                       wtest   = (userDomainZmin - Cz) * RzInv
                       proceed = (wtest > 0.0) .and. (wtest <= w)

                       if (proceed) then
                           Px = Cx + wtest * Rx
                           Py = Cy + wtest * Ry
                           onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                         .and. (Py >= compDomainYmin) .and. (Py <= compDomainYmax)
                           if (onDomainFace) then
                               w = wtest
                               posZ = userDomainZmin
                               posX = min (  Px, userDomainXmax)
                               posX = max (posX, userDomainXmin)
                               posY = min (  Py, userDomainYmax)
                               posY = max (posY, userDomainYmin)
                               protonMissesDomain = .false.
                           end if
                       end if

                       wtest   = (userDomainZmax - Cz) * RzInv
                       proceed = (wtest > 0.0) .and. (wtest <= w)

                       if (proceed) then
                           Px = Cx + wtest * Rx
                           Py = Cy + wtest * Ry
                           onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                         .and. (Py >= compDomainYmin) .and. (Py <= compDomainYmax)
                           if (onDomainFace) then
                               w = wtest
                               posZ = userDomainZmax
                               posX = min (  Px, userDomainXmax)
                               posX = max (posX, userDomainXmin)
                               posY = min (  Py, userDomainYmax)
                               posY = max (posY, userDomainYmin)
                               protonMissesDomain = .false.
                           end if
                       end if

                   end if
!
!
!     ...Catch the situation in which the proton does not hit the domain. At this moment,
!        we treat this as an error. For the future we might think about this simply as
!        a proton traveling not through the domain (we might then have still to check, if
!        the proton is incident on the recording screen).
!
!
                  if (protonMissesDomain) then
                      call Driver_abortFlash("pi_createProtons3DRec: proton does not hit domain boundary!")
                  end if
!
!
!     ...loop over all blocks and see, if the proton is contained in one of it. As soon as
!        that block is found, exit the block loop, since each proton can be assigned to only
!        one block. A proton belongs to a block, if its x,y,z coordinate is such that:
!
!                        x,y,z block lower limit  <=  proton x,y,z  <  x,y,z block upper limit
!
!        except when any of the upper limits of the block coincide with the domain boundaries.
!        In that case the less '<' sign on the right must be replaced by a <= sign, otherwise
!        protons will dissapear and not accounted for, resulting in proton loss.
!
!
                  do block = 1, blockCount

                     xminBlock = blockBndBox (1,block)
                     xmaxBlock = blockBndBox (2,block)
                     yminBlock = blockBndBox (3,block)
                     ymaxBlock = blockBndBox (4,block)
                     zminBlock = blockBndBox (5,block)
                     zmaxBlock = blockBndBox (6,block)

                     xminLimit = (posX >= xminBlock)
                     yminLimit = (posY >= yminBlock)
                     zminLimit = (posZ >= zminBlock)
                     xmaxLimit = (posX <  xmaxBlock) .or. ((posX == xmaxBlock) .and. (xmaxBlock == userDomainXmax))
                     ymaxLimit = (posY <  ymaxBlock) .or. ((posY == ymaxBlock) .and. (ymaxBlock == userDomainYmax))
                     zmaxLimit = (posZ <  zmaxBlock) .or. ((posZ == zmaxBlock) .and. (zmaxBlock == userDomainZmax))

                     inBlock  =      xminLimit &
                               .and. xmaxLimit &
                               .and. yminLimit &
                               .and. ymaxLimit &
                               .and. zminLimit &
                               .and. zmaxLimit

                     if (inBlock) then

                         distToFaceMinX = abs (xminBlock - posX)
                         distToFaceMaxX = abs (xmaxBlock - posX)
                         distToFaceMinY = abs (yminBlock - posY)
                         distToFaceMaxY = abs (ymaxBlock - posY)
                         distToFaceMinZ = abs (zminBlock - posZ)
                         distToFaceMaxZ = abs (zmaxBlock - posZ)

                         blockFaceMinX = (distToFaceMinX < cellWallThicknessHalf)
                         blockFaceMaxX = (distToFaceMaxX < cellWallThicknessHalf)
                         blockFaceMinY = (distToFaceMinY < cellWallThicknessHalf)
                         blockFaceMaxY = (distToFaceMaxY < cellWallThicknessHalf)
                         blockFaceMinZ = (distToFaceMinZ < cellWallThicknessHalf)
                         blockFaceMaxZ = (distToFaceMaxZ < cellWallThicknessHalf)

                         blockReflectMinX = (blockBoundary (1,block) == REFLECTING)
                         blockReflectMaxX = (blockBoundary (2,block) == REFLECTING)
                         blockReflectMinY = (blockBoundary (3,block) == REFLECTING)
                         blockReflectMaxY = (blockBoundary (4,block) == REFLECTING)
                         blockReflectMinZ = (blockBoundary (5,block) == REFLECTING)
                         blockReflectMaxZ = (blockBoundary (6,block) == REFLECTING)

                         protonReflects =     (blockFaceMinX .and. blockReflectMinX) &
                                         .or. (blockFaceMaxX .and. blockReflectMaxX) &
                                         .or. (blockFaceMinY .and. blockReflectMinY) &
                                         .or. (blockFaceMaxY .and. blockReflectMaxY) &
                                         .or. (blockFaceMinZ .and. blockReflectMinZ) &
                                         .or. (blockFaceMaxZ .and. blockReflectMaxZ)

                         if (ignoreBoundary .or. .not.protonReflects) then

                             pi_protonCount = pi_protonCount + 1        ! this count is per processor, not block

                             if (pi_protonCount > pi_maxProtonCount) then
                                 call Driver_abortFlash ("pi_createProtons3DRec: No storage left for proton array")
                             end if

                             blockID = blockList (block)

                             RsizeInv = 1.0 / sqrt (Rx * Rx + Ry * Ry + Rz * Rz)

                             uRx = Rx * RsizeInv                        ! unit vector along proton vector
                             uRy = Ry * RsizeInv
                             uRz = Rz * RsizeInv

                             velX = uRx * initialProtonSpeed            !
                             velY = uRy * initialProtonSpeed            ! in cm/s
                             velZ = uRz * initialProtonSpeed            !
!
!
!     ...We now need to make sure, that the proton is placed inside the block and not within the block
!        wall. If this is not done, then initial cell index determination can suffer from roundoff
!        errors with placement of the proton in a wrong cell ouside the current block. Note, that even
!        after pushing the proton inside the block is done, there is still the possibility, that the proton
!        sits exactly on a cell face boundary, but this will happen inside the block, so there is no
!        danger of placing the proton initially in the wrong block. The appropriate (initial) nudging for
!        the cell will be done in the corresponding proton tracing routine, as here we do not have the
!        individual cell information.
!
!
                             if (blockFaceMinX) then
                                 posX = xminBlock + cellWallThicknessHalf
                             else if (blockFaceMaxX) then
                                 posX = xmaxBlock - cellWallThicknessHalf
                             end if

                             if (blockFaceMinY) then
                                 posY = yminBlock + cellWallThicknessHalf
                             else if (blockFaceMaxY) then
                                 posY = ymaxBlock - cellWallThicknessHalf
                             end if

                             if (blockFaceMinZ) then
                                 posZ = zminBlock + cellWallThicknessHalf
                             else if (blockFaceMaxZ) then
                                 posZ = zmaxBlock - cellWallThicknessHalf
                             end if

                             pi_protons (PROTON_POSX,pi_protonCount) = posX
                             pi_protons (PROTON_POSY,pi_protonCount) = posY
                             pi_protons (PROTON_POSZ,pi_protonCount) = posZ
                             pi_protons (PROTON_VELX,pi_protonCount) = velX
                             pi_protons (PROTON_VELY,pi_protonCount) = velY
                             pi_protons (PROTON_VELZ,pi_protonCount) = velZ
                             pi_protons (PROTON_BLCK,pi_protonCount) = real (blockID)
                             pi_protons (PROTON_PROC,pi_protonCount) = real (pi_globalMe)
                             pi_protons (PROTON_BEAM,pi_protonCount) = real (beam)
                             pi_protons (PROTON_DETC,pi_protonCount) = real (detector)
                             pi_protons (PROTON_DGJV,pi_protonCount) = 0.0
                             pi_protons (PROTON_DGKX,pi_protonCount) = 0.0
                             pi_protons (PROTON_DGKY,pi_protonCount) = 0.0
                             pi_protons (PROTON_DGKZ,pi_protonCount) = 0.0

                             exit                                 ! exits the block loop

                         end if                                   ! proton reflection condition
                     end if                                       ! proton in block
                  end do                                          ! block loop
               end do                                             ! individual protons in capsule
            end do                                                ! chunks of capsule protons loop
         end do                                                   ! beam proton creation loop
     end if                                                       ! active beam condition
  end do                                                          ! beam loop
!
!
!     ...Deallocate all allocated arrays.
!
!
  deallocate (blockBndBox)
  deallocate (blockBoundary)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_createProtons3DRec
