!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonImaging_init
!!
!! NAME
!!  
!!  ProtonImaging_init
!!
!! SYNOPSIS
!! 
!!  call ProtonImaging_init ()
!!
!! DESCRIPTION
!!
!!  Perform various initializations for the ProtonImaging unit.
!!
!! ARGUMENTS
!!
!!***

subroutine ProtonImaging_init ()

  use ProtonImaging_data

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Logfile_interface,           ONLY : Logfile_stampMessage

  use Grid_interface,              ONLY : Grid_getMinCellSizes, &
                                          Grid_getGeometry

  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getComm,    &
                                          Driver_getMype,    &
                                          Driver_getNumProcs

  use pi_interface,                ONLY : pi_printBeamsData,          &
                                          pi_IOinit,                  &
                                          pi_printDetectorsData,      &
                                          pi_printMainData,           &
                                          pi_setupBeams,              &
                                          pi_setupDetectors,          &
                                          pi_setupDiskProtons,        &
                                          pi_setupProtons,            &
                                          pi_setupScreenProtons,      &
                                          pi_statisticalInit,         &
                                          pi_statisticalSetSeed

#include "ProtonImaging.h"
#include "Flash.h"
#include "constants.h"

#ifdef FLASH_GRID_PARTICLES
#include "GridParticles.h"
#endif

  implicit none

  character (len = MAX_STRING_LENGTH) :: xyFormat, dgFormat

  logical :: grid1D
  logical :: grid2D
  logical :: grid3D
  logical :: gridCartesian
  logical :: gridCylindrical
  logical :: gridSpherical
  logical :: gridPolar

  integer :: rpGeometry

  real    :: domainSizeX
  real    :: domainSizeY
  real    :: domainSizeZ
  real    :: shortestEdge

  real    :: minCellSizes (1:MDIM)
!
!
!     ...Check, if proton imaging is needed at all. If not, return at once.
!
!
  call RuntimeParameters_get ("useProtonImaging",     pi_useProtonImaging)

  if(.not. pi_useProtonImaging) return
!
!
!     ...Get the needed external data.
!
!
  call Driver_getMype        (GLOBAL_COMM,                      pi_globalMe                  )
  call Driver_getComm        (GLOBAL_COMM,                      pi_globalComm                )
  call Driver_getNumProcs    (GLOBAL_COMM,                      pi_globalNumProcs            )
  call Driver_getMype        (  MESH_COMM,                      pi_meshMe                    )
  call Driver_getComm        (  MESH_COMM,                      pi_meshComm                  )
  call Driver_getNumProcs    (  MESH_COMM,                      pi_meshNumProcs              )

  call RuntimeParameters_get ("basenm",                         pi_baseName                  )
  call RuntimeParameters_get ("pi_3Din2D",                      pi_3Din2D                    )
  call RuntimeParameters_get ("pi_3Din2DwedgeAngle",            pi_3Din2DwedgeAngle          )
  call RuntimeParameters_get ("pi_cellStepTolerance",           pi_cellStepTolerance         )
  call RuntimeParameters_get ("pi_cellWallThicknessFactor",     pi_cellWallThicknessFactor   )
  call RuntimeParameters_get ("pi_detectorDGwriteFormat",       pi_detectorDGwriteFormat     )
  call RuntimeParameters_get ("pi_detectorFileNameTimeStamp",   pi_detectorFileNameTimeStamp )
  call RuntimeParameters_get ("pi_detectorXYwriteFormat",       pi_detectorXYwriteFormat     )
  call RuntimeParameters_get ("pi_ignoreElectricalField",       pi_ignoreElectricalField     )
  call RuntimeParameters_get ("pi_IOaddDetectorScreens",        pi_IOaddDetectorScreens      )
  call RuntimeParameters_get ("pi_IOaddProtonsCapsule2Domain",  pi_IOaddProtonsCapsule2Domain)
  call RuntimeParameters_get ("pi_IOaddProtonsDomain2Screen",   pi_IOaddProtonsDomain2Screen )
  call RuntimeParameters_get ("pi_IOmaxBlockCrossingNumber",    pi_IOmaxBlockCrossingNumber  )
  call RuntimeParameters_get ("pi_IOnumberOfProtons2Plot",      pi_IOnumberOfProtons2Plot    )
  call RuntimeParameters_get ("pi_flagDomainMissingProtons",    pi_flagDomainMissingProtons  )
  call RuntimeParameters_get ("pi_maxProtonCount",              pi_maxProtonCount            )
  call RuntimeParameters_get ("pi_numberOfBeams",               pi_numberOfBeams             )
  call RuntimeParameters_get ("pi_numberOfDetectors",           pi_numberOfDetectors         )
  call RuntimeParameters_get ("pi_opaqueBoundaries",            pi_opaqueBoundaries          )
  call RuntimeParameters_get ("pi_printBeams",                  pi_printBeams                )
  call RuntimeParameters_get ("pi_printDetectors",              pi_printDetectors            )
  call RuntimeParameters_get ("pi_printMain",                   pi_printMain                 )
  call RuntimeParameters_get ("pi_printProtons",                pi_printProtons              )
  call RuntimeParameters_get ("pi_protonDeterminism",           pi_protonDeterminism         )
  call RuntimeParameters_get ("pi_randomNumberSeedIncrement",   pi_randomNumberSeedIncrement )
  call RuntimeParameters_get ("pi_randomNumberSeedInitial",     pi_randomNumberSeedInitial   )
  call RuntimeParameters_get ("pi_recalculateCellData",         pi_recalculateCellData       )
  call RuntimeParameters_get ("pi_recordOffScreenProtons",      pi_recordOffScreenProtons    )
  call RuntimeParameters_get ("pi_RungeKuttaMethod",            pi_RungeKuttaMethod          )
  call RuntimeParameters_get ("pi_screenProtonBucketSize",      pi_screenProtonBucketSize    )
  call RuntimeParameters_get ("pi_screenProtonDiagnostics",     pi_screenProtonDiagnostics   )
  call RuntimeParameters_get ("pi_timeResolvedProtonImaging",   pi_timeResolvedProtonImaging )
  call RuntimeParameters_get ("pi_useIOprotonPlot",             pi_useIOprotonPlot           )
  call RuntimeParameters_get ("pi_useParabolicApproximation",   pi_useParabolicApproximation )
  call RuntimeParameters_get ("threadProtonTrace",              pi_threadProtonTrace         )
  call RuntimeParameters_get ("xmin",                           pi_xminDomain                )
  call RuntimeParameters_get ("xmax",                           pi_xmaxDomain                )
  call RuntimeParameters_get ("ymin",                           pi_yminDomain                )
  call RuntimeParameters_get ("ymax",                           pi_ymaxDomain                )
  call RuntimeParameters_get ("zmin",                           pi_zminDomain                )
  call RuntimeParameters_get ("zmax",                           pi_zmaxDomain                )

  call PhysicalConstants_get ("Avogadro",                       pi_Avogadro                  )
  call PhysicalConstants_get ("Boltzmann",                      pi_Boltzmann                 )
  call PhysicalConstants_get ("speed of light",                 pi_speedOfLight              )
  call PhysicalConstants_get ("electron charge",                pi_protonCharge              )
  call PhysicalConstants_get ("proton mass",                    pi_protonMass                )
!
!
!     ...Check for problems.
!
!
  if (pi_numberOfBeams < 1) then
      call Driver_abortFlash ("ProtonImaging_init: No proton imaging beam defined!")
  end if

  if (pi_numberOfBeams > PI_MAXBEAMS) then
      call Logfile_stampMessage ("ProtonImaging_init: ERROR")
      call Logfile_stampMessage ("# of beams > maximum # of beam runtime parameters!")
      call Logfile_stampMessage ("Not enough beam runtime parameters created.")
      call Logfile_stampMessage ("Increase the maximum number of beam runtime parameters")
      call Logfile_stampMessage ("using the pi_maxBeams=<number of beams> setup option.")
      call Driver_abortFlash    ("ProtonImaging_init: Not enough beam runtime parameters. See Log File!")
  end if

  if (pi_numberOfDetectors < 1) then
      call Driver_abortFlash ("ProtonImaging_init: No proton detector defined!")
  end if

  if (pi_numberOfDetectors > PI_MAXDETECTORS) then
      call Logfile_stampMessage ("ProtonImaging_init: ERROR")
      call Logfile_stampMessage ("# of detectors > maximum # of detector runtime parameters!")
      call Logfile_stampMessage ("Not enough detector runtime parameters created.")
      call Logfile_stampMessage ("Increase the maximum number of detector runtime parameters")
      call Logfile_stampMessage ("using the pi_maxDetectors=<number of detectors> setup option.")
      call Driver_abortFlash    ("ProtonImaging_init: Not enough detector runtime parameters. See Log File!")
  end if
!
!
!     ...Set some needed data. Specify the infinite time and the thickness of
!        the cell walls. The cell wall thickness is dependent of the length of the
!        shortest possible cell edge. The infinite time and velocity is set as
!        the largest possible floating point number representable on the machine.
!
!        Set also the computational domain. This is the domain obtained from the user
!        defined domain + the margin for rounding errors. Hence the computational
!        domain is slightly larger than the user defined domain. The use of the
!        computational domain is necessary for judging if protons hit the domain or not.
!        The computational domain is defined as follows, depending on the size
!        of each domain dimension (only X coordinate is shown, the others are the same):
!
!             domainErrorMarginX = (pi_xmaxDomain - pi_xminDomain) * pi_domainTolerance
!
!             computational domain xmin  =  pi_xminDomain - domainErrorMarginX
!             computational domain xmax  =  pi_xminDomain + domainErrorMarginX
!
!
  call Grid_getMinCellSizes (minCellSizes)

  domainSizeX  = pi_xmaxDomain - pi_xminDomain
  domainSizeY  = pi_ymaxDomain - pi_yminDomain
  domainSizeZ  = pi_zmaxDomain - pi_zminDomain
  shortestEdge = minval (minCellSizes (1:NDIM))

  pi_baseName               = adjustl (pi_baseName)                           ! to get ready to use 'trim'
  pi_cellWallThickness      = shortestEdge * pi_cellWallThicknessFactor       ! in cm
  pi_domainErrorMarginX     = domainSizeX  * pi_domainTolerance               ! in cm
  pi_domainErrorMarginY     = domainSizeY  * pi_domainTolerance               ! in cm
  pi_domainErrorMarginZ     = domainSizeZ  * pi_domainTolerance               ! in cm
  pi_infiniteTime           = huge (1.0)                                      ! in s
  pi_infiniteSpeed          = huge (1.0)                                      ! in cm/s
  pi_largestPositiveInteger = huge (1)                                        ! no units
  pi_largestPositiveReal    = huge (1.0)                                      ! no units
  pi_notSetInteger          = - huge (1)                                      ! no units
  pi_notSetReal             = - huge (1.0)                                    ! no units
  pi_protonChargePerMass    = pi_protonCharge / pi_protonMass                 ! in esu/g
  pi_speedOfLightInv        = 1.0 / pi_speedOfLight                           ! in s/cm
  pi_speedOfLightSquared    = pi_speedOfLight * pi_speedOfLight               ! in cm^2/s^2
  pi_squareRoot4Pi          = sqrt (4.0 * PI)                                 ! no units
  pi_unitRoundoff           = epsilon (1.0)                                   ! no units

  if (pi_3Din2D) then
      pi_3Din2DwedgeAngle  = pi_3Din2DwedgeAngle * pi_degrees2rad        ! convert to radians
      pi_3Din2DwedgeCosine = cos (pi_3Din2DwedgeAngle)                   ! no units
      pi_3Din2DwedgeSine   = sin (pi_3Din2DwedgeAngle)                   ! no units
      pi_3Din2DwedgeSlope  = tan (pi_3Din2DwedgeAngle * 0.5)             ! no units
  end if
!
!
!    ...Create a handle to the current geometry.
!
!
  grid3D = (NDIM == 3)
  grid2D = (NDIM == 2)
  grid1D = (NDIM == 1)

  call Grid_getGeometry (rpGeometry)

  gridCartesian   = (rpGeometry == CARTESIAN)
  gridCylindrical = (rpGeometry == CYLINDRICAL)
  gridSpherical   = (rpGeometry == SPHERICAL)
  gridPolar       = (rpGeometry == POLAR)

  if (grid3D .and. gridCartesian)   pi_gridGeometry = GRID_3DCARTESIAN
  if (grid2D .and. gridCartesian)   pi_gridGeometry = GRID_2DCARTESIAN
  if (grid1D .and. gridCartesian)   pi_gridGeometry = GRID_1DCARTESIAN
  if (grid3D .and. gridCylindrical) pi_gridGeometry = GRID_3DCYLINDRICAL
  if (grid2D .and. gridCylindrical) pi_gridGeometry = GRID_2DCYLINDRICAL
  if (grid1D .and. gridCylindrical) pi_gridGeometry = GRID_1DCYLINDRICAL
  if (grid3D .and. gridSpherical)   pi_gridGeometry = GRID_3DSPHERICAL
  if (grid2D .and. gridSpherical)   pi_gridGeometry = GRID_2DSPHERICAL
  if (grid1D .and. gridSpherical)   pi_gridGeometry = GRID_1DSPHERICAL
  if (grid3D .and. gridPolar)       pi_gridGeometry = GRID_3DPOLAR
  if (grid2D .and. gridPolar)       pi_gridGeometry = GRID_2DPOLAR
  if (grid1D .and. gridPolar)       pi_gridGeometry = GRID_1DPOLAR
!
!
!    ...Before proceeding, catch the unsupported geometries (not needed or not
!       yet implemented). Only 3D-type geometries are supported.
!
!
  if (pi_gridGeometry == GRID_2DCARTESIAN   .or. &
      pi_gridGeometry == GRID_1DCARTESIAN   .or. &
      pi_gridGeometry == GRID_3DCYLINDRICAL .or. &
      pi_gridGeometry == GRID_1DCYLINDRICAL .or. &
      pi_gridGeometry == GRID_3DSPHERICAL   .or. &
      pi_gridGeometry == GRID_2DSPHERICAL   .or. &
      pi_gridGeometry == GRID_1DSPHERICAL   .or. &
      pi_gridGeometry == GRID_3DPOLAR       .or. &
      pi_gridGeometry == GRID_2DPOLAR       .or. &
      pi_gridGeometry == GRID_1DPOLAR) then

      call Driver_abortFlash ('[ProtonImaging_init] ERROR: unsupported geometry')
  end if
!
!
!    ...Catch an inconsistent geometry with respect to 3D in 2D proton tracing.
!       Is the wedge opening angle set?
!
!
  if (pi_gridGeometry == GRID_2DCYLINDRICAL .and. .not.pi_3Din2D) then
      call Driver_abortFlash ('[ProtonImaging_init] ERROR: Pure 2D cylindrical grid geometry not allowed')
  end if

  if (pi_3Din2D) then

      if (pi_gridGeometry /= GRID_2DCYLINDRICAL) then
          call Driver_abortFlash ('[ProtonImaging_init] ERROR: Bad 3D in 2D grid geometry')
      end if

      if (pi_3Din2DwedgeAngle <= 0.0 .or. pi_3Din2DwedgeAngle >= 180.0) then
          call Driver_abortFlash ('[ProtonImaging_init] ERROR: Bad wedge angle for 3D in 2D tracing')
      end if

  end if
!
!
!     ...The following data is needed for moving the protons from block to block using
!        the grid particles unit. Note the switch in the particles index list in case of
!        3D in 2D proton tracing, which mapps the proton's z-coordinate to the 2D cylindrical
!        grid y-coordinate.
!
!
#ifdef FLASH_GRID_PARTICLES
  pi_particleIndexCount = GRPT_ALL

  pi_particleIndexList (1:GRPT_ALL) = GRPT_RESET

  if (pi_3Din2D) then
      pi_particleIndexList (GRPT_POSX_IND) = PROTON_POSX
      pi_particleIndexList (GRPT_POSY_IND) = PROTON_POSZ
      pi_particleIndexList (GRPT_POSZ_IND) = PROTON_POSY  ! necessary -> gives (unused) index within allowed range
  else
      pi_particleIndexList (GRPT_POSX_IND) = PROTON_POSX
      pi_particleIndexList (GRPT_POSY_IND) = PROTON_POSY
      pi_particleIndexList (GRPT_POSZ_IND) = PROTON_POSZ
  end if

  pi_particleIndexList (GRPT_BLK_IND ) = PROTON_BLCK
  pi_particleIndexList (GRPT_PROC_IND) = PROTON_PROC
  pi_particleIndexList (GRPT_TAG_IND ) = PROTON_TAGS
#endif
!
!
!     ...Initialize the random number generator.
!
!
  call pi_statisticalInit    ()
  call pi_statisticalSetSeed (pi_randomNumberSeedInitial)
!
!
!     ...Allocate the arrays that will describe proton positions in each beam traget and capsule.
!
!
  allocate (pi_xCircle (1:BEAM_GRIDARRAYSIZE))
  allocate (pi_yCircle (1:BEAM_GRIDARRAYSIZE))
  allocate (pi_xSphere (1:BEAM_GRIDARRAYSIZE))
  allocate (pi_ySphere (1:BEAM_GRIDARRAYSIZE))
  allocate (pi_zSphere (1:BEAM_GRIDARRAYSIZE))
!
!
!     ...Set up the beams, detectors and the proton arrays from input file.
!        Setup of detectors must be after the beams, since detectors might
!        be aligned along specific beams. Set up disk protons only if time
!        resolved proton imaging is going to be performed.
!
!
  call pi_setupBeams         ()
  call pi_setupDetectors     ()
  call pi_setupProtons       ()
  call pi_setupScreenProtons ()

  if (pi_timeResolvedProtonImaging) then
      call pi_setupDiskProtons   ()
  end if
!
!
!    ...Initialize IO for plotting protons. The IO must be initialized after! the proton
!       beams, because we need to know exactly the total number of protons that will be
!       emitted over all beams.
!
!
  call pi_IOinit ()
!
!
!    ...Set the file name for monitoring the proton imaging progress.
!
!
  pi_monitorFileName = trim (pi_baseName) // "ProtonImagingMonitor.dat"
!
!
!     ...Prepare the format for writing lines of data to the detector(s). Allocate the
!        needed arrays on the master processor to be used for flushing screen protons
!        to disk (allocate single sized arrays on the other nodes to avoid crashes in
!        debug mode).
!
!
  write (xyFormat,'(i1,a)') 2, pi_detectorXYwriteFormat
  write (dgFormat,'(i1,a)') pi_numberOfDiagnostics, pi_detectorDGwriteFormat

  if (pi_screenProtonDiagnostics) then
      pi_detectorLNwriteFormat   = "(" // trim (xyFormat) // "," // trim (dgFormat) // ")"
      pi_detectorLNwriteFormat   = trim (pi_detectorLNwriteFormat)
      pi_screenProtonRecordCount = 2 + pi_numberOfDiagnostics         ! currently: x,y,Jv,Kx,Ky,Kz
  else
      pi_detectorLNwriteFormat   = "(" // trim (xyFormat) // ")"
      pi_detectorLNwriteFormat   = trim (pi_detectorLNwriteFormat)
      pi_screenProtonRecordCount = 2                                  ! just: x,y
  end if

  if (pi_globalMe == MASTER_PE) then
      allocate (pi_screenProtonCountOffsets (0:pi_globalNumProcs-1))
      allocate (pi_screenProtonCountProcs   (0:pi_globalNumProcs-1))
      allocate (pi_screenProtonBucketCount  (1:pi_numberOfDetectors))
      allocate (pi_screenProtonBuckets      (1:pi_screenProtonRecordCount, &
                                             1:pi_screenProtonBucketSize,  &
                                             1:pi_numberOfDetectors        ))
  else
      allocate (pi_screenProtonCountOffsets (0:0))
      allocate (pi_screenProtonCountProcs   (0:0))
      allocate (pi_screenProtonBucketCount  (1:1))
      allocate (pi_screenProtonBuckets      (1:1,1:1,1:1))
  end if
!
!
!     ...If time resolved proton imaging is going to be performed, do preparations
!        for accumulating and writing disk protons to disk.
!
!
  if (pi_timeResolvedProtonImaging) then
      if (pi_globalMe == MASTER_PE) then
          allocate (pi_diskProtonCountOffsets (0:pi_globalNumProcs-1))
          allocate (pi_diskProtonCountProcs   (0:pi_globalNumProcs-1))
      else
          allocate (pi_diskProtonCountOffsets (0:0))
          allocate (pi_diskProtonCountProcs   (0:0))
      end if
  end if
!
!
!     ...If requested, print some data.
!
!  
  if (pi_printMain)      call pi_printMainData      ()
  if (pi_printBeams)     call pi_printBeamsData     ()
  if (pi_printDetectors) call pi_printDetectorsData ()
!
!
!    ...Ready!
!
!
  return
end subroutine ProtonImaging_init
