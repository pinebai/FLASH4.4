!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/EnergyDeposition_init
!!
!! NAME
!!  
!!  EnergyDeposition_init
!!
!! SYNOPSIS
!! 
!!  call EnergyDeposition_init ()
!!
!! DESCRIPTION
!!
!!  Perform various initializations for the EnergyDeposition unit.
!!
!! ARGUMENTS
!!
!!***

subroutine EnergyDeposition_init ()

  use EnergyDeposition_data

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Grid_interface,              ONLY : Grid_getMinCellSizes, &
                                          Grid_getGeometry

  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getComm,    &
                                          Driver_getMype,    &
                                          Driver_getNumProcs

  use Logfile_interface,           ONLY : Logfile_stampVarMask

  use Simulation_interface,        ONLY : Simulation_getVarnameType,&
                                          Simulation_mapStrToInt

  use ed_interface,                ONLY : ed_laserIOInit,     &
                                          ed_printBeamsData,  &
                                          ed_printMainData,   &
                                          ed_printPulsesData, &
                                          ed_setupBeams,      &
                                          ed_setupPulses,     &
                                          ed_setupRays

  use ed_commInterface,            ONLY : ed_commInit

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

#ifdef FLASH_GRID_PARTICLES
# include "GridParticles.h"
#endif

  implicit none

  logical :: grid1D
  logical :: grid2D
  logical :: grid3D
  logical :: gridCartesian
  logical :: gridCylindrical
  logical :: gridSpherical
  logical :: gridPolar

  integer :: ut_getFreeFileUnit
  integer :: rpGeometry
  integer :: depoVarType

  real    :: domainSizeX
  real    :: domainSizeY
  real    :: domainSizeZ
  real    :: shortestEdge

  real    :: minCellSizes (1:MDIM)
!
!
!     ...Check, if energy deposition is needed at all. If not, return at once.
!
!
  call RuntimeParameters_get ("useEnergyDeposition",          ed_useEnergyDeposition)

  if(.not. ed_useEnergyDeposition) return
!
!
!     ...Get the needed external data.
!
!
  call Driver_getMype        (GLOBAL_COMM,                      ed_globalMe                  )
  call Driver_getComm        (GLOBAL_COMM,                      ed_globalComm                )
  call Driver_getNumProcs    (GLOBAL_COMM,                      ed_globalNumProcs            )
  call Driver_getMype        (  MESH_COMM,                      ed_meshMe                    )
  call Driver_getComm        (  MESH_COMM,                      ed_meshComm                  )
  call Driver_getNumProcs    (  MESH_COMM,                      ed_meshNumProcs              )

  call RuntimeParameters_get ("basenm",                         ed_baseName                  )
  call RuntimeParameters_get ("ed_cellStepTolerance",           ed_cellStepTolerance         )
  call RuntimeParameters_get ("ed_cellWallThicknessFactor",     ed_cellWallThicknessFactor   )
  call RuntimeParameters_get ("ed_computeGradNeleX",            ed_computeGradNeleX          )
  call RuntimeParameters_get ("ed_computeGradNeleY",            ed_computeGradNeleY          )
  call RuntimeParameters_get ("ed_computeGradNeleZ",            ed_computeGradNeleZ          )
  call RuntimeParameters_get ("ed_cubicInterpolationZeroDerv",  ed_cubicInterpolationZeroDerv)
  call RuntimeParameters_get ("ed_enforcePositiveNele",         ed_enforcePositiveNele       )
  call RuntimeParameters_get ("ed_enforcePositiveTele",         ed_enforcePositiveTele       )
  call RuntimeParameters_get ("ed_gradOrder",                   ed_gradOrder                 )
  call RuntimeParameters_get ("ed_laser3Din2D",                 ed_laser3Din2D               )
  call RuntimeParameters_get ("ed_laser3Din2DwedgeAngle",       ed_laser3Din2DwedgeAngle     )
  call RuntimeParameters_get ("ed_maxRayCount",                 ed_maxRayCount               )
  call RuntimeParameters_get ("ed_numberOfBeams",               ed_numberOfBeams             )
  call RuntimeParameters_get ("ed_numberOfPulses",              ed_numberOfPulses            )
  call RuntimeParameters_get ("ed_powerStepTolerance",          ed_powerStepTolerance        )
  call RuntimeParameters_get ("ed_printBeams",                  ed_printBeams                )
  call RuntimeParameters_get ("ed_printMain",                   ed_printMain                 )
  call RuntimeParameters_get ("ed_printPulses",                 ed_printPulses               )
  call RuntimeParameters_get ("ed_printRays",                   ed_printRays                 )
  call RuntimeParameters_get ("ed_rayDeterminism",              ed_rayDeterminism            )
  call RuntimeParameters_get ("ed_rayZeroPower",                ed_rayZeroPower              )
  call RuntimeParameters_get ("ed_RungeKuttaMethod",            ed_RungeKuttaMethod          )
  call RuntimeParameters_get ("ed_saveOutOfDomainRays",         ed_saveOutOfDomainRays       )
  call RuntimeParameters_get ("threadRayTrace",                 ed_threadRayTrace            )
  call RuntimeParameters_get ("xmin",                           ed_xminDomain                )
  call RuntimeParameters_get ("xmax",                           ed_xmaxDomain                )
  call RuntimeParameters_get ("ymin",                           ed_yminDomain                )
  call RuntimeParameters_get ("ymax",                           ed_ymaxDomain                )
  call RuntimeParameters_get ("zmin",                           ed_zminDomain                )
  call RuntimeParameters_get ("zmax",                           ed_zmaxDomain                )

  call RuntimeParameters_get ("ed_depoReuseMaxSteps",           ed_depoReuseMaxSteps         ) 
  call RuntimeParameters_get ("ed_depoVarName",                 ed_depoVarName               ) 

  call Simulation_mapStrToInt(ed_depoVarName, ed_depoVar, MAPBLOCK_UNK)

  if (ed_depoVar < 1 .OR. ed_depoVar > NUNK_VARS) then
     if (ed_globalMe==MASTER_PE) then
        print*,'[EnergyDeposition_init] ERROR: ed_depoVar=',ed_depoVar,'ed_depoVarName=',ed_depoVarName
     end if
     call Driver_abortFlash("[EnergyDeposition_init] ERROR: ed_depoVarName is invalid.")
  end if
  call Simulation_getVarnameType(ed_depoVar, depoVarType)
  if (depoVarType == VARTYPE_ERROR) then
     call Driver_abortFlash("[EnergyDeposition_init] ERROR: faliled to determine variable type of depoVar.")
  end if
  ed_depoVarIsPerMass = (depoVarType .NE. VARTYPE_PER_VOLUME)


  call RuntimeParameters_get ("ed_irradVarName",                ed_irradVarName              ) 
  call Simulation_mapStrToInt(ed_irradVarName, ed_irradVar, MAPBLOCK_UNK)

  if (ed_irradVar < 1 .OR. ed_irradVar > NUNK_VARS) then
     ed_irradVar =  -1
  else
     if (ed_globalMe==MASTER_PE) then
        print*,'[EnergyDeposition_init] INFO: Using ed_irradVar=',ed_irradVar,'ed_irradVarName=',ed_irradVarName
     end if
  end if


  !! Allow selective guardcell fill calls ---------------------------------------
  ed_gcMaskSize = NUNK_VARS
  allocate(ed_gcMask(ed_gcMaskSize))
  ed_gcMask = .FALSE.
#ifdef TELE_VAR
  ed_gcMask(TELE_VAR) = .TRUE.
#endif
#ifdef NELE_VAR
  ed_gcMask(NELE_VAR) = .TRUE.
#endif
#ifdef SUMY_MSCALAR
  ed_gcMask(SUMY_MSCALAR) = .TRUE.
#endif
#ifdef YE_MSCALAR
  ed_gcMask(YE_MSCALAR) = .TRUE.
#endif

#if NSPECIES > 0
#ifdef SPECIES_BEGIN
  ed_gcMask(SPECIES_BEGIN:SPECIES_END) = .TRUE.
#endif
#endif

  call Logfile_stampVarMask(ed_gcMask, .FALSE., '[EnergyDeposition_init]', 'gcNeed')





  call PhysicalConstants_get ("Avogadro",                     ed_Avogadro               )
  call PhysicalConstants_get ("Boltzmann",                    ed_Boltzmann              )
  call PhysicalConstants_get ("speed of light",               ed_speedOfLight           )
  call PhysicalConstants_get ("electron charge",              ed_electronCharge         )
  call PhysicalConstants_get ("electron mass",                ed_electronMass           ) 
      
!
!
!     ...Set some needed data. Specify the infinite time and the thickness of
!        the cell walls. The cell wall thickness is dependent of the length of the
!        shortest possible cell edge. The infinite time, velocity and power is set as
!        the largest possible floating point number representable on the machine.
!
!        Set also the computational domain. This is the domain obtained from the user
!        defined domain + the margin for rounding errors. Hence the computational
!        domain is slightly larger than the user defined domain. The use of the
!        computational domain is necessary for judging if rays hit the domain or not.
!        The computational domain is defined as follows, depending on the size
!        of each domain dimension (only X coordinate is shown, the others are the same):
!
!             domainErrorMarginX = (ed_xmaxDomain - ed_xminDomain) * ed_domainTolerance
!
!             computational domain xmin  =  ed_xminDomain - domainErrorMarginX
!             computational domain xmax  =  ed_xminDomain + domainErrorMarginX
!
!
  call Grid_getMinCellSizes (minCellSizes)

  domainSizeX  = ed_xmaxDomain - ed_xminDomain
  domainSizeY  = ed_ymaxDomain - ed_yminDomain
  domainSizeZ  = ed_zmaxDomain - ed_zminDomain
  shortestEdge = minval (minCellSizes (1:NDIM))

  ed_baseName               = adjustl (ed_baseName)                           ! to get ready to use 'trim'
  ed_cellWallThickness      = shortestEdge * ed_cellWallThicknessFactor       ! in cm
  ed_domainErrorMarginX     = domainSizeX  * ed_domainTolerance               ! in cm
  ed_domainErrorMarginY     = domainSizeY  * ed_domainTolerance               ! in cm
  ed_domainErrorMarginZ     = domainSizeZ  * ed_domainTolerance               ! in cm
  ed_energyInTotal          = 0.0                                             ! in erg
  ed_energyOutTotal         = 0.0                                             ! in erg
  ed_infinitePower          = huge (1.0)                                      ! in erg/s
  ed_infiniteTime           = huge (1.0)                                      ! in s
  ed_infiniteSpeed          = huge (1.0)                                      ! in cm/s
  ed_largestPositiveInteger = huge (1)                                        ! no units
  ed_largestPositiveReal    = huge (1.0)                                      ! no units
  ed_notSetInteger          = - huge (1)                                      ! no units
  ed_notSetReal             = - huge (1.0)                                    ! no units
  ed_speedOfLightSquared    = ed_speedOfLight * ed_speedOfLight               ! in cm^2/s^2
  ed_unitRoundoff           = epsilon (1.0)                                   ! no units

  if (ed_laser3Din2D) then
      ed_laser3Din2DwedgeAngle  = ed_laser3Din2DwedgeAngle * ed_degrees2rad   ! convert to radians
      ed_laser3Din2DwedgeCosine = cos (ed_laser3Din2DwedgeAngle)              ! no units
      ed_laser3Din2DwedgeSine   = sin (ed_laser3Din2DwedgeAngle)              ! no units
      ed_laser3Din2DwedgeSlope  = tan (ed_laser3Din2DwedgeAngle * 0.5)        ! no units
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

  if (grid3D .and. gridCartesian)   ed_gridGeometry = GRID_3DCARTESIAN
  if (grid2D .and. gridCartesian)   ed_gridGeometry = GRID_2DCARTESIAN
  if (grid1D .and. gridCartesian)   ed_gridGeometry = GRID_1DCARTESIAN
  if (grid3D .and. gridCylindrical) ed_gridGeometry = GRID_3DCYLINDRICAL
  if (grid2D .and. gridCylindrical) ed_gridGeometry = GRID_2DCYLINDRICAL
  if (grid1D .and. gridCylindrical) ed_gridGeometry = GRID_1DCYLINDRICAL
  if (grid3D .and. gridSpherical)   ed_gridGeometry = GRID_3DSPHERICAL
  if (grid2D .and. gridSpherical)   ed_gridGeometry = GRID_2DSPHERICAL
  if (grid1D .and. gridSpherical)   ed_gridGeometry = GRID_1DSPHERICAL
  if (grid3D .and. gridPolar)       ed_gridGeometry = GRID_3DPOLAR
  if (grid2D .and. gridPolar)       ed_gridGeometry = GRID_2DPOLAR
  if (grid1D .and. gridPolar)       ed_gridGeometry = GRID_1DPOLAR
!
!
!    ...Before proceeding, catch the unsupported geometries (not needed or not
!       yet implemented).
!
!
  if (ed_gridGeometry == GRID_3DCYLINDRICAL .or. &
      ed_gridGeometry == GRID_1DCYLINDRICAL .or. &
      ed_gridGeometry == GRID_3DSPHERICAL   .or. &
      ed_gridGeometry == GRID_2DSPHERICAL   .or. &
      ed_gridGeometry == GRID_1DSPHERICAL   .or. &
      ed_gridGeometry == GRID_3DPOLAR       .or. &
      ed_gridGeometry == GRID_2DPOLAR       .or. &
      ed_gridGeometry == GRID_1DPOLAR) then

      call Driver_abortFlash ('[EnergyDeposition_init] ERROR: unsupported geometry')
  end if
!
!
!    ...Catch an inconsistent geometry with respect to 3D in 2D ray tracing.
!       Is the wedge opening angle set?
!
!
  if (ed_laser3Din2D) then

      if (ed_gridGeometry /= GRID_2DCYLINDRICAL) then
          call Driver_abortFlash ('[EnergyDeposition_init] ERROR: Bad 3D in 2D grid geometry')
      end if

      if (ed_laser3Din2DwedgeAngle <= 0.0 .or. ed_laser3Din2DwedgeAngle >= 180.0) then
          call Driver_abortFlash ('[EnergyDeposition_init] ERROR: Bad wedge angle for 3D in 2D ray trace')
      end if

  end if
!
!
!     ...The following data is needed for moving the rays from block to block using
!        the grid particles unit. The rays are essentially treated as particles.
!        Note the switch in the particles index list in case of laser 3D in 2D ray tracing,
!        which mapps the ray's z-coordinate to the 2D cylindrical grid y-coordinate.
!
!
#ifdef FLASH_GRID_PARTICLES
  ed_particleIndexCount = GRPT_ALL

  ed_particleIndexList (1:GRPT_ALL) = GRPT_RESET

  if (ed_laser3Din2D) then
      ed_particleIndexList (GRPT_POSX_IND) = RAY_POSX
      ed_particleIndexList (GRPT_POSY_IND) = RAY_POSZ
      ed_particleIndexList (GRPT_POSZ_IND) = RAY_POSY  ! necessary -> gives (unused) index within allowed range
  else
      ed_particleIndexList (GRPT_POSX_IND) = RAY_POSX
      ed_particleIndexList (GRPT_POSY_IND) = RAY_POSY
      ed_particleIndexList (GRPT_POSZ_IND) = RAY_POSZ
  end if

  ed_particleIndexList (GRPT_BLK_IND ) = RAY_BLCK
  ed_particleIndexList (GRPT_PROC_IND) = RAY_PROC
  ed_particleIndexList (GRPT_TAG_IND ) = RAY_TAGS
#endif
!
!
!     ...Set up the pulses and beams from input file. Also setup the ray array.
!
!
  call ed_setupPulses ()
  call ed_setupBeams  ()
  call ed_setupRays   ()
!
!
!    ...Initialize laser IO. Laser IO must be initialized after the beams, because
!       we need to know the total number of rays over all beams.
!
!
  call ed_laserIOInit ()
!
!
!    ...Open the file that will contain the energy stamp printout for each timestep.
!       This is done only on the master processor.
!
!
  ed_energyProfileFileUnit = ut_getFreeFileUnit ()
  ed_energyProfileFileName = trim (ed_baseName) // "LaserEnergyProfile.dat"

  if (ed_globalMe == MASTER_PE) then
      open (unit = ed_energyProfileFileUnit, file = ed_energyProfileFileName)
  end if
!
!
!     ...If requested, print some data.
!
!  
  if (ed_printMain)   call ed_printMainData   ()
  if (ed_printPulses) call ed_printPulsesData ()
  if (ed_printBeams)  call ed_printBeamsData  ()
!
!
!    ...Ready!
!
!
  call ed_CommInit()
  return
end subroutine EnergyDeposition_init
