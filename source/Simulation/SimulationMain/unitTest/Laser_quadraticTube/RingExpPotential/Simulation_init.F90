!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data, for the laser quadratic
!!  tube unit test.
!!
!!***

subroutine Simulation_init()

  use Simulation_data

  use Grid_interface,              ONLY : Grid_getGeometry,    &
                                          Grid_getMinCellSizes

  use EnergyDeposition_data,       ONLY : ed_Boltzmann,      &
                                          ed_electronCharge, &
                                          ed_electronMass,   &
                                          ed_speedOfLight

  use Driver_interface,            ONLY : Driver_abortFlash,  &
                                          Driver_getComm,     &
                                          Driver_getMype,     &
                                          Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use ed_interface,                ONLY : ed_inverseBremsstrahlungRate
  use ed_extractBeamsData,         ONLY : ed_extractBeamData
  use ed_overrideBeamsData,        ONLY : ed_overrideBeamData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: geometry
  integer :: lrefineMin, lrefineMax

  real    :: ratio
  real    :: nc,e,c,kB,lnF,Tfocal,Iaa

  real    :: minCellSizes (1:MDIM)
!
!
!     ...Get the needed data.
!
!
  call Driver_getComm        (GLOBAL_COMM,                  sim_globalComm         )
  call Driver_getMype        (GLOBAL_COMM,                  sim_globalMe           )
  call Driver_getNumProcs    (GLOBAL_COMM,                  sim_globalNumProcs     )

  call RuntimeParameters_get ("lrefine_min",                lrefineMin             )
  call RuntimeParameters_get ("lrefine_max",                lrefineMax             )
  call RuntimeParameters_get ("basenm",                     sim_baseName           )
  call RuntimeParameters_get ("sim_printBlockVariables",    sim_printBlockVariables)
  call RuntimeParameters_get ("sim_totalNumberOfBoxes",     sim_totalNumberOfBoxes )
  call RuntimeParameters_get ("sim_alpha",                  sim_alpha              )
  call RuntimeParameters_get ("sim_exitPowerFraction",      sim_exitPowerFraction  )
  call RuntimeParameters_get ("ymax",                       sim_yfocal             )

  if (lrefineMin /= lrefineMax) then
      call Driver_abortFlash ('[Simulation_init] ERROR: lrefine_min must match lrefine_max')
  end if

  if (lrefineMax > sim_refinementLevelMax) then
      call Driver_abortFlash ('[Simulation_init] ERROR: lrefine_max too large')
  end if

  sim_refinementLevel = lrefineMax

  call Grid_getGeometry (geometry)

  if ((NDIM /= 3) .or. (geometry /= CARTESIAN)  ) then
      call Driver_abortFlash ('Laser quadratic tube ring problem: geometry must be 3D cartesian!')
  endif

  if (sim_alpha < 1 .or. sim_alpha > 20) then
      call Driver_abortFlash ('[Simulation_init] ERROR: potential exponent alpha out of range!')
  end if

  if (sim_yfocal /= 2.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: The focal point yF must be = 2.0!')
  end if

  if (sim_exitPowerFraction >= 1.0 .or. sim_exitPowerFraction <= 0.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: exit power fraction out of range!')
  end if
!
!
!       ...Set the Boltzman constant locally.
!
!
  sim_Boltzmann = ed_Boltzmann
!
!
!       ...Get the cell sizes for the entire domain. Since the maximum and
!          minimum refinement level are the same, there is only one cell size.
!          The cell sizes for the x- and z-direction must be equal.
!
!
  call Grid_getMinCellSizes (minCellSizes)

  sim_cellSizeX = minCellSizes (IAXIS)
  sim_cellSizeY = minCellSizes (JAXIS)
  sim_cellSizeZ = minCellSizes (KAXIS)

  if (sim_cellSizeX /= sim_cellSizeZ) then
      call Driver_abortFlash ('[Simulation_init] ERROR: cell sizes X and Z not equal!')
  end if
!
!
!     ...Allocate all necessary box arrays.
!
!
  allocate (sim_powersBox (1:sim_totalNumberOfBoxes))
  allocate (sim_radialBox (1:sim_totalNumberOfBoxes))
  allocate (sim_xvalueBox (1:sim_totalNumberOfBoxes))
  allocate (sim_boxPowers (1:sim_totalNumberOfBoxes))
  allocate (sim_boxRadius (1:sim_totalNumberOfBoxes))
  allocate (sim_boxXvalue (1:sim_totalNumberOfBoxes))
!
!
!     ...Retrieve the following data of the beam:
!
!           1) the number of rays that will be launched for the beam
!           2) the beam's radius at the target
!
!
  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'numberOfRays',        &
                             dataValue  = sim_numRaysLaunched    )

  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'targetSemiAxisMajor', &
                             dataValue  = sim_r0                 )

  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'targetX',             &
                             dataValue  = sim_xTubeCenter        )

  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'targetZ',             &
                             dataValue  = sim_zTubeCenter        )

  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'frequency',           &
                             dataValue  = sim_beamFrequency      )

  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'wavelength',          &
                             dataValue  = sim_beamWavelength     )

  if (sim_r0 /= 1.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: The laser ring radius r0 must be = 1.0!')
  end if
!
!
!     ...Calculate the critical electron number density and set the electron number density
!        and the electron temperature (in Kelvin) at the center of the tube.
!
!
  ratio   = sim_beamFrequency / ed_electronCharge
  sim_nc  = ed_electronMass * PI * ratio * ratio

  nc      = sim_nc
  e       = ed_electronCharge
  c       = ed_speedOfLight
  kB      = ed_Boltzmann
  lnF     = log (sim_exitPowerFraction)
  Tfocal  = sim_FocalTime   (sim_alpha)
  Iaa     = sim_IaaIntegral (sim_alpha)

  sim_Te0 = - sqrt ((PI + PI)/ed_electronMass)
  sim_Te0 = sim_Te0 * (nc * e * e * e * e) / (6 * c * lnF * (kB ** 1.5))
  sim_Te0 = sim_Te0 * (Tfocal + Iaa)
  sim_Te0 = sim_Te0 ** (2.0 / 3.0)
!
!
!     ...Calculate the initial ray velocity (in units of c) and store it into the beam
!        (before launching).
!
!
  sim_vy = 1.0 / Tfocal

  call ed_overrideBeamdata    (1,"initialRaySpeed",  sim_vy)

  return
end subroutine Simulation_init
