!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingCubicPath/Simulation_init
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

  use Grid_interface,              ONLY : Grid_getGeometry

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

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: geometry
  integer :: lrefineMin, lrefineMax

  real    :: integral
  real    :: lnLambda
  real    :: ratio
  real    :: tcross
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
  call RuntimeParameters_get ("xmax",                       sim_xc                 )
  call RuntimeParameters_get ("ymax",                       sim_yfocal             )
  call RuntimeParameters_get ("zmax",                       sim_zc                 )

  call ed_extractBeamdata    (1,"frequency",                sim_beamFrequency      )
  call ed_extractBeamdata    (1,"wavelength",               sim_beamWavelength     )

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
!
!
!     ...Allocate all necessary box arrays.
!
!
  allocate (sim_radialBox (1:sim_totalNumberOfBoxes))
  allocate (sim_xvalueBox (1:sim_totalNumberOfBoxes))
  allocate (sim_boxRadius (1:sim_totalNumberOfBoxes))
  allocate (sim_boxXvalue (1:sim_totalNumberOfBoxes))
!
!
!     ...Set the tube center position and the radius around the center on which
!        the rays will be launched.
!
!
  sim_xw =  5.0
  sim_zw =  5.0
  sim_R  =  3.0
!
!
!     ...Calculate the critical electron number density and set the electron number density
!        and the electron temperature (in Kelvin) at the center of the tube. Set the average
!        Z number throughout the tube as 1 and define the quadratic potential factor P.
!
!
  ratio      = sim_beamFrequency / ed_electronCharge
  sim_nc     = ed_electronMass * PI * ratio * ratio
  sim_nw     = sim_nc
  sim_Tw     = 10.0 * sim_keV2Kelvin

  return
end subroutine Simulation_init
