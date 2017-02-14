!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/Simulation_init
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

  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getComm,    &
                                          Driver_getMype

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
  call Driver_getMype        (GLOBAL_COMM,                  sim_globalMe           )
  call Driver_getComm        (GLOBAL_COMM,                  sim_globalComm         )

  call RuntimeParameters_get ("lrefine_min",                lrefineMin             )
  call RuntimeParameters_get ("lrefine_max",                lrefineMax             )
  call RuntimeParameters_get ("basenm",                     sim_baseName           )
  call RuntimeParameters_get ("sim_lasersOrientation",      sim_lasersOrientation  )
  call RuntimeParameters_get ("sim_printBlockVariables",    sim_printBlockVariables)
  call RuntimeParameters_get ("xmax",                       sim_domainXmax         )
  call RuntimeParameters_get ("ymax",                       sim_domainYmax         )
  call RuntimeParameters_get ("zmax",                       sim_domainZmax         )

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

  if (     (NDIM == 3) .and. (geometry == CARTESIAN)  ) then
      sim_geometry = GRID_3DCARTESIAN
  else if ((NDIM == 2) .and. (geometry == CARTESIAN)  ) then
      sim_geometry = GRID_2DCARTESIAN
  else if ((NDIM == 2) .and. (geometry == CYLINDRICAL)) then
      sim_geometry = GRID_2DCYLINDRICAL
  else
      call Driver_abortFlash ('Laser quadratic tube unit test I : unsupported geometry!')
  endif
!
!
!     ...Catch geometry / lasers orientation mismatch.
!
!
  if (sim_geometry == GRID_2DCYLINDRICAL  .and.  sim_lasersOrientation == 'Y') then
      call Driver_abortFlash ('[Simulation_init] ERROR: 2D cylindrical / lasers Y orientation not possible!')
  end if
!
!
!     ...Set the simulation error bars.
!
!
  call sim_setErrorBars ()
!
!
!     ...Set the tube center position and the radius around the center on which
!        the rays will be launched. Set also the tube dimensions.
!
!
  sim_R  =  3.0
  sim_xw =  5.0
  sim_yw =  5.0
  sim_zw =  5.0
  sim_xc = 10.0
  sim_yc = 10.0
  sim_zc = 10.0
!
!
!     ...Calculate the critical electron number density and set the electron number density
!        and the electron temperature (in Kelvin) at the center of the tube. Set the average
!        Z number throughout the tube as 1 and define the quadratic potential factor P.
!
!
  ratio  = sim_beamFrequency / ed_electronCharge
  sim_nc = ed_electronMass * PI * ratio * ratio
  sim_nw = sim_nc * 0.5
  sim_Tw = 10.0 * sim_keV2Kelvin
  sim_Z  = 1.0
  sim_A  = (sim_nc - sim_nw) * 0.04     ! 0.04 is the value of 1/(xc-xw)^2 in our case
!
!
!     ...Calculate the inverse-bremsstrahlung rate at the center of the tube.
!        The Coulomb Factor is set equal to 1..
!
!
  lnLambda = 1.0

  sim_nuw = ed_inverseBremsstrahlungRate (sim_Z,             &
                                          ed_electronCharge, &
                                          ed_electronMass,   &
                                          ed_Boltzmann,      &
                                          sim_Tw,            &
                                          sim_nw,            &
                                          sim_nc,            &
                                          lnLambda           )
!
!
!     ...Calculate the analytical decay factor for the ray's power.
!
!
  tcross   = (PI / (ed_speedOfLight + ed_speedOfLight)) * sqrt (sim_nc / sim_A)
  integral = sim_nuw * tcross * (1.0 + (sim_A * sim_R * sim_R) / (2.0 * sim_nw))

  sim_powerDecayFactor = exp (- integral)

  return
end subroutine Simulation_init
