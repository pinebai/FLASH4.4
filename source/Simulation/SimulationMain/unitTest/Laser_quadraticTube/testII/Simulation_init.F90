!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/Simulation_init
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

  if (     (NDIM == 3) .and. (geometry == CARTESIAN)  ) then
      sim_geometry = GRID_3DCARTESIAN
  else if ((NDIM == 2) .and. (geometry == CARTESIAN)  ) then
      sim_geometry = GRID_2DCARTESIAN
  else if ((NDIM == 2) .and. (geometry == CYLINDRICAL)) then
      sim_geometry = GRID_2DCYLINDRICAL
  else
      call Driver_abortFlash ('Laser quadratic tube unit test: unsupported geometry!')
  endif
!
!
!     ...Set the simulation error bars.
!
!
!  call sim_setErrorBars ()
!
!
!     ...Set the tube center position and the radius around the center on which
!        the rays will be launched.
!
!
  select case (sim_geometry)

  case (GRID_2DCARTESIAN , GRID_2DCYLINDRICAL)

      sim_xw =  5.0
      sim_zw =  0.0

  case (GRID_3DCARTESIAN)

      sim_xw =  5.0
      sim_zw =  5.0

  case default

    call Driver_abortFlash ('Laser quadratic tube: unsupported geometry (Simulation_init)!')

  end select
!
!
!     ...Calculate the critical electron number density and set the electron number density
!        and the electron temperature (in Kelvin) at the center of the tube. Set the average
!        Z number throughout the tube as 1 and define the quadratic potential factor P.
!
!
  ratio      = sim_beamFrequency / ed_electronCharge
  sim_nc     = ed_electronMass * PI * ratio * ratio
  sim_nw     = sim_nc * 0.5
  sim_Tw     = 10.0 * sim_keV2Kelvin
  sim_Z      = 1.0
  sim_A      = (sim_nc - sim_nw) / ((sim_xc - sim_xw) ** 2)
  sim_tcross = (PI / (ed_speedOfLight + ed_speedOfLight)) * sqrt (sim_nc / sim_A)
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
!     ...Calculate the power decay factor.
!
!
  sim_powerDecayFactor = exp (- sim_nuw * sim_tcross)

  return
end subroutine Simulation_init
