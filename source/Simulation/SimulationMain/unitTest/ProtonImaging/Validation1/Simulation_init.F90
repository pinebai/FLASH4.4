!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/Simulation_init
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
!!  Initializes all the data specified in Simulation_data for the proton imaging
!!  unit test.
!!
!!***

subroutine Simulation_init()

  use Simulation_data

  use Grid_interface,              ONLY : Grid_getGeometry,    &
                                          Grid_getMinCellSizes

  use Driver_interface,            ONLY : Driver_abortFlash,  &
                                          Driver_getComm,     &
                                          Driver_getMype,     &
                                          Driver_getNumProcs

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: geometry
  integer :: lrefineMin, lrefineMax

  real    :: minCellSizes (1:MDIM)
!
!
!     ...Get the needed data.
!
!
  call Driver_getComm        (GLOBAL_COMM,                  sim_globalComm           )
  call Driver_getMype        (GLOBAL_COMM,                  sim_globalMe             )
  call Driver_getNumProcs    (GLOBAL_COMM,                  sim_globalNumProcs       )

  call RuntimeParameters_get ("lrefine_min",                lrefineMin               )
  call RuntimeParameters_get ("lrefine_max",                lrefineMax               )
  call RuntimeParameters_get ("basenm",                     sim_baseName             )
  call RuntimeParameters_get ("sim_xCenter",                sim_xCenter              )
  call RuntimeParameters_get ("sim_zCenter",                sim_zCenter              )
  call RuntimeParameters_get ("sim_clockwiseB",             sim_clockwiseB           )
  call RuntimeParameters_get ("sim_magneticFluxDensity",    sim_magneticFluxDensity  )
  call RuntimeParameters_get ("sim_printBlockVariables",    sim_printBlockVariables  )

  call PhysicalConstants_get ("speed of light",             sim_speedOfLight         )
  call PhysicalConstants_get ("electron charge",            sim_protonCharge         )
  call PhysicalConstants_get ("proton mass",                sim_protonMass           )
!
!
!       ...Catch bad data.
!
!
  if (lrefineMin /= lrefineMax) then
      call Driver_abortFlash ('[Simulation_init] ERROR: lrefine_min must match lrefine_max')
  end if

  if (sim_magneticFluxDensity < 0.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: magnetic flux density B < 0')
  end if

  sim_refinementLevel = lrefineMax

  call Grid_getGeometry (geometry)

  if ((NDIM /= 3) .or. (geometry /= CARTESIAN)  ) then
      call Driver_abortFlash ('Proton Imaging: the geometry must be 3D cartesian!')
  endif
!
!
!       ...Do some needed chores.
!
!
  sim_baseName = adjustl (sim_baseName)    ! to get ready to use 'trim'
!
!
!       ...Get the cell sizes for the entire domain. Since the maximum and
!          minimum refinement level are the same, there is only one cell size.
!
!
  call Grid_getMinCellSizes (minCellSizes)

  sim_cellSizeX = minCellSizes (IAXIS)
  sim_cellSizeY = minCellSizes (JAXIS)
  sim_cellSizeZ = minCellSizes (KAXIS)
!
!
!     ...Ready!
!
!
  return
end subroutine Simulation_init
