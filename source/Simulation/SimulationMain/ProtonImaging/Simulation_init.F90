!!****if* source/Simulation/SimulationMain/ProtonImaging/Simulation_init
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
!!  Initializes all the data specified in Simulation_data for the proton
!!  imaging run.
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
  call RuntimeParameters_get ("sim_printBlockVariables",    sim_printBlockVariables  )
!
!
!       ...Check and print some data.
!
!
  if (sim_globalMe == MASTER_PE) then
      write (*,*) ' lrefine min = ',lrefineMin
      write (*,*) ' lrefine max = ',lrefineMax
  end if

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
!       ...Get the cell sizes for the entire domain.
!
!
  call Grid_getMinCellSizes (minCellSizes)

  sim_cellSizeX = minCellSizes (IAXIS)
  sim_cellSizeY = minCellSizes (JAXIS)
  sim_cellSizeZ = minCellSizes (KAXIS)

  if (sim_globalMe == MASTER_PE) then
      write (*,*) ' min cell X = ',sim_cellSizeX
      write (*,*) ' min cell Y = ',sim_cellSizeY
      write (*,*) ' min cell Z = ',sim_cellSizeZ
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine Simulation_init
