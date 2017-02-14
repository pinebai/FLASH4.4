!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_doAnalysis3DRec
!!
!!  NAME 
!!
!!   sim_doAnalysis3DRec
!!
!!  SYNOPSIS
!!
!!   sim_doAnalysis3DRec (logical (out) :: perfect)
!!
!!  DESCRIPTION
!!
!!   This routine performs the analysis on the launched rays for 3D geometries.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis3DRec (perfect)

  use Simulation_data,             ONLY : sim_globalComm,             &
                                          sim_globalMe,               &
                                          sim_nFocalRays,             &
                                          sim_nRaysMax,               &
                                          sim_powerDecayFactor,       &
                                          sim_rayPexit,               &
                                          sim_rayXexit,               &
                                          sim_rayYexit,               &
                                          sim_rayZexit,               &
                                          sim_rayFexitPercentError1,  &
                                          sim_rayFexitPercentError12, &
                                          sim_rayPexitPercentError1,  &
                                          sim_rayPexitPercentError12, &
                                          sim_symmetryTolerance,      &
                                          sim_xw,                     &
                                          sim_yw,                     &
                                          sim_zw,                     &
                                          sim_lasersOrientation

  use Driver_interface,            ONLY : Driver_abortFlash

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  include "Flash_mpi.h"

  logical, intent (out) :: perfect

  integer  :: error
  integer  :: ray

  real     :: FocalPercentErrorX (1:sim_nRaysMax)
  real     :: FocalPercentErrorY (1:sim_nRaysMax)
  real     :: FocalPercentErrorZ (1:sim_nRaysMax)
  real     :: PowerPercentError  (1:sim_nRaysMax)
!
!
!     ...Calculate the percentage error between the actual and analytical solution
!        for the rays exit coordinates and the rays exit power. Compare to the error
!        bars.
!
!        The 'perfect' indicator will be set initially to true as a default on all
!        processors. Only the master processor perfroms the testing and has thus the
!        possibility to change its status to false. The logical status of this indicator
!        will be broadcast from the master to all processors after the analysis has been
!        performed.
!
!
  perfect = .true.

  if (sim_globalMe == MASTER_PE) then

      do ray = 1,8
         FocalPercentErrorX (ray) = abs (sim_rayXexit (ray) - sim_xw) * 100.0 / sim_xw
         FocalPercentErrorY (ray) = abs (sim_rayYexit (ray) - sim_yw) * 100.0 / sim_yw
         FocalPercentErrorZ (ray) = abs (sim_rayZexit (ray) - sim_zw) * 100.0 / sim_zw
         PowerPercentError  (ray) = abs (sim_rayPexit (ray) - sim_powerDecayFactor) * 100.0 / sim_powerDecayFactor
      end do

      select case (sim_lasersOrientation)

        case ('XY')

          do ray = 1,4
             perfect = perfect .and. (FocalPercentErrorX (ray)  <= sim_rayFexitPercentError1)
             perfect = perfect .and. (FocalPercentErrorY (ray)  <= sim_rayFexitPercentError1)
             perfect = perfect .and. (PowerPercentError  (ray)  <= sim_rayPexitPercentError1)
          end do

          do ray = 5,8
             perfect = perfect .and. (FocalPercentErrorX (ray)  <= sim_rayFexitPercentError12)
             perfect = perfect .and. (FocalPercentErrorY (ray)  <= sim_rayFexitPercentError12)
             perfect = perfect .and. (PowerPercentError  (ray)  <= sim_rayPexitPercentError12)
          end do

        case ('XZ')

          do ray = 1,4
             perfect = perfect .and. (FocalPercentErrorX (ray)  <= sim_rayFexitPercentError1)
             perfect = perfect .and. (FocalPercentErrorZ (ray)  <= sim_rayFexitPercentError1)
             perfect = perfect .and. (PowerPercentError  (ray)  <= sim_rayPexitPercentError1)
          end do

          do ray = 5,8
             perfect = perfect .and. (FocalPercentErrorX (ray)  <= sim_rayFexitPercentError12)
             perfect = perfect .and. (FocalPercentErrorZ (ray)  <= sim_rayFexitPercentError12)
             perfect = perfect .and. (PowerPercentError  (ray)  <= sim_rayPexitPercentError12)
          end do

        case ('YZ')

          do ray = 1,4
             perfect = perfect .and. (FocalPercentErrorY (ray)  <= sim_rayFexitPercentError1)
             perfect = perfect .and. (FocalPercentErrorZ (ray)  <= sim_rayFexitPercentError1)
             perfect = perfect .and. (PowerPercentError  (ray)  <= sim_rayPexitPercentError1)
          end do

          do ray = 5,8
             perfect = perfect .and. (FocalPercentErrorY (ray)  <= sim_rayFexitPercentError12)
             perfect = perfect .and. (FocalPercentErrorZ (ray)  <= sim_rayFexitPercentError12)
             perfect = perfect .and. (PowerPercentError  (ray)  <= sim_rayPexitPercentError12)
          end do

        case default

          call Driver_abortFlash ('[sim_doAnalysis3DRec] ERROR: Lasers orientation bad!')

      end select

  end if
!
!
!     ...Broadcast the 'perfect' indicator.
!
!
  call MPI_Bcast (perfect,        &
                  1,              &
                  MPI_LOGICAL,    &
                  MASTER_PE,      &
                  sim_globalComm, &
                  error           )
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis3DRec
