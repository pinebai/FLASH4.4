!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_doAnalysis
!!
!! NAME 
!!
!!  sim_doAnalysis
!!
!! SYNOPSIS
!!
!!  sim_doAnalysis (logical (out) :: perfect)
!!
!! DESCRIPTION
!!
!!  This routine analyzes the ray exit results.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis (perfect)

  use Simulation_Data,         ONLY : sim_beamTargetRadius,     &
                                      sim_focalPointRadius,     &
                                      sim_globalComm,           &
                                      sim_globalMe,             &
                                      sim_maxDeltaPZ,           &
                                      sim_maxDeltaPZErrorBar,   &
                                      sim_maxDeltaRsqr,         &
                                      sim_maxDeltaRsqrErrorBar, &
                                      sim_numRaysAnalyzed,      &
                                      sim_numRaysLaunched,      &
                                      sim_percentFocus,         &
                                      sim_powerDecayFactor,     &
                                      sim_powerPartition,       &
                                      sim_refinementLevel,      &
                                      sim_xw,                   &
                                      sim_yfocal,               &
                                      sim_zw

  use EnergyDeposition_Data,   ONLY : ed_numberOfSavedRays,  &
                                      ed_raysSaved

  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  logical, intent (out) :: perfect

  integer  :: error
  integer  :: ray

  real     :: deltaRsqr, deltaPZ
  real     :: maxDeltaRsqr, maxDeltaPZ
  real     :: maxDeltaRsqrErrorBar, maxDeltaPZErrorBar
  real     :: rayPower
  real     :: rayX, rayY, rayZ
  real     :: Rb, Rf

  real     :: recvarray (1:2)
  real     :: sendarray (1:2)
!
!
!     ...Calculate the maximum square radial deviation and the maximum absolute partition
!        normalized power deviation for the current processor.
!
!
  maxDeltaRsqr = 0.0
  maxDeltaPZ   = 0.0

  do ray = 1,ed_numberOfSavedRays

     rayX     = ed_raysSaved (ray) % rayX
     rayY     = ed_raysSaved (ray) % rayY
     rayZ     = ed_raysSaved (ray) % rayZ
     rayPower = ed_raysSaved (ray) % rayPower

     if (rayY /= sim_yfocal) then
         call Driver_abortFlash ('[sim_doAnalysis] ERROR: ray Y exit different from yfocal')
     end if

     deltaRsqr = (rayX - sim_xw) ** 2 + (rayZ - sim_zw) ** 2
     deltaPZ   = abs (rayPower * sim_powerPartition - sim_powerDecayFactor)

     maxDeltaRsqr = max (deltaRsqr , maxDeltaRsqr)
     maxDeltaPZ   = max (deltaPZ   , maxDeltaPZ  )

  end do
!
!
!     ...Check, if these two values are below the error bars.
!
!
  maxDeltaRsqrErrorBar = sim_maxDeltaRsqrErrorBar (sim_refinementLevel)
  maxDeltaPZErrorBar   = sim_maxDeltaPZErrorBar   (sim_refinementLevel)

  perfect =      (maxDeltaRsqr <= maxDeltaRsqrErrorBar) &
           .and. (maxDeltaPZ   <= maxDeltaPZErrorBar  )
!
!
!     ...Determine the overall maximum values between all the processors and
!        communicate these values to all processor.
!
!
  sendarray (1) = maxDeltaRsqr
  sendarray (2) = maxDeltaPZ

  call MPI_Allreduce (sendarray,      &     ! the sending array on current processor P
                      recvarray,      &     ! the receiving array at the master
                      2,              &     ! number of sending elements on current processor P
                      FLASH_REAL,     &     ! of type real
                      MPI_Max,        &     ! kind of reduce operation
                      sim_globalComm, &     ! global communicator
                      error           )     ! error handle

  sim_maxDeltaRsqr = recvarray (1)
  sim_maxDeltaPZ   = recvarray (2)
!
!
!     ...Determine the total number of rays analyzed and give this info to all the
!        processors.
!
!
  call MPI_AllReduce (ed_numberOfSavedRays, &     ! the sending value on current processor P
                      sim_numRaysAnalyzed,  &     ! the receiving value at the master
                      1,                    &     ! number of sending elements on current processor P
                      FLASH_INTEGER,        &     ! of type real
                      MPI_Sum,              &     ! kind of reduce operation
                      sim_globalComm,       &     ! global communicator
                      error                 )     ! error handle

  perfect = perfect .and. (sim_numRaysAnalyzed == sim_numRaysLaunched)
!
!
!     ...Calculate the radius of the focal point. This is the maximum radius obtained for
!        a ray at the focal point. Also calculate the percent focus, defined as:
!
!                             (Rb - Rf)                   Rb = radius of beam target
!                             --------- x 100
!                                 Rb                      Rf = radius of focal point
!
!
  sim_focalPointRadius  = sqrt (sim_maxDeltaRsqr)

  Rb = sim_beamTargetRadius
  Rf = sim_focalPointRadius

  sim_percentFocus = ((Rb - Rf) / Rb) * 100.0
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
