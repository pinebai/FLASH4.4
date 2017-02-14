!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/Ring/sim_doAnalysis
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
                                      sim_maxDeltaP,            &
                                      sim_maxDeltaRsqr,         &
                                      sim_numRaysAnalyzed,      &
                                      sim_numRaysLaunched,      &
                                      sim_powerDecayFactor,     &
                                      sim_refinementLevel,      &
                                      sim_xw,                   &
                                      sim_yfocal,               &
                                      sim_zw,                   &
                                      sim_R,                    &
                                      sim_totalNumberOfBoxes,   &
                                      sim_boxPower,             &
                                      sim_boxRadius,            &
                                      sim_powerBox,             &
                                      sim_radialBox

  use EnergyDeposition_Data,   ONLY : ed_numberOfSavedRays,  &
                                      ed_raysSaved

  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  logical, intent (out) :: perfect

  integer  :: box, boxR, boxP
  integer  :: error
  integer  :: ray

  real     :: boxesPerUnitP, boxesPerUnitR      
  real     :: deltaRsqr, deltaP, deltaR
  real     :: maxDeltaRsqr, maxDeltaP
  real     :: maxMeasuredDeltaP
  real     :: powerBoxWidth, powerBoxWidthHalf  
  real     :: radialBoxWidth, radialBoxWidthHalf        
  real     :: rayPower
  real     :: rayX, rayY, rayZ

  real     :: recvarray (1:2)
  real     :: sendarray (1:2)
!
!
!     ...Calculate the maximum square radial deviation and the maximum absolute partition
!        normalized power deviation for the current processor.
!
!
  sim_powerBox  (:) = 0
  sim_radialBox (:) = 0
  sim_boxPower  (:) = 0.0
  sim_boxRadius (:) = 0.0

  maxMeasuredDeltaP   = 0.003

  powerBoxWidth       = maxMeasuredDeltaP / real (sim_totalNumberOfBoxes)
  powerBoxWidthHalf   = powerBoxWidth * 0.5
  boxesPerUnitP       = real (sim_totalNumberOfBoxes) / maxMeasuredDeltaP

  radialBoxWidth      = sim_R / real (sim_totalNumberOfBoxes)
  radialBoxWidthHalf  = radialBoxWidth * 0.5
  boxesPerUnitR       = real (sim_totalNumberOfBoxes) / sim_R

  do box = 1, sim_totalNumberOfBoxes
     sim_boxPower  (box) = powerBoxWidthHalf  + (box - 1) * powerBoxWidth
     sim_boxRadius (box) = radialBoxWidthHalf + (box - 1) * radialBoxWidth
  end do

  maxDeltaRsqr = 0.0
  maxDeltaP    = 0.0

  do ray = 1,ed_numberOfSavedRays

     rayX     = ed_raysSaved (ray) % rayX
     rayY     = ed_raysSaved (ray) % rayY
     rayZ     = ed_raysSaved (ray) % rayZ
     rayPower = ed_raysSaved (ray) % rayPower

     if (rayY /= sim_yfocal) then
         call Driver_abortFlash ('[sim_doAnalysis] ERROR: ray Y exit different from yfocal')
     end if

     deltaRsqr = (rayX - sim_xw) ** 2 + (rayZ - sim_zw) ** 2
     deltaP    = abs (rayPower - sim_powerDecayFactor)
     deltaR    = sqrt (deltaRsqr)

     maxDeltaRsqr = max (deltaRsqr , maxDeltaRsqr)
     maxDeltaP    = max (deltaP    , maxDeltaP   )

     boxP = int (ceiling (boxesPerUnitP * deltaP))
     boxP = max (1,boxP)                             ! safety net against box Nr = 0
     boxP = min (sim_totalNumberOfBoxes , boxP)      ! safety net against box Nr out of maximum range

     boxR = int (ceiling (boxesPerUnitR * deltaR))
     boxR = max (1,boxR)                             ! safety net against box Nr = 0
     boxR = min (sim_totalNumberOfBoxes , boxR)      ! safety net against box Nr out of maximum range

     sim_powerBox  (boxP) = sim_powerBox  (boxP) + 1
     sim_radialBox (boxR) = sim_radialBox (boxR) + 1

  end do
!
!
!     ...Determine the overall maximum values between all the processors and
!        communicate these values to all processor.
!
!
  sendarray (1) = maxDeltaRsqr
  sendarray (2) = maxDeltaP

  call MPI_Allreduce (sendarray,      &     ! the sending array on current processor P
                      recvarray,      &     ! the receiving array at the master
                      2,              &     ! number of sending elements on current processor P
                      FLASH_REAL,     &     ! of type real
                      MPI_Max,        &     ! kind of reduce operation
                      sim_globalComm, &     ! global communicator
                      error           )     ! error handle

  sim_maxDeltaRsqr = recvarray (1)
  sim_maxDeltaP    = recvarray (2)

  sim_focalPointRadius = sqrt (sim_maxDeltaRsqr)
!
!
!     ...Determine the total number of rays analyzed and give this info to all the
!        processors.
!
!
  call MPI_AllReduce (ed_numberOfSavedRays, &     ! the sending value on current processor P
                      sim_numRaysAnalyzed,  &     ! the receiving value at the master
                      1,                    &     ! number of sending elements on current processor P
                      FLASH_INTEGER,        &     ! of type integer
                      MPI_Sum,              &     ! kind of reduce operation
                      sim_globalComm,       &     ! global communicator
                      error                 )     ! error handle

  perfect = (sim_numRaysAnalyzed == sim_numRaysLaunched)
!
!
!     ...Sum up all power and radial box contents on all processors.
!
!
  call MPI_AllReduce (MPI_IN_PLACE,           &     ! the sending buffer on current processor P
                      sim_powerBox,           &     ! the receiving buffer (in place)
                      sim_totalNumberOfBoxes, &     ! number of sending elements on current processor P
                      FLASH_INTEGER,          &     ! of type integer
                      MPI_Sum,                &     ! kind of reduce operation
                      sim_globalComm,         &     ! global communicator
                      error                   )     ! error handle

  call MPI_AllReduce (MPI_IN_PLACE,           &     ! the sending buffer on current processor P
                      sim_radialBox,          &     ! the receiving buffer (in place)
                      sim_totalNumberOfBoxes, &     ! number of sending elements on current processor P
                      FLASH_INTEGER,          &     ! of type integer
                      MPI_Sum,                &     ! kind of reduce operation
                      sim_globalComm,         &     ! global communicator
                      error                   )     ! error handle
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
