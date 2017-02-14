!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/sim_doAnalysis
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

  use Simulation_Data,         ONLY : sim_focalPointRadius,     &
                                      sim_globalComm,           &
                                      sim_globalMe,             &
                                      sim_maxDeltaRsqr,         &
                                      sim_numRaysAnalyzed,      &
                                      sim_numRaysLaunched,      &
                                      sim_refinementLevel,      &
                                      sim_xTubeCenter,          &
                                      sim_yfocal,               &
                                      sim_zTubeCenter,          &
                                      sim_r0,                   &
                                      sim_totalNumberOfBoxes,   &
                                      sim_cellSizeX,            &
                                      sim_cellSizeY,            &
                                      sim_cellSizeZ,            &
                                      sim_boxPowers,            &
                                      sim_boxRadius,            &
                                      sim_boxXvalue,            &
                                      sim_powersBox,            &
                                      sim_radialBox,            &
                                      sim_xvalueBox,            &
                                      sim_avgDeltaX,            &
                                      sim_rmsDeltaX,            &
                                      sim_sigmaDeltaX,          &
                                      sim_avgDeltaR,            &
                                      sim_rmsDeltaR,            &
                                      sim_sigmaDeltaR

  use EnergyDeposition_Data,   ONLY : ed_numberOfSavedRays,  &
                                      ed_raysSaved

  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  logical, intent (out) :: perfect

  integer  :: box, boxP, boxR, boxX
  integer  :: error
  integer  :: m, n, nTotal
  integer  :: ray, rayTag

  real     :: averageR, averageRsqr
  real     :: averageX, averageXsqr
  real     :: boxesPerUnitP, boxesPerUnitR      
  real     :: deltaRsqr, deltaR, deltaP, deltaX
  real     :: maxDeltaRsqr
  real     :: maxMeasuredDeltaP, maxMeasuredDeltaR
  real     :: R, Rsqr
  real     :: X, Xsqr
  real     :: powersBoxWidth, powersBoxWidthHalf        
  real     :: radialBoxWidth, radialBoxWidthHalf        
  real     :: rayX, rayY, rayZ
  real     :: rayPower
  real     :: sumR, sumRsqr
  real     :: sumX, sumXsqr

  real     :: recvarray (1:2)
  real     :: sendarray (1:2)
!
!
!     ...Calculate the maximum square radial deviation for the current processor.
!
!
  sim_xvalueBox (:) = 0
  sim_powersBox (:) = 0
  sim_radialBox (:) = 0
  sim_boxPowers (:) = 0.0
  sim_boxRadius (:) = 0.0
  sim_boxXvalue (:) = 0.0

  maxMeasuredDeltaP   = 0.01
  maxMeasuredDeltaR   = sim_r0

  powersBoxWidth      = maxMeasuredDeltaP / real (sim_totalNumberOfBoxes)
  powersBoxWidthHalf  = powersBoxWidth * 0.5
  boxesPerUnitP       = real (sim_totalNumberOfBoxes) / maxMeasuredDeltaP
  radialBoxWidth      = maxMeasuredDeltaR / real (sim_totalNumberOfBoxes)
  radialBoxWidthHalf  = radialBoxWidth * 0.5
  boxesPerUnitR       = real (sim_totalNumberOfBoxes) / maxMeasuredDeltaR

  do box = 1, sim_totalNumberOfBoxes
     sim_boxPowers (box) = powersBoxWidthHalf + (box - 1) * powersBoxWidth
     sim_boxRadius (box) = radialBoxWidthHalf + (box - 1) * radialBoxWidth
     sim_boxXvalue (box) = radialBoxWidthHalf + (box - 1) * radialBoxWidth
  end do

  maxDeltaRsqr = 0.0

  do ray = 1,ed_numberOfSavedRays

     rayX     = ed_raysSaved (ray) % rayX
     rayY     = ed_raysSaved (ray) % rayY
     rayZ     = ed_raysSaved (ray) % rayZ
     rayPower = ed_raysSaved (ray) % rayPower
     rayTag   = ed_raysSaved (ray) % rayTag

!     if (rayY /= sim_yfocal) then
!         write (*,'(a,es24.16)') ' rayY       = ',rayY
!         write (*,'(a,es24.16)') ' sim_yfocal = ',sim_yfocal
!         call Driver_abortFlash ('[sim_doAnalysis] ERROR: ray Y exit different from yfocal')
!     end if

     deltaRsqr = (rayX - sim_xTubeCenter) ** 2 + (rayZ - sim_zTubeCenter) ** 2
     deltaR    = sqrt (deltaRsqr)
     deltaX    = abs (rayX - sim_xTubeCenter)
     deltaP    = abs (rayPower - 0.5)

     maxDeltaRsqr = max (deltaRsqr , maxDeltaRsqr)

     boxP = int (ceiling (boxesPerUnitP * deltaP))
     boxP = max (1,boxP)                             ! safety net against box Nr = 0
     boxP = min (sim_totalNumberOfBoxes , boxP)      ! safety net against box Nr out of maximum range

     boxR = int (ceiling (boxesPerUnitR * deltaR))
     boxR = max (1,boxR)                             ! safety net against box Nr = 0
     boxR = min (sim_totalNumberOfBoxes , boxR)      ! safety net against box Nr out of maximum range

     boxX = int (ceiling (boxesPerUnitR * deltaX))
     boxX = max (1,boxX)                             ! safety net against box Nr = 0
     boxX = min (sim_totalNumberOfBoxes , boxX)      ! safety net against box Nr out of maximum range

!     if (deltaR/sim_cellSizeX < 1.2) then
         sim_powersBox (boxP) = sim_powersBox (boxP) + 1
!     end if

     sim_radialBox (boxR) = sim_radialBox (boxR) + 1
     sim_xvalueBox (boxX) = sim_xvalueBox (boxX) + 1

  end do
!
!
!     ...Determine the overall maximum values between all the processors and
!        communicate these values to all processor.
!
!
  sendarray (1) = maxDeltaRsqr
  sendarray (2) = 0.0

  call MPI_Allreduce (sendarray,      &     ! the sending array on current processor P
                      recvarray,      &     ! the receiving array at the master
                      2,              &     ! number of sending elements on current processor P
                      FLASH_REAL,     &     ! of type real
                      MPI_Max,        &     ! kind of reduce operation
                      sim_globalComm, &     ! global communicator
                      error           )     ! error handle

  sim_maxDeltaRsqr = recvarray (1)

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
!     ...Sum up all power box contents on all processors.
!
!
  call MPI_AllReduce (MPI_IN_PLACE,           &     ! the sending buffer on current processor P
                      sim_powersBox,          &     ! the receiving buffer (in place)
                      sim_totalNumberOfBoxes, &     ! number of sending elements on current processor P
                      FLASH_INTEGER,          &     ! of type integer
                      MPI_Sum,                &     ! kind of reduce operation
                      sim_globalComm,         &     ! global communicator
                      error                   )     ! error handle
!
!
!     ...Sum up all radial box contents on all processors.
!
!
  call MPI_AllReduce (MPI_IN_PLACE,           &     ! the sending buffer on current processor P
                      sim_radialBox,          &     ! the receiving buffer (in place)
                      sim_totalNumberOfBoxes, &     ! number of sending elements on current processor P
                      FLASH_INTEGER,          &     ! of type integer
                      MPI_Sum,                &     ! kind of reduce operation
                      sim_globalComm,         &     ! global communicator
                      error                   )     ! error handle
!
!
!     ...Sum up all xvalue box contents on all processors.
!
!
  call MPI_AllReduce (MPI_IN_PLACE,           &     ! the sending buffer on current processor P
                      sim_xvalueBox,          &     ! the receiving buffer (in place)
                      sim_totalNumberOfBoxes, &     ! number of sending elements on current processor P
                      FLASH_INTEGER,          &     ! of type integer
                      MPI_Sum,                &     ! kind of reduce operation
                      sim_globalComm,         &     ! global communicator
                      error                   )     ! error handle
!
!
!     ...Calculate the root mean square and standard deviation values.
!
!
  nTotal  = 0
  sumX    = 0.0
  sumXsqr = 0.0
  sumR    = 0.0
  sumRsqr = 0.0

  do box = 1, sim_totalNumberOfBoxes
     m       = sim_xvalueBox (box)
     n       = sim_radialBox (box)
     nTotal  = nTotal + n
     R       = sim_boxRadius (box)
     Rsqr    = R * R
     X       = sim_boxXvalue (box)
     Xsqr    = X * X
     sumX    = sumX    + m * X
     sumXsqr = sumXsqr + m * Xsqr
     sumR    = sumR    + n * R
     sumRsqr = sumRsqr + n * Rsqr
  end do

  averageX    = sumX    / real (nTotal)
  averageXsqr = sumXsqr / real (nTotal)
  averageR    = sumR    / real (nTotal)
  averageRsqr = sumRsqr / real (nTotal)

  sim_avgDeltaR = averageR
  sim_avgDeltaX = averageX
  sim_rmsDeltaR = sqrt (averageRsqr)
  sim_rmsDeltaX = sqrt (averageXsqr)

  sumRsqr = 0.0
  sumXsqr = 0.0

  do box = 1, sim_totalNumberOfBoxes
     n       = sim_radialBox (box)
     R       = sim_boxRadius (box)
     sumRsqr = sumRsqr + n * (R - averageR) * (R - averageR)
     m       = sim_xvalueBox (box)
     X       = sim_boxXvalue (box)
     sumXsqr = sumXsqr + m * (X - averageX) * (X - averageX)
  end do

  averageRsqr = sumRsqr / real (nTotal)
  averageXsqr = sumXsqr / real (nTotal)

  sim_sigmaDeltaR = sqrt (averageRsqr)
  sim_sigmaDeltaX = sqrt (averageXsqr)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
