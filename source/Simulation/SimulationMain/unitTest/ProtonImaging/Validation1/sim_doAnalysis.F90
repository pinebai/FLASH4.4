!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/sim_doAnalysis
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
!!  This routine analyzes the proton imaging results. It gets the data recorded
!!  on the proton detector screen and compares it to the theoretically expected
!!  results.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis (perfect)

  use Simulation_Data,         ONLY : sim_baseName,            &
                                      sim_globalComm,          &
                                      sim_globalMe,            &
                                      sim_magneticFluxDensity, &
                                      sim_screenX,             &
                                      sim_screenY,             &
                                      sim_xCenter,             &
                                      sim_zCenter
                                      
  use ProtonImaging_Data,      ONLY : pi_beams,                 &
                                      pi_detectorLNwriteFormat, &
                                      pi_MeV2erg,               &
                                      pi_protonChargePerMass,   &
                                      pi_protonMass,            &
                                      pi_speedOfLight

  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  logical, intent (out) :: perfect

  character (len = MAX_STRING_LENGTH) :: fileName

  logical :: fileExists

  integer :: error
  integer :: fileRecordLength
  integer :: fileUnit
  integer :: n, nProtons
  integer :: ut_getFreeFileUnit

  real    :: c,m,E,v

  real    :: B, Blimit
  real    :: deltaR
  real    :: Dy
  real    :: minDeltaR, maxDeltaR
  real    :: r,x,y
  real    :: Qm
  real    :: sumDeltaR, thrDeltaR
  real    :: v0
!
!
!     ...The 'perfect' indicator will be set initially to false as a default on all
!        processors. Only the master processor perfroms the analysis and has thus the
!        possibility to change its status to true. The logical status of this indicator
!        will be broadcast from the master to all processors after the analysis has been
!        performed.
!
!
  perfect = .false.

  if (sim_globalMe == MASTER_PE) then

      fileName = trim (sim_baseName) // "ProtonDetectorFile01"

      inquire (file = fileName, exist = fileExists)

      if (fileExists) then
          fileUnit = ut_getFreeFileUnit ()
          open (fileUnit, file = fileName)

          nProtons = 0
          do
            read (fileUnit, pi_detectorLNwriteFormat, end=10) x,y
            nProtons = nProtons + 1
          end do

 10       allocate (sim_screenX (1:nProtons))
          allocate (sim_screenY (1:nProtons))

          rewind (fileUnit)

          do n = 1,nProtons
             read (fileUnit, pi_detectorLNwriteFormat) x,y
             sim_screenX (n) = x
             sim_screenY (n) = y
          end do
      else
          call Driver_abortFlash ('Proton detector file does not exist!')
      endif
!
!
!     ...Analyze the proton (x,y) pairs.
!
!
      perfect = .false.

      sumDeltaR = 0.0
      minDeltaR = huge (1.0)
      maxDeltaR = 0.0

      do n = 1,nProtons
         x = sim_screenX (n) * 3.0     ! screen magnification
         y = sim_screenY (n) * 3.0     ! screen magnification
         x = x - sim_xCenter
         y = y - sim_zCenter
         r = sqrt (x * x + y * y)
         deltaR = r - 0.5              ! subtract the beam radius

         sumDeltaR = sumDeltaR + deltaR
         minDeltaR = min (DeltaR, minDeltaR)
         maxDeltaR = max (DeltaR, maxDeltaR)
      end do

      c  = pi_speedOfLight
      m  = pi_protonMass
      B  = sim_magneticFluxDensity
      Qm = pi_protonChargePerMass
      E  = pi_beams (1) % protonEnergy            ! MeV
      E  = E * pi_MeV2erg                         ! erg
      v0 = pi_beams (1) % initialProtonSpeed      ! in cm/s
      Dy = 1.0

      Blimit = c * v0 / (Dy * Qm)
      write (*,'(a,es14.6)') ' B limit = ',Blimit

      thrDeltaR = c * v0 * (1.0 - sqrt (1.0 - (Dy * B * Qm / (c * v0)) ** 2)) / (B * Qm)

      write (*,'(a,es14.5)') ' deltaR (max) = ',maxDeltaR
      write (*,'(a,es14.5)') ' deltaR (min) = ',minDeltaR
      write (*,'(a,es14.5)') ' deltaR (avg) = ',sumDeltaR / real (nProtons)
      write (*,'(a,es14.5)') ' deltaR (thr) = ',thrDeltaR

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
end subroutine sim_doAnalysis
