!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/sim_analyzeDetectorsData
!!
!! NAME 
!!
!!  sim_analyzeDetectorsData
!!
!! SYNOPSIS
!!
!!  sim_analyzeDetectorsData ()
!!
!! DESCRIPTION
!!
!!  This routine performs some analysis on the proton detectors data.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_analyzeDetectorsData ()

  use Simulation_Data,         ONLY : sim_baseName,            &
                                      sim_globalComm,          &
                                      sim_globalMe,            &
                                      sim_magneticFluxDensity, &
                                      sim_protonCharge,        &
                                      sim_protonMass,          &
                                      sim_screenX,             &
                                      sim_screenY,             &
                                      sim_speedOfLight,        &
                                      sim_xCenter,             &
                                      sim_zCenter
                                      
  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  character (len = MAX_STRING_LENGTH) :: fileName
  character (len = 10               ) :: fileLineFormat = '(2es20.10)'

  logical :: fileExists

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
!     ...Only the master does the analysis.
!
!
  if (sim_globalMe /= MASTER_PE) then
      return
  end if
!
!
!     ...Get the data from the detector files.
!
!
  fileName = trim (sim_baseName) // "ProtonDetectorFile01"

  inquire (file = fileName, exist = fileExists)

  if (fileExists) then
      fileUnit = ut_getFreeFileUnit ()
      open (fileUnit, file = fileName)

      nProtons = 0
      do
        read (fileUnit,fileLineFormat,end=10) x,y
        nProtons = nProtons + 1
      end do

   10 write (*,*) ' Number of proton (x,y) pairs on detector screen = ',nProtons

      allocate (sim_screenX (1:nProtons))
      allocate (sim_screenY (1:nProtons))

      rewind (fileUnit)

      do n = 1,nProtons
         read (fileUnit,fileLineFormat) x,y
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

  c = sim_speedOfLight
  m = sim_protonMass
  B  = sim_magneticFluxDensity
  Qm = sim_protonCharge / sim_protonMass
  E  = 3.0               ! MeV
  E  = E * 1.60217657e-06 ! erg
  v0 = sqrt (1.0 - (m*c*c/(E + m*c*c))**2)  ! in units of c
  write (*,'(a,es14.5)') ' speed of protons (in c) = ',v0
  v0 = v0 * c
  Dy = 1.0

  Blimit = c * v0 / (Dy * Qm)
  write (*,'(a,es14.6)') ' B limit = ',Blimit

  thrDeltaR = c * v0 * (1.0 - sqrt (1.0 - (Dy * B * Qm / (c * v0)) ** 2)) / (B * Qm)

  write (*,'(a,es14.5)') ' deltaR (max) = ',maxDeltaR
  write (*,'(a,es14.5)') ' deltaR (min) = ',minDeltaR
  write (*,'(a,es14.5)') ' deltaR (avg) = ',sumDeltaR / real (nProtons)
  write (*,'(a,es14.5)') ' deltaR (thr) = ',thrDeltaR

!  do n = 1,50
!     E = real (n) * 1.60217657e-06
!     v = sqrt (1.0 - (m*c*c/(E + m*c*c))**2)
!     write (*,'(a,i2,es14.6)') ' MeV, v (in c) = ',n,v
!  end do

!
!
!     ...Ready!
!
!
  return
end subroutine sim_analyzeDetectorsData
