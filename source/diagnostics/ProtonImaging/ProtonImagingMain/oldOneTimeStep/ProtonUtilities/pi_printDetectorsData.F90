!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_printDetectorsData
!!
!! NAME
!!
!!  pi_printDetectorsData
!!
!! SYNOPSIS
!!
!!  call pi_printDetectorsData ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the generated data
!!  for all detectors to a text file. The information is written out to a file named
!!  <basenm>ProtonImagingDetectors.txt, where <basenm> is the runtime parameter for
!!  output file names.
!!
!!***

subroutine pi_printDetectorsData ()

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_baseName,         &
                                  pi_detectors,        &
                                  pi_globalMe,         &
                                  pi_numberOfDetectors

  implicit none
   
#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  character (len =  MAX_STRING_LENGTH) :: fileName
  character (len =  1                ) :: tiltAxis

  integer :: alignbeamNr
  integer :: detector
  integer :: fileUnit
  integer :: ut_getFreeFileUnit

  real    :: centerX, centerY, centerZ
  real    :: dist2BC
  real    :: maxDev
  real    :: nx, ny, nz
  real    :: pinholeDist, pinholeRadius
  real    :: pinholeX, pinholeY, pinholeZ
  real    :: sideLength
  real    :: tiltAngle
  real    :: uXx, uXy, uXz
  real    :: uYx, uYy, uYz

  real    :: comp (1:9)
!
!
!   ...Do the printout only on the master processor.
!
!
  if (pi_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (pi_baseName) // "ProtonImagingDetectors.txt"

  open (fileUnit, file = fileName)
!
!
!     ...Loop over all detectors.
!
!
  do detector = 1, pi_numberOfDetectors

     alignbeamNr   = pi_detectors (detector) % alignWRTbeamNr
     uXx           = pi_detectors (detector) % axisXunitX
     uXy           = pi_detectors (detector) % axisXunitY
     uXz           = pi_detectors (detector) % axisXunitZ
     uYx           = pi_detectors (detector) % axisYunitX
     uYy           = pi_detectors (detector) % axisYunitY
     uYz           = pi_detectors (detector) % axisYunitZ
     centerX       = pi_detectors (detector) % centerX
     centerY       = pi_detectors (detector) % centerY
     centerZ       = pi_detectors (detector) % centerZ
     dist2BC       = pi_detectors (detector) % distance2beamCapsule
     nx            = pi_detectors (detector) % normalX
     ny            = pi_detectors (detector) % normalY
     nz            = pi_detectors (detector) % normalZ
     pinholeDist   = pi_detectors (detector) % pinholeDist2Det
     pinholeRadius = pi_detectors (detector) % pinholeRadius
     pinholeX      = pi_detectors (detector) % pinholeX
     pinholeY      = pi_detectors (detector) % pinholeY
     pinholeZ      = pi_detectors (detector) % pinholeZ
     sideLength    = pi_detectors (detector) % sideLength
     tiltAngle     = pi_detectors (detector) % sideTiltingAngle
     tiltAxis      = pi_detectors (detector) % sideTiltingAxis

     write (fileUnit,'(/)'        )
     write (fileUnit,'(a,i2)'     ) "               PROTON DETECTOR NR ",detector
     write (fileUnit,'(/)'        )

     if (alignbeamNr > 0) then

     write (fileUnit,'(a,i20)'    ) " Detector is aligned wrt beam nr                = ", alignbeamNr
     write (fileUnit,'(a,es20.12)') " Detector distance from beam capsule center     = ", dist2BC

     end if

     write (fileUnit,'(a,es20.12)') " Center global x-coordinate                     = ", centerX
     write (fileUnit,'(a,es20.12)') " Center global y-coordinate                     = ", centerY
     write (fileUnit,'(a,es20.12)') " Center global z-coordinate                     = ", centerZ
     write (fileUnit,'(a,es20.12)') " Detector square side length (cm)               = ", sideLength
     write (fileUnit,'(a,es20.12)') " Side tilting angle from tilting axis (degrees) = ", tiltAngle
     write (fileUnit,'(a,a1)'     ) " Global coordinate axis used for side tilting   = ", tiltAxis
     write (fileUnit,'(a,es20.12)') " Normal unit vector local x-coordinate          = ", nx
     write (fileUnit,'(a,es20.12)') " Normal unit vector local y-coordinate          = ", ny
     write (fileUnit,'(a,es20.12)') " Normal unit vector local z-coordinate          = ", nz
     write (fileUnit,'(a,es20.12)') " Axis X unit vector local x-coordinate          = ", uXx
     write (fileUnit,'(a,es20.12)') " Axis X unit vector local y-coordinate          = ", uXy
     write (fileUnit,'(a,es20.12)') " Axis X unit vector local z-coordinate          = ", uXz
     write (fileUnit,'(a,es20.12)') " Axis Y unit vector local x-coordinate          = ", uYx
     write (fileUnit,'(a,es20.12)') " Axis Y unit vector local y-coordinate          = ", uYy
     write (fileUnit,'(a,es20.12)') " Axis Y unit vector local z-coordinate          = ", uYz

     if (pinholeRadius > 0.0) then

     write (fileUnit,'(a,es20.12)') " Pinhole center global x-coordinate             = ", pinholeX
     write (fileUnit,'(a,es20.12)') " Pinhole center global y-coordinate             = ", pinholeY
     write (fileUnit,'(a,es20.12)') " Pinhole center global z-coordinate             = ", pinholeZ
     write (fileUnit,'(a,es20.12)') " Pinhole center distance to detector center     = ", pinholeDist
     write (fileUnit,'(a,es20.12)') " Pinhole radius (cm)                            = ", pinholeRadius

     end if
!
!
!     ...Print out maximum deviation from right hand rule vector triple (uX,uY,n).
!
!        We should have: uY x uX - n = 0 (3 components)
!                        uX x n - uY = 0 (3 components)
!                        uY x n + uX = 0 (3 components)
!
!        The absolute maximum deviation from 0 is recorded for all 9 components and printed.
!
!
     comp (1) = abs (uYy * uXz - uYz * uXy - nx)
     comp (2) = abs (uYz * uXx - uYx * uXz - ny)
     comp (3) = abs (uYx * uXy - uYy * uXx - nz)

     comp (4) = abs (uXy * nz - uXz * ny - uYx)
     comp (5) = abs (uXz * nx - uXx * nz - uYy)
     comp (6) = abs (uXx * ny - uXy * nx - uYz)

     comp (7) = abs (uYy * nz - uYz * ny + uXx)
     comp (8) = abs (uYz * nx - uYx * nz + uXy)
     comp (9) = abs (uYx * ny - uYy * nx + uXz)

     maxDev = maxval (comp (1:9))

     write (fileUnit,'(a,es20.12)') " Max dev of (uX,uY,n) from right hand rule      = ", maxDev

  end do
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!  
  return
end subroutine pi_printDetectorsData
