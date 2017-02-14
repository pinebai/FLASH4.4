!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_printMainData
!!
!! NAME
!!
!!  pi_printMainData
!!
!! SYNOPSIS
!!
!!  call pi_printMainData ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints the main data of the proton imaging unit after initialization.
!!  The information is written out to a file named <basenm>ProtonImagingMainPrint.txt, where
!!  <basenm> is the runtime parameter for output file names.
!!
!!***

subroutine pi_printMainData ()

  use ProtonImaging_data,  ONLY : pi_3Din2DwedgeAngle,         &
                                  pi_3Din2DwedgeSlope,         &
                                  pi_badTiltingAxis,           &
                                  pi_baseName,                 &
                                  pi_cellStepTolerance,        &
                                  pi_cellWallThickness,        &
                                  pi_detectorLNwriteFormat,    &
                                  pi_domainErrorMarginX,       &
                                  pi_domainErrorMarginY,       &
                                  pi_domainErrorMarginZ,       &
                                  pi_domainTolerance,          &
                                  pi_flagDomainMissingProtons, &
                                  pi_globalMe,                 &
                                  pi_infiniteTime,             &
                                  pi_infiniteSpeed,            &
                                  pi_IOmaxBlockCrossingNumber, &
                                  pi_IOmaxPointsPerBlock,      &
                                  pi_IOmaxProtonCount,         &
                                  pi_IOnumberOfProtons2Plot,   &
                                  pi_IOprotonWriteModulo,      &
                                  pi_Joule2erg,                &
                                  pi_largestPositiveInteger,   &
                                  pi_largestPositiveReal,      &
                                  pi_microns2cm,               &
                                  pi_normalizedTolerance,      &
                                  pi_notSetInteger,            &
                                  pi_notSetReal,               &
                                  pi_orthogonalTolerance,      &
                                  pi_protonCharge,             &
                                  pi_protonChargePerMass,      &
                                  pi_protonMass,               &
                                  pi_screenProtonDiagnostics,  &
                                  pi_speedOfLight,             &
                                  pi_totalProtons2Launch,      &
                                  pi_unitRoundoff,             &
                                  pi_useIOprotonPlot,          &
                                  pi_xminDomain,               &
                                  pi_xmaxDomain,               &
                                  pi_yminDomain,               &
                                  pi_ymaxDomain,               &
                                  pi_zminDomain,               &
                                  pi_zmaxDomain

  implicit none
   
#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer :: fileUnit
  integer :: ut_getFreeFileUnit
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
  fileName = trim (pi_baseName) // "ProtonImagingMainPrint.txt"

  open (fileUnit, file = fileName)
!
!
!   ...Print out the title. 
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "   PROTON IMAGING MAIN PRINTOUT"
  write (fileUnit,'(/)')
!
!
!     ...Print the main data.
!
!
  write (fileUnit,'(1x,a,a)'      ) "  Simulation base name for output files = ",trim (pi_baseName)
  write (fileUnit,'(1x,a,a)'      ) "  Format for data lines for detector(s) = ",pi_detectorLNwriteFormat
  write (fileUnit,'(1x,a,L1)'     ) "         Diagnostics for screen protons = ",pi_screenProtonDiagnostics
  write (fileUnit,'(1x,a,L1)'     ) "  Abort run, if protons miss the domain = ",pi_flagDomainMissingProtons
  write (fileUnit,'(1x,a,L1)'     ) "      Are protons to be plotted via IO? = ",pi_useIOprotonPlot
  write (fileUnit,'(1x,a,i20)'    ) "         Largest positive integer value = ",pi_largestPositiveInteger
  write (fileUnit,'(1x,a,i20)'    ) "                  Not set integer value = ",pi_notSetInteger
  write (fileUnit,'(1x,a,i20)'    ) " Total number of protons to be launched = ",pi_totalProtons2Launch
  write (fileUnit,'(1x,a,i20)'    ) "       IO maximum block crossing number = ",pi_IOmaxBlockCrossingNumber
  write (fileUnit,'(1x,a,i20)'    ) "   IO maximum plotting points per block = ",pi_IOmaxPointsPerBlock
  write (fileUnit,'(1x,a,i20)'    ) "    IO maximum # of IO protons per node = ",pi_IOmaxProtonCount
  write (fileUnit,'(1x,a,i20)'    ) "     IO number of protons to be plotted = ",pi_IOnumberOfProtons2Plot
  write (fileUnit,'(1x,a,i20)'    ) "        IO every x-th proton is plotted = ",pi_IOprotonWriteModulo
  write (fileUnit,'(1x,a,es20.10)') "                     Not set real value = ",pi_notSetReal
  write (fileUnit,'(1x,a,es20.10)') "                            Domain xmin = ",pi_xminDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain xmax = ",pi_xmaxDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain ymin = ",pi_yminDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain ymax = ",pi_ymaxDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain zmin = ",pi_zminDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain zmax = ",pi_zmaxDomain
  write (fileUnit,'(1x,a,es20.10)') "    Domain relative computational error = ",pi_domainTolerance
  write (fileUnit,'(1x,a,es20.10)') "    Domain error margin in x-coordinate = ",pi_domainErrorMarginX
  write (fileUnit,'(1x,a,es20.10)') "    Domain error margin in y-coordinate = ",pi_domainErrorMarginY
  write (fileUnit,'(1x,a,es20.10)') "    Domain error margin in z-coordinate = ",pi_domainErrorMarginZ
  write (fileUnit,'(1x,a,es20.10)') "                    Cell wall thickness = ",pi_cellWallThickness
  write (fileUnit,'(1x,a,es20.10)') "  Cell step tolerance (unit: cell edge) = ",pi_cellStepTolerance
  write (fileUnit,'(1x,a,es20.10)') "             Bad tilting axis criterion = ",pi_badTiltingAxis
  write (fileUnit,'(1x,a,es20.10)') "      Normalization tolerance criterion = ",pi_normalizedTolerance
  write (fileUnit,'(1x,a,es20.10)') "         Orthogonal tolerance criterion = ",pi_orthogonalTolerance
  write (fileUnit,'(1x,a,es20.10)') "      Energy conversion factor J -> erg = ",pi_Joule2erg
  write (fileUnit,'(1x,a,es20.10)') " Length conversion factor microns -> cm = ",pi_microns2cm
  write (fileUnit,'(1x,a,es20.10)') "                  Speed of light (cm/s) = ",pi_speedOfLight
  write (fileUnit,'(1x,a,es20.10)') "                     Mass of proton (g) = ",pi_protonMass
  write (fileUnit,'(1x,a,es20.10)') "                 Charge of proton (esu) = ",pi_protonCharge
  write (fileUnit,'(1x,a,es20.10)') "      Charge per mass of proton (esu/g) = ",pi_protonChargePerMass
  write (fileUnit,'(1x,a,es20.10)') "                      Infinite time (s) = ",pi_infiniteTime
  write (fileUnit,'(1x,a,es20.10)') "                  Infinite speed (cm/s) = ",pi_infiniteSpeed
  write (fileUnit,'(1x,a,es20.10)') "        Unit roundoff (machine epsilon) = ",pi_unitRoundoff
  write (fileUnit,'(1x,a,es20.10)') "            Largest positive real value = ",pi_largestPositiveReal
  write (fileUnit,'(1x,a,es20.10)') "         3D in 2D wedge angle (radians) = ",pi_3Din2DwedgeAngle
  write (fileUnit,'(1x,a,es20.10)') "     3D in 2D wedge slope (m in y = mx) = ",pi_3Din2DwedgeSlope
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
end subroutine pi_printMainData
