!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printMainData
!!
!! NAME
!!
!!  ed_printMainData
!!
!! SYNOPSIS
!!
!!  call ed_printMainData ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints the main data of the laser unit after initialization.
!!  The information is written out to a file named <basenm>LaserMainDataPrint.txt, where
!!  <basenm> is the runtime parameter for output file names.
!!
!!***

subroutine ed_printMainData ()

  use EnergyDeposition_data,  ONLY : ed_badTorsionAxis,         &
                                     ed_baseName,               &
                                     ed_cellStepTolerance,      &
                                     ed_cellWallThickness,      &
                                     ed_domainErrorMarginX,     &
                                     ed_domainErrorMarginY,     &
                                     ed_domainErrorMarginZ,     &
                                     ed_domainTolerance,        &
                                     ed_electronCharge,         &
                                     ed_electronMass,           &
                                     ed_energyProfileFileName,  &
                                     ed_globalMe,               &
                                     ed_Joule2erg,              &
                                     ed_infinitePower,          &
                                     ed_infiniteTime,           &
                                     ed_infiniteSpeed,          &
                                     ed_largestPositiveInteger, &
                                     ed_largestPositiveReal,    &
                                     ed_laser3Din2DwedgeAngle,  &
                                     ed_laser3Din2DwedgeSlope,  &
                                     ed_microns2cm,             &
                                     ed_normalizedTolerance,    &
                                     ed_notSetInteger,          &
                                     ed_notSetReal,             &
                                     ed_orthogonalTolerance,    &
                                     ed_rayZeroPower,           &
                                     ed_speedOfLight,           &
                                     ed_unitRoundoff,           &
                                     ed_xminDomain,             &
                                     ed_xmaxDomain,             &
                                     ed_yminDomain,             &
                                     ed_ymaxDomain,             &
                                     ed_zminDomain,             &
                                     ed_zmaxDomain

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
  if (ed_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (ed_baseName) // "LaserMainDataPrint.txt"

  open (fileUnit, file = fileName)
!
!
!   ...Print out the title. 
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "   LASER MAIN PRINTOUT"
  write (fileUnit,'(/)')
!
!
!     ...Print the main data.
!
!
  write (fileUnit,'(1x,a,es20.10)') "                  Speed of light (cm/s) = ",ed_speedOfLight
  write (fileUnit,'(1x,a,es20.10)') "                   Mass of electron (g) = ",ed_electronMass
  write (fileUnit,'(1x,a,es20.10)') "               Charge of electron (esu) = ",ed_electronCharge
  write (fileUnit,'(1x,a,es20.10)') "                      Infinite time (s) = ",ed_infiniteTime
  write (fileUnit,'(1x,a,es20.10)') "                  Infinite speed (cm/s) = ",ed_infiniteSpeed
  write (fileUnit,'(1x,a,es20.10)') "                 Infinite power (erg/s) = ",ed_infinitePower
  write (fileUnit,'(1x,a,es20.10)') "        Unit roundoff (machine epsilon) = ",ed_unitRoundoff
  write (fileUnit,'(1x,a,es20.10)') "            Largest positive real value = ",ed_largestPositiveReal
  write (fileUnit,'(1x,a,i20)'    ) "         Largest positive integer value = ",ed_largestPositiveInteger
  write (fileUnit,'(1x,a,i20)'    ) "                  Not set integer value = ",ed_notSetInteger
  write (fileUnit,'(1x,a,es20.10)') "                     Not set real value = ",ed_notSetReal
  write (fileUnit,'(1x,a,es20.10)') "                            Domain xmin = ",ed_xminDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain xmax = ",ed_xmaxDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain ymin = ",ed_yminDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain ymax = ",ed_ymaxDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain zmin = ",ed_zminDomain
  write (fileUnit,'(1x,a,es20.10)') "                            Domain zmax = ",ed_zmaxDomain
  write (fileUnit,'(1x,a,es20.10)') "    Domain relative computational error = ",ed_domainTolerance
  write (fileUnit,'(1x,a,es20.10)') "    Domain error margin in x-coordinate = ",ed_domainErrorMarginX
  write (fileUnit,'(1x,a,es20.10)') "    Domain error margin in y-coordinate = ",ed_domainErrorMarginY
  write (fileUnit,'(1x,a,es20.10)') "    Domain error margin in z-coordinate = ",ed_domainErrorMarginZ
  write (fileUnit,'(1x,a,es20.10)') "                    Cell wall thickness = ",ed_cellWallThickness
  write (fileUnit,'(1x,a,es20.10)') "  Cell step tolerance (unit: cell edge) = ",ed_cellStepTolerance
  write (fileUnit,'(1x,a,es20.10)') "             Bad torsion axis criterion = ",ed_badTorsionAxis
  write (fileUnit,'(1x,a,es20.10)') "      Normalization tolerance criterion = ",ed_normalizedTolerance
  write (fileUnit,'(1x,a,es20.10)') "         Orthogonal tolerance criterion = ",ed_orthogonalTolerance
  write (fileUnit,'(1x,a,es20.10)') "      Energy conversion factor J -> erg = ",ed_Joule2erg
  write (fileUnit,'(1x,a,es20.10)') " Length conversion factor microns -> cm = ",ed_microns2cm
  write (fileUnit,'(1x,a,es20.10)') "     Ray power is zero if below (erg/s) = ",ed_rayZeroPower
  write (fileUnit,'(1x,a,a)'      ) "  Simulation base name for output files = ",trim (ed_baseName)
  write (fileUnit,'(1x,a,a)'      ) "               Energy profile file name = ",trim (ed_energyProfileFileName)
  write (fileUnit,'(1x,a,es20.10)') "   Laser 3D in 2D wedge angle (radians) = ",ed_laser3Din2DwedgeAngle
  write (fileUnit,'(1x,a,es20.10)') " Laser 3D in 2D wedge slope (in y = mx) = ",ed_laser3Din2DwedgeSlope
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
end subroutine ed_printMainData
