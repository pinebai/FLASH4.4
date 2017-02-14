!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/sim_printResults
!!
!!  NAME 
!!
!!   sim_printResults
!!
!!  SYNOPSIS
!!
!!   sim_printResults ()
!!
!!  DESCRIPTION
!!
!!   This routine prints some general info about the simulation. It prints the obtained maximum
!!   square radial deviation to a .dat file. Only the master processor does the printout.
!!
!!***

subroutine sim_printResults ()

  use Simulation_data,        ONLY : sim_alpha,              &
                                     sim_baseName,           &
                                     sim_focalPointRadius,   &
                                     sim_globalMe,           &
                                     sim_maxDeltaRsqr,       &
                                     sim_nc,                 &
                                     sim_r0,                 &
                                     sim_Te0,                &
                                     sim_numRaysAnalyzed,    &
                                     sim_numRaysLaunched,    &
                                     sim_refinementLevel,    &
                                     sim_xTubeCenter,        &
                                     sim_yfocal,             &
                                     sim_zTubeCenter,        &
                                     sim_cellSizeX,          &
                                     sim_cellSizeY,          &
                                     sim_cellSizeZ,          &
                                     sim_vy,                 &
                                     sim_totalNumberOfBoxes, &
                                     sim_boxPowers,          &
                                     sim_boxRadius,          &
                                     sim_boxXvalue,          &
                                     sim_powersBox,          &
                                     sim_radialBox,          &
                                     sim_xvalueBox,          &
                                     sim_avgDeltaX,          &
                                     sim_rmsDeltaX,          &
                                     sim_sigmaDeltaX,        &
                                     sim_avgDeltaR,          &
                                     sim_rmsDeltaR,          &
                                     sim_sigmaDeltaR

  use EnergyDeposition_data,  ONLY : ed_Boltzmann,      &
                                     ed_electronCharge, &
                                     ed_electronMass,   &
                                     ed_speedOfLight
  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer  :: box
  integer  :: fileUnit
  integer  :: addEntries, nEntries, nTotalEntries
  integer  :: ray
  integer  :: ut_getFreeFileUnit

  real     :: power
  real     :: radius
  real     :: xvalue
!
!
!   ...Do the printout only on the master processor.
!
!
  if (sim_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (sim_baseName) // "Results.dat"

  open (fileUnit, file = fileName)
!
!
!   ...Print out the main title. 
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "      LASER RING: EXPONENTIAL TUBE POTENTIAL SIMULATION PRINTOUT"
  write (fileUnit,'(/)')
!
!
!     ...Print out the general info.
!
!
  write (fileUnit,'(1x,a,es20.12)') "                       Tube center x-coordinate (cm) = ",sim_xTubeCenter
  write (fileUnit,'(1x,a,es20.12)') "                       Tube center z-coordinate (cm) = ",sim_zTubeCenter
  write (fileUnit,'(1x,a,es20.12)') "                              Laser ring radius (cm) = ",sim_r0
  write (fileUnit,'(1x,a,es20.12)') "                             Focal y-coordinate (cm) = ",sim_yfocal
  write (fileUnit,'(1x,a,es20.12)') "                         Cell size x-coordinate (cm) = ",sim_cellSizeX
  write (fileUnit,'(1x,a,es20.12)') "                         Cell size y-coordinate (cm) = ",sim_cellSizeY
  write (fileUnit,'(1x,a,es20.12)') "                         Cell size z-coordinate (cm) = ",sim_cellSizeZ
  write (fileUnit,'(1x,a,es20.12)') "                  Critical electron density (#/cm^3) = ",sim_nc
  write (fileUnit,'(1x,a,a)'      ) "            Electron density at tube center (#/cm^3) = ","1/4 critical"
  write (fileUnit,'(1x,a,a)'      ) "            Electron density at tube edge   (#/cm^3) = ","1/2 critical"
  write (fileUnit,'(1x,a,es20.12)') "             Electron temperature at tube center (K) = ",sim_Te0
  write (fileUnit,'(1x,a,es20.12)') "                                  Boltzmann constant = ",ed_Boltzmann
  write (fileUnit,'(1x,a,es20.12)') "                                   Electron mass (g) = ",ed_electronMass
  write (fileUnit,'(1x,a,es20.12)') "                               Electron charge (esu) = ",ed_electronCharge
  write (fileUnit,'(1x,a,es20.12)') "                               Speed of Light (cm/s) = ",ed_speedOfLight
  write (fileUnit,'(1x,a,es20.12)') "                   Initial Ray Speed (in units of c) = ",sim_vy
  write (fileUnit,'(1x,a,a)'      ) "         k value in Tube Potential (in units of c^2) = ","1/8"
  write (fileUnit,'(1x,a,i2)'     ) "                    Exponent value in Tube Potential = ",sim_alpha
  write (fileUnit,'(1x,a,i2)'     ) "                         Simulation refinement level = ",sim_refinementLevel
  write (fileUnit,'(1x,a,i10)'    ) "                             Number of rays launched = ",sim_numRaysLaunched
  write (fileUnit,'(1x,a,i10)'    ) "                             Number of rays analyzed = ",sim_numRaysAnalyzed
!
!
!   ...Print out the results title. 
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                     RESULTS"
  write (fileUnit,'(/)')
!
!
!     ...Print out the results.
!
!
  write (fileUnit,'(1x,a,es20.10)') "                             Focal point radius (cm) = ",sim_focalPointRadius
  write (fileUnit,'(1x,a,es20.10)') "              Maximum square radial deviation (cm^2) = ",sim_maxDeltaRsqr
  write (fileUnit,'(1x,a,es20.10)') "                       Average radial deviation (cm) = ",sim_avgDeltaR
  write (fileUnit,'(1x,a,es20.10)') "           Root mean square of radial deviation (cm) = ",sim_rmsDeltaR
  write (fileUnit,'(1x,a,es20.10)') "         Standard deviation of radial deviation (cm) = ",sim_sigmaDeltaR
  write (fileUnit,'(1x,a,es20.10)') "                            Average x deviation (cm) = ",sim_avgDeltaX
  write (fileUnit,'(1x,a,es20.10)') "                Root mean square of x deviation (cm) = ",sim_rmsDeltaX
  write (fileUnit,'(1x,a,es20.10)') "              Standard deviation of x deviation (cm) = ",sim_sigmaDeltaX
!
!
!     ...Print out the power boxes.
!
!
  if (1 == 0) then

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    POWER BOX CONTENT"
  write (fileUnit,'(/)')

  nTotalEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_powersBox (box)
     nTotalEntries = nTotalEntries + nEntries

     if (nEntries > 0) then
!         write (fileUnit,'(1x,a,i8,a,i8,   a)') " Box Nr ",box," has ",nEntries," entries"
     else
!         write (fileUnit,'(1x,a,i8,a,i2,6x,a)') " Box Nr ",box," has ",nEntries," entries"
     end if
  end do

  write (fileUnit,'(/)')
  write (fileUnit,'(1x,a,i12)') " Total number of entries in all boxes = ",nTotalEntries

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    POWER BOX CONTENT PRINTOUT FOR DATA FILE (Power [erg]/ Count pairs)"
  write (fileUnit,'(/)')

  addEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_powersBox (box)
     addEntries = addEntries + nEntries

     power    = sim_boxPowers (box)

     if (real (addEntries) / real (nTotalEntries) > 0.9) then
         write (fileUnit,'(1x,a,es20.10      )') " power > 90% ",power
         addEntries = 0
     end if

     if (nEntries > 0) then
         write (fileUnit,'(1x,a,es20.10,i8   )') " ",power , nEntries
     else
!         write (fileUnit,'(1x,a,es20.10,i2,6x)') " ",power , nEntries
     end if
  end do

  end if
!
!
!     ...Print out the radial boxes.
!
!
!  if (1 == 0) then

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    RADIAL BOX CONTENT"
  write (fileUnit,'(/)')

  nTotalEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_radialBox (box)
     nTotalEntries = nTotalEntries + nEntries

     if (nEntries > 0) then
!         write (fileUnit,'(1x,a,i8,a,i8,   a)') " Box Nr ",box," has ",nEntries," entries"
     else
!         write (fileUnit,'(1x,a,i8,a,i2,6x,a)') " Box Nr ",box," has ",nEntries," entries"
     end if
  end do

  write (fileUnit,'(/)')
  write (fileUnit,'(1x,a,i12)') " Total number of entries in all boxes = ",nTotalEntries

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    RADIAL BOX CONTENT PRINTOUT FOR DATA FILE (Radius [units of cellsize]/ Count pairs)"
  write (fileUnit,'(/)')

  addEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_radialBox (box)
     addEntries = addEntries + nEntries

     radius   = sim_boxRadius (box)

     if (mod (nEntries,8) /= 0) then
         write (fileUnit,'(1x,a              )') " Unexpected number of rays! Symmetry violation! "
     end if

     if (real (addEntries) / real (nTotalEntries) > 0.9) then
         write (fileUnit,'(1x,a,es20.10      )') " radius > 90% ",radius/sim_cellSizeX
         addEntries = 0
     end if

     if (nEntries > 0) then
         write (fileUnit,'(1x,a,es20.10,i8   )') " ",radius/sim_cellSizeX , nEntries
     else
!         write (fileUnit,'(1x,a,es20.10,i2,6x)') " ",radius/sim_cellSizeX , nEntries
     end if
  end do

!  end if
!
!
!     ...Print out the x value boxes.
!
!
  if (1 == 0) then

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    X VALUE BOX CONTENT"
  write (fileUnit,'(/)')

  nTotalEntries = 0

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_xvalueBox (box)
     nTotalEntries = nTotalEntries + nEntries
     write (fileUnit,'(1x,a,i8,a,i8,a)') " Box Nr ",box," has ",nEntries," entries"
  end do
  write (fileUnit,'(/)')
  write (fileUnit,'(1x,a,i12)') " Total number of entries in all boxes = ",nTotalEntries

  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "    X VALUE BOX CONTENT PRINTOUT FOR DATA FILE (X / Count pairs)"
  write (fileUnit,'(/)')

  do box = 1, sim_totalNumberOfBoxes
     nEntries = sim_xvalueBox (box)
     xvalue   = sim_boxXvalue (box)
     write (fileUnit,'(1x,a,es20.10,i8)') " ",xvalue,nEntries
  end do

  end if
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_printResults
