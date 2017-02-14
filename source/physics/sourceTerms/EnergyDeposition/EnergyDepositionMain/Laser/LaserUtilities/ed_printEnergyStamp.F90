!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printEnergyStamp
!!
!! NAME
!!
!!  ed_printEnergyStamp
!!
!! SYNOPSIS
!!
!!  call ed_printEnergyStamp (real, intent (in) :: timeStep,
!!                            real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Prints out the current energy statistics of the laser:
!!
!!        i) total laser energy pumped into the domain so far
!!       ii) total laser energy exiting the domain unused so far
!!      iii) laser energy pumped into the domain at the current timestep
!!       iv) laser energy exiting the domain unused at the current timestep
!!
!! ARGUMENTS
!!
!!  timeStep       : Current time step value
!!  timeSimulation : Current simulation time
!!
!! NOTES
!!
!!  Only the master processor writes to the energy profile file. Before this is
!!  done, a check is made, if the numbers make sense. There can be no more energy
!!  leaving the system than has been put into.
!!
!!***

subroutine ed_printEnergyStamp (timeStep, timeSimulation)

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_currentStepNumber,     &
                                     ed_energyInTimestep,      &
                                     ed_energyInTotal,         &
                                     ed_energyOutTimestep,     &
                                     ed_energyOutTotal,        &
                                     ed_energyProfileFileUnit, &
                                     ed_globalComm,            &
                                     ed_globalMe,              &
                                     ed_normalizedTolerance

  use Logfile_interface,      ONLY : Logfile_stampMessage

  implicit none
   
#include "Flash.h"
#include "constants.h"

  include "Flash_mpi.h"

  real, intent (in) :: timeStep
  real, intent (in) :: timeSimulation

  logical       :: badTimestepEnergies
  logical       :: badTotalEnergies
  logical, save :: firstTime = .true.

  integer :: error
  integer :: fileUnit

  real    :: energyLocal (1:2)
  real    :: energySum   (1:2)

  character(len=MAX_STRING_LENGTH) :: errmsg
!
!
!   ...Collect the time step energies on the master.
!
!
  energyLocal (1) = ed_energyInTimestep
  energyLocal (2) = ed_energyOutTimestep

  call MPI_Reduce (energyLocal,   &
                   energySum,     &
                   2,             &
                   FLASH_REAL,    &
                   MPI_Sum,       &
                   MASTER_PE,     &
                   ed_globalComm, &
                   error          )
!
!
!   ...Do the energy stamp printout only on the master processor.
!
!
  if (ed_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Update the accumulation variables.
!
!
  ed_energyInTimestep  = energySum (1)
  ed_energyOutTimestep = energySum (2)
  ed_energyInTotal     = ed_energyInTotal  + ed_energyInTimestep
  ed_energyOutTotal    = ed_energyOutTotal + ed_energyOutTimestep
!
!
!   ...Check the energy numbers for energy violation.
!
!
  badTimestepEnergies = ((ed_energyOutTimestep - ed_energyInTimestep)/ed_energyInTimestep > ed_normalizedTolerance)
  badTotalEnergies    = ((ed_energyOutTotal - ed_energyInTotal)/ed_energyInTotal > ed_normalizedTolerance)

  if (badTimestepEnergies) then
      call Logfile_stampMessage("[ed_printEnergyStamp] Error: Energy out of domain > Energy into domain!")
      write(errmsg,*) "[ed_printEnergyStamp]   ed_energyOutTimestep = ", ed_energyOutTimestep
      call Logfile_stampMessage(trim(errmsg))
      write(errmsg,*) "[ed_printEnergyStamp]   ed_energyInTimestep  = ", ed_energyInTimestep
      call Logfile_stampMessage(trim(errmsg))

      call Driver_abortFlash ("ed_printEnergyStamp: dE out of domain > dE into domain!")
  end if

  if (badTotalEnergies) then
      call Logfile_stampMessage("[ed_printEnergyStamp] Error: Energy out of domain > Energy into domain!")
      write(errmsg,*) "[ed_printEnergyStamp]   ed_energyOutTotal = ", ed_energyOutTotal
      call Logfile_stampMessage(trim(errmsg))
      write(errmsg,*) "[ed_printEnergyStamp]   ed_energyInTotal  = ", ed_energyInTotal
      call Logfile_stampMessage(trim(errmsg))

      call Driver_abortFlash ("ed_printEnergyStamp: Energy out of domain > Energy into domain!")
  end if
!
!
!     ...Print the header if first time calling the routine.
!
!
  fileUnit = ed_energyProfileFileUnit

  if (firstTime) then
      write (fileUnit,'(7a)') "#Step "              , &
                              "    Time (s)"        , &
                              "     dt (s) "        , &
                              "     Energy in (erg)", &
                              "    Energy out (erg)", &
                              "       dE in (erg)  ", &
                              "      dE out (erg)  "
      write (fileUnit,'(a)')  "#"
      firstTime = .false.
  end if
!
!
!     ...Print the energy stamp.
!
!
  write (fileUnit,'(i5,1x,2es12.4,4es20.10)') ed_currentStepNumber, &
                                              timeSimulation,       &
                                              timeStep,             &
                                              ed_energyInTotal,     &
                                              ed_energyOutTotal,    &
                                              ed_energyInTimestep,  &
                                              ed_energyOutTimestep

!
!
!    ...Ready!
!
!  
  return
end subroutine ed_printEnergyStamp
