!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_checkReuseDepo
!!
!! NAME
!!
!!  ed_checkReuseDepo
!!
!! SYNOPSIS
!!
!!  call ed_checkReuseDepo (logical(OUT) :: reuseDepo)
!!
!! DESCRIPTION
!!
!!  This routine is called each cycle to check to see whether the result of a previous
!!  energy deposition computation, saved in the "depo" variable, should be reused for the
!!  current time step.
!!
!!  It is time to write ray data if:
!!
!!        1. The run-time parameter ed_useLaserIO is true
!!        2. This is NOT the second half of the time step when split driver is in use
!!        3. This is the cycle immediately following a plot
!!
!! ARGUMENTS
!!
!!  passSplitDriver : indicates first/second half of time step for split driver
!!
!!***

subroutine ed_checkReuseDepo (reuseDepo,laserIsOn)

  use EnergyDeposition_data,  ONLY : ed_laserIOWrite, &
                                     ed_currentStepNumber, &
                                     ed_depoReuseMaxSteps, &
                                     ed_depoVarValid,      &
                                     ed_savedDepoStepNumber

  implicit none

  logical, intent(OUT) :: reuseDepo
  logical, intent(in)  :: laserIsOn

  reuseDepo = .FALSE.

  if (.NOT. ed_depoVarValid) return

  ed_depoVarValid = .FALSE.

  if (.NOT. laserIsOn) return
  if (ed_depoReuseMaxSteps < 0) return
  if (ed_laserIOWrite) return

  if (ed_currentStepNumber - ed_savedDepoStepNumber .LE. ed_depoReuseMaxSteps) then
     ed_depoVarValid = .TRUE.
     reuseDepo = .TRUE.
  end if

  return
end subroutine ed_checkReuseDepo
