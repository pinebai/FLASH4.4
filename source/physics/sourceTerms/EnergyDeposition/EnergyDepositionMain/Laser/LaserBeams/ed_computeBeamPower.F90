!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_computeBeamPower
!!
!! NAME
!!
!!  ed_computeBeamPower
!!
!! SYNOPSIS
!!
!!  call ed_computeBeamPower (real,    intent (in)  :: timeStep,
!!                            real,    intent (in)  :: timeSimulation,
!!                            integer, intent (in)  :: pulseNumber, 
!!                            real,    intent (out) :: beamPower)
!!
!! DESCRIPTION
!!
!!  Computes the average power of the beam between the time frame:
!!
!!             timeSimulation --> timeSimulation + timeStep
!!
!!  from the pulse data supplied. Zero beam power is returned, if the time
!!  frame is outside the pulse defining time range.
!!
!! ARGUMENTS
!!
!!  timeStep       : current time step value
!!  timeSimulation : current simulation time
!!  pulseNumber    : pulse identification number
!!  beamPower      : the average beam power obtained
!!
!! NOTES
!!
!!***

subroutine ed_computeBeamPower (timeStep, timeSimulation, pulseNumber, beamPower)

  use EnergyDeposition_data,    ONLY : ed_numberOfPulses,        &
                                       ed_pulseNumberOfSections, &
                                       ed_pulses

  use Driver_interface,         ONLY : Driver_abortFlash
  
  implicit none

#include "Flash.h"
#include "constants.h"

  real,    intent (in)  :: timeStep
  real,    intent (in)  :: timeSimulation
  integer, intent (in)  :: pulseNumber
  real,    intent (out) :: beamPower

  logical :: noPulse

  integer :: n
  integer :: numberOfSections

  real    :: area
  real    :: energy
  real    :: p1, p2, pa, pb
  real    :: slope
  real    :: t1, t2, ta, tb
  real    :: timeMax, timeMaxPulse
  real    :: timeMin, timeMinPulse
!
!
!     ...Catch bad pulse ID number and set time window limits.
!
!
  if ((pulseNumber < 1) .or. (pulseNumber > ed_numberOfPulses)) then
       call Driver_abortFlash ("ed_computeBeamPower: invalid pulse number!")
  end if

  timeMin = timeSimulation
  timeMax = timeSimulation + timeStep
!
!
!     ...Check, if time window is out of pulse time boundaries. Return with zero
!        beam power in this case.
!
!
!
!             P |
!               |
!               |                /o------o\
!               |               / |      | \
!               |              /  |      |  \
!               |             /   |      |   \     /o
!               |            o    |      |    \   / |
!               |            |    |      |     \ /  |
!               |            |    |      |      o   |
!               |            |    |      |      |   |
!                ------------------------------------------------
!                           t1   t2     t3     t4  t5           T
!                  ^      ^                            ^      ^
!                  |      |                            |      |
!                  |      |                            |      |
!                  |   timeMax                         |   timeMax
!               timeMin                             timeMin
!
!
!
  numberOfSections = ed_pulseNumberOfSections    (pulseNumber)
  timeMinPulse     = ed_pulses (1               , pulseNumber) % pulseTime
  timeMaxPulse     = ed_pulses (numberOfSections, pulseNumber) % pulseTime

  noPulse = (timeMax <= timeMinPulse) .or. (timeMin >= timeMaxPulse)

  if (noPulse) then
      beamPower = 0.0
      return
  end if
!
!
!     ...Some of the pulse shape trapezoids overlap with the time window.
!        Loop over all trapezoids, get the correct 'integration' limits,
!        find the relevant area of the trapezoid and add it to the beam
!        energy.
!
!
  energy = 0.0

  do n = 2, numberOfSections

     t1 = ed_pulses (n-1,pulseNumber) % pulseTime
     t2 = ed_pulses (n,  pulseNumber) % pulseTime

     if (t1 >= timeMax) then
         exit
     end if

     if ((t1 < timeMax) .and. (t2 > timeMin)) then

         ta = max (t1,timeMin)
         tb = min (t2,timeMax)

         p1 = ed_pulses (n-1,pulseNumber) % pulsePower
         p2 = ed_pulses (n  ,pulseNumber) % pulsePower

         slope = (p2 - p1) / (t2 - t1)

         pa = p1 + slope * (ta - t1)
         pb = p2 - slope * (t2 - tb)

         area = 0.5 * (pa + pb) * (tb - ta)

         energy = energy + area

     end if

  end do
!
!
!     ...Calculate the average beam power.
!
!
  beamPower = energy / timeStep
!
!
!     ...Ready!
!
!
  return
end subroutine ed_computeBeamPower
