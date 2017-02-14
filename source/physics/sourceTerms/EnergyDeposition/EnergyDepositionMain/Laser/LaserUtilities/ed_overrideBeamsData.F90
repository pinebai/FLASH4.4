!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_overrideBeamsData
!!
!! NAME
!!
!!  ed_overrideBeamsData
!!
!! SYNOPSIS
!!
!!  use ed_overrideBeamsData
!!
!! DESCRIPTION
!!
!!  This module allows the external user to override data of all the beams.
!!  The character keyword necessary to identify the type of data to be overridden
!!  must be the same as the name of the data variable under which the beam type data
!!  is stored. The main function is overloaded with specific functions for each
!!  data type. If any of the beam number ID or character keyword does not match
!!  with the beams that have been set up, a message is printed and the calculation
!!  aborted.
!!
!! ARGUMENTS
!!
!!  beamID     : the ID number of the beam
!!  entryField : a character string identifying the data
!!  dataValue  : the overriding value of the data
!!
!!***

Module ed_overrideBeamsData

implicit none

interface ed_overrideBeamData
   module procedure ed_overrideBeamDataCharacter
   module procedure ed_overrideBeamDataLogical
   module procedure ed_overrideBeamDataInteger
   module procedure ed_overrideBeamDataReal
end interface

contains
!
!
!     ...The character version.
!
!
subroutine ed_overrideBeamDataCharacter (beamID,     &
                                         entryField, &
                                         dataValue   )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in) :: beamID
  character (len=*), intent (in) :: entryField
  character (len=*), intent (in) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_overrideBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_overrideBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("crossSectionFunctionType")
         ed_beams (beamID) % crossSectionFunctionType = dataValue
  case ("gridType")
         ed_beams (beamID) % gridType                 = dataValue
  case ("semiAxisMajorTorsionAxis")
         ed_beams (beamID) % semiAxisMajorTorsionAxis = dataValue
  case default
         call Driver_abortFlash ("ed_overrideBeamData: No such character entry field!")

  end select

  return
end subroutine ed_overrideBeamDataCharacter
!
!
!     ...The logical version.
!
!
subroutine ed_overrideBeamDataLogical (beamID,     &
                                       entryField, &
                                       dataValue   )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in) :: beamID
  character (len=*), intent (in) :: entryField
  logical,           intent (in) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_overrideBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_overrideBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("ignoreBoundaryCondition")
         ed_beams (beamID) % ignoreBoundaryCondition = dataValue
  case default
         call Driver_abortFlash ("ed_overrideBeamData: No such logical entry field!")

  end select

  return
end subroutine ed_overrideBeamDataLogical
!
!
!     ...The integer version.
!
!
subroutine ed_overrideBeamDataInteger (beamID,     &
                                       entryField, &
                                       dataValue   )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in) :: beamID
  character (len=*), intent (in) :: entryField
  integer,           intent (in) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_overrideBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_overrideBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("dimensionality")
         ed_beams (beamID) % dimensionality   = dataValue
  case ("gridnTics1stDim")
         ed_beams (beamID) % gridnTics1stDim  = dataValue
  case ("gridnTics2ndDim")
         ed_beams (beamID) % gridnTics2ndDim  = dataValue
  case ("gridSeed")
         ed_beams (beamID) % gridSeed         = dataValue
  case ("gridSeedMaximum")
         ed_beams (beamID) % gridSeedMaximum  = dataValue
  case ("gridSeedStepping")
         ed_beams (beamID) % gridSeedStepping = dataValue
  case ("numberOfRays")
         ed_beams (beamID) % numberOfRays     = dataValue
  case ("pulseNumber")
         ed_beams (beamID) % pulseNumber      = dataValue
  case default
         call Driver_abortFlash ("ed_overrideBeamData: No such integer entry field!")

  end select

  return
end subroutine ed_overrideBeamDataInteger
!
!
!     ...The real version.
!
!
subroutine ed_overrideBeamDataReal (beamID,     &
                                    entryField, &
                                    dataValue   )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in) :: beamID
  character (len=*), intent (in) :: entryField
  real,              intent (in) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_overrideBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_overrideBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("distanceLens2Target")
         ed_beams (beamID) % distanceLens2Target       = dataValue
  case ("frequency")
         ed_beams (beamID) % frequency                 = dataValue
  case ("gaussianCenterMajor")
         ed_beams (beamID) % gaussianCenterMajor       = dataValue
  case ("gaussianCenterMinor")
         ed_beams (beamID) % gaussianCenterMinor       = dataValue
  case ("gaussianExponent")
         ed_beams (beamID) % gaussianExponent          = dataValue
  case ("gaussianRadiusMajor")
         ed_beams (beamID) % gaussianRadiusMajor       = dataValue
  case ("gaussianRadiusMinor")
         ed_beams (beamID) % gaussianRadiusMinor       = dataValue
  case ("gridDelta1stDim")
         ed_beams (beamID) % gridDelta1stDim           = dataValue
  case ("gridDelta2ndDim")
         ed_beams (beamID) % gridDelta2ndDim           = dataValue
  case ("gridFirstTic1stDim")
         ed_beams (beamID) %gridFirstTic1stDim         = dataValue
  case ("gridFirstTic2ndDim")
         ed_beams (beamID) %gridFirstTic2ndDim         = dataValue
  case ("gridWeight")
         ed_beams (beamID) % gridWeight                = dataValue
  case ("initialRaySpeed")
         ed_beams (beamID) % initialRaySpeed           = dataValue
  case ("lensSemiAxisMajor")
         ed_beams (beamID) % lensSemiAxisMajor         = dataValue
  case ("lensSemiAxisMinor")
         ed_beams (beamID) % lensSemiAxisMinor         = dataValue
  case ("lensX")
         ed_beams (beamID) % lensX                     = dataValue
  case ("lensY")
         ed_beams (beamID) % lensY                     = dataValue
  case ("lensZ")
         ed_beams (beamID) % lensZ                     = dataValue
  case ("pulseStartingTime")
         ed_beams (beamID) % pulseStartingTime         = dataValue
  case ("pulseEndingTime")
         ed_beams (beamID) % pulseEndingTime           = dataValue
  case ("semiAxisMajorTorsionAngle")
         ed_beams (beamID) % semiAxisMajorTorsionAngle = dataValue
  case ("semiAxisUnitMajorX")
         ed_beams (beamID) % semiAxisUnitMajorX        = dataValue
  case ("semiAxisUnitMajorY")
         ed_beams (beamID) % semiAxisUnitMajorY        = dataValue
  case ("semiAxisUnitMajorZ")
         ed_beams (beamID) % semiAxisUnitMajorZ        = dataValue
  case ("semiAxisUnitMinorX")
         ed_beams (beamID) % semiAxisUnitMinorX        = dataValue
  case ("semiAxisUnitMinorY")
         ed_beams (beamID) % semiAxisUnitMinorY        = dataValue
  case ("semiAxisUnitMinorZ")
         ed_beams (beamID) % semiAxisUnitMinorZ        = dataValue
  case ("target2LensMagnification")
         ed_beams (beamID) % target2LensMagnification  = dataValue
  case ("targetSemiAxisMajor")
         ed_beams (beamID) % targetSemiAxisMajor       = dataValue
  case ("targetSemiAxisMinor")
         ed_beams (beamID) % targetSemiAxisMinor       = dataValue
  case ("targetX")
         ed_beams (beamID) % targetX                   = dataValue
  case ("targetY")
         ed_beams (beamID) % targetY                   = dataValue
  case ("targetZ")
         ed_beams (beamID) % targetZ                   = dataValue
  case ("wavelength")
         ed_beams (beamID) % wavelength                = dataValue
  case default
         call Driver_abortFlash ("ed_overrideBeamData: No such real entry field!")

  end select

  return
end subroutine ed_overrideBeamDataReal

end Module ed_overrideBeamsData
