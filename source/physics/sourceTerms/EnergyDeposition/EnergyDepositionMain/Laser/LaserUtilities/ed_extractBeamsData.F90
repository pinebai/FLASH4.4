!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_extractBeamsData
!!
!! NAME
!!
!!  ed_extractBeamsData
!!
!! SYNOPSIS
!!
!!  use ed_extractBeamsData
!!
!! DESCRIPTION
!!
!!  This module allows the external user to extract all data from all the beams.
!!  The character keyword necessary to identify the type of data wanted must
!!  be the same as the name of the data variable under which the beam type data
!!  is stored. The main function is overloaded with specific functions for each
!!  data type. If any of the beam number ID or character keyword does not match
!!  with the beams that have been set up, a message is printed and the calculation
!!  aborted.
!!
!! ARGUMENTS
!!
!!  beamID     : the ID number of the beam
!!  entryField : a character string identifying the data
!!  dataValue  : the returned value of the data requested
!!
!!***

Module ed_extractBeamsData

implicit none

interface ed_extractBeamData
   module procedure ed_extractBeamDataCharacter
   module procedure ed_extractBeamDataLogical
   module procedure ed_extractBeamDataInteger
   module procedure ed_extractBeamDataReal
end interface

contains
!
!
!     ...The character version.
!
!
subroutine ed_extractBeamDataCharacter (beamID,               &
                                        entryField,           &
                                                    dataValue )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in)  :: beamID
  character (len=*), intent (in)  :: entryField
  character (len=*), intent (out) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_extractBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_extractBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("crossSectionFunctionType")
         dataValue = ed_beams (beamID) % crossSectionFunctionType
  case ("gridType")
         dataValue = ed_beams (beamID) % gridType
  case ("semiAxisMajorTorsionAxis")
         dataValue = ed_beams (beamID) % semiAxisMajorTorsionAxis
  case default
         call Driver_abortFlash ("ed_extractBeamData: No such character entry field!")

  end select

  return
end subroutine ed_extractBeamDataCharacter
!
!
!     ...The logical version.
!
!
subroutine ed_extractBeamDataLogical (beamID,               &
                                      entryField,           &
                                                  dataValue )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in)  :: beamID
  character (len=*), intent (in)  :: entryField
  logical,           intent (out) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_extractBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_extractBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("ignoreBoundaryCondition")
         dataValue = ed_beams (beamID) % ignoreBoundaryCondition
  case default
         call Driver_abortFlash ("ed_extractBeamData: No such logical entry field!")

  end select

  return
end subroutine ed_extractBeamDataLogical
!
!
!     ...The integer version.
!
!
subroutine ed_extractBeamDataInteger (beamID,               &
                                      entryField,           &
                                                  dataValue )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in)  :: beamID
  character (len=*), intent (in)  :: entryField
  integer,           intent (out) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_extractBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_extractBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("dimensionality")
         dataValue = ed_beams (beamID) % dimensionality
  case ("gridnTics1stDim")
         dataValue = ed_beams (beamID) % gridnTics1stDim
  case ("gridnTics2ndDim")
         dataValue = ed_beams (beamID) % gridnTics2ndDim
  case ("gridSeed")
         dataValue = ed_beams (beamID) % gridSeed
  case ("gridSeedMaximum")
         dataValue = ed_beams (beamID) % gridSeedMaximum
  case ("gridSeedStepping")
         dataValue = ed_beams (beamID) % gridSeedStepping
  case ("numberOfRays")
         dataValue = ed_beams (beamID) % numberOfRays
  case ("pulseNumber")
         dataValue = ed_beams (beamID) % pulseNumber
  case default
         call Driver_abortFlash ("ed_extractBeamData: No such integer entry field!")

  end select

  return
end subroutine ed_extractBeamDataInteger
!
!
!     ...The real version.
!
!
subroutine ed_extractBeamDataReal (beamID,               &
                                   entryField,           &
                                               dataValue )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_beams,         &
                                     ed_beamsAreSetup, &
                                     ed_numberOfBeams

  implicit none

  integer,           intent (in)  :: beamID
  character (len=*), intent (in)  :: entryField
  real,              intent (out) :: dataValue

  if (.not.ed_beamsAreSetup) then
       call Driver_abortFlash ("ed_extractBeamData: No beams are set up!")
  end if

  if ((beamID < 1) .or. (beamID > ed_numberOfBeams) ) then
       call Driver_abortFlash ("ed_extractBeamData: Beam ID out of range!")
  end if

  select case (entryField)

  case ("distanceLens2Target")
         dataValue = ed_beams (beamID) % distanceLens2Target
  case ("frequency")
         dataValue = ed_beams (beamID) % frequency
  case ("gaussianCenterMajor")
         dataValue = ed_beams (beamID) % gaussianCenterMajor
  case ("gaussianCenterMinor")
         dataValue = ed_beams (beamID) % gaussianCenterMinor
  case ("gaussianExponent")
         dataValue = ed_beams (beamID) % gaussianExponent
  case ("gaussianRadiusMajor")
         dataValue = ed_beams (beamID) % gaussianRadiusMajor
  case ("gaussianRadiusMinor")
         dataValue = ed_beams (beamID) % gaussianRadiusMinor
  case ("gridDelta1stDim")
         dataValue = ed_beams (beamID) % gridDelta1stDim
  case ("gridDelta2ndDim")
         dataValue = ed_beams (beamID) % gridDelta2ndDim
  case ("gridFirstTic1stDim")
         dataValue = ed_beams (beamID) % gridFirstTic1stDim
  case ("gridFirstTic2ndDim")
         dataValue = ed_beams (beamID) % gridFirstTic2ndDim
  case ("gridWeight")
         dataValue = ed_beams (beamID) % gridWeight
  case ("initialRaySpeed")
         dataValue = ed_beams (beamID) % initialRaySpeed
  case ("lensSemiAxisMajor")
         dataValue = ed_beams (beamID) % lensSemiAxisMajor
  case ("lensSemiAxisMinor")
         dataValue = ed_beams (beamID) % lensSemiAxisMinor
  case ("lensX")
         dataValue = ed_beams (beamID) % lensX
  case ("lensY")
         dataValue = ed_beams (beamID) % lensY
  case ("lensZ")
         dataValue = ed_beams (beamID) % lensZ
  case ("pulseStartingTime")
         dataValue = ed_beams (beamID) % pulseStartingTime
  case ("pulseEndingTime")
         dataValue = ed_beams (beamID) % pulseEndingTime
  case ("semiAxisMajorTorsionAngle")
         dataValue = ed_beams (beamID) % semiAxisMajorTorsionAngle
  case ("semiAxisUnitMajorX")
         dataValue = ed_beams (beamID) % semiAxisUnitMajorX
  case ("semiAxisUnitMajorY")
         dataValue = ed_beams (beamID) % semiAxisUnitMajorY
  case ("semiAxisUnitMajorZ")
         dataValue = ed_beams (beamID) % semiAxisUnitMajorZ
  case ("semiAxisUnitMinorX")
         dataValue = ed_beams (beamID) % semiAxisUnitMinorX
  case ("semiAxisUnitMinorY")
         dataValue = ed_beams (beamID) % semiAxisUnitMinorY
  case ("semiAxisUnitMinorZ")
         dataValue = ed_beams (beamID) % semiAxisUnitMinorZ
  case ("target2LensMagnification")
         dataValue = ed_beams (beamID) % target2LensMagnification
  case ("targetSemiAxisMajor")
         dataValue = ed_beams (beamID) % targetSemiAxisMajor
  case ("targetSemiAxisMinor")
         dataValue = ed_beams (beamID) % targetSemiAxisMinor
  case ("targetX")
         dataValue = ed_beams (beamID) % targetX
  case ("targetY")
         dataValue = ed_beams (beamID) % targetY
  case ("targetZ")
         dataValue = ed_beams (beamID) % targetZ
  case ("wavelength")
         dataValue = ed_beams (beamID) % wavelength
  case default
         call Driver_abortFlash ("ed_extractBeamData: No such real entry field!")

  end select

  return
end subroutine ed_extractBeamDataReal

end Module ed_extractBeamsData
