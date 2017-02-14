!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonBeams/pi_setupBeams
!!
!! NAME
!!
!!  pi_setupBeams
!!
!! SYNOPSIS
!!
!!  call pi_setupBeams ()
!!
!! DESCRIPTION
!!
!!  Sets up and characterizes the beams to be used for proton imaging.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  All needed info is read in as runtime parameters.
!!
!!***

subroutine pi_setupBeams ()

  use Driver_interface,            ONLY : Driver_abortFlash

  use pi_interface,                ONLY : pi_beamsCheck, &
                                          pi_beamsInfo
  
  use ProtonImaging_data,          ONLY : pi_beams,             &
                                          pi_beamsAreSetup,     &
                                          pi_degrees2rad,       &
                                          pi_gridGeometry,      &
                                          pi_3Din2D,            &
                                          pi_microns2cm,        &
                                          pi_numberOfBeams,     &
                                          pi_numberOfDetectors, &
                                          pi_speedOfLight

  use Logfile_interface,           ONLY : Logfile_stampMessage

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "ProtonImaging.h"
#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: parameterString

  logical :: beams1D, beams2D, beams3D
  logical :: noBoundaryCondition

  integer :: beam
  integer :: capsuleGrainLevel
  integer :: detector
  integer :: numberOfProtons

  real    :: apertureAngle
  real    :: capsuleRadius
  real    :: capsuleX
  real    :: capsuleY
  real    :: capsuleZ
  real    :: protonEnergy
  real    :: targetRadius
  real    :: targetX
  real    :: targetY
  real    :: targetZ
  real    :: time2Launch
!
!
!     ...Determine the geometrical nature of the beams and catch bad choices.
!
!  
  beams1D =     (pi_gridGeometry == GRID_1DCARTESIAN)   &
           .or. (pi_gridGeometry == GRID_1DSPHERICAL)   &
           .or. (pi_gridGeometry == GRID_1DPOLAR)       &
           .or. (pi_gridGeometry == GRID_1DCYLINDRICAL  )

  beams2D =     (pi_gridGeometry == GRID_2DCARTESIAN)   &
           .or. (pi_gridGeometry == GRID_2DSPHERICAL)   &
           .or. (pi_gridGeometry == GRID_2DPOLAR)       &
           .or. (pi_gridGeometry == GRID_2DCYLINDRICAL .and. .not.pi_3Din2D)

  beams3D =     (pi_gridGeometry == GRID_3DCARTESIAN)   &
           .or. (pi_gridGeometry == GRID_3DCYLINDRICAL) &
           .or. (pi_gridGeometry == GRID_3DSPHERICAL)   &
           .or. (pi_gridGeometry == GRID_3DPOLAR)       &
           .or. (pi_gridGeometry == GRID_2DCYLINDRICAL .and. pi_3Din2D)

  if (beams1D .or. beams2D) then
      call Driver_abortFlash ("pi_setupBeams: 1D or 2D geometry! Proton imaging makes no sense!")
  end if
!
!
!     ...Allocate the beams array and read the beam runtime parameters.
!
!  
  allocate (pi_beams (1:pi_numberOfBeams))

  do beam = 1, pi_numberOfBeams

     write (parameterString,'(a,i0)') "pi_beamCapsuleX_", beam
     call RuntimeParameters_get (parameterString, capsuleX)

     write (parameterString,'(a,i0)') "pi_beamCapsuleY_", beam
     call RuntimeParameters_get (parameterString, capsuleY)

     write (parameterString,'(a,i0)') "pi_beamCapsuleZ_", beam
     call RuntimeParameters_get (parameterString, capsuleZ)

     write (parameterString,'(a,i0)') "pi_beamCapsuleRadius_", beam
     call RuntimeParameters_get (parameterString, capsuleRadius)

     write (parameterString,'(a,i0)') "pi_beamCapsuleGrainLevel_", beam
     call RuntimeParameters_get (parameterString, capsuleGrainLevel)

     write (parameterString,'(a,i0)') "pi_beamTargetX_", beam
     call RuntimeParameters_get (parameterString, targetX)

     write (parameterString,'(a,i0)') "pi_beamTargetY_", beam
     call RuntimeParameters_get (parameterString, targetY)

     write (parameterString,'(a,i0)') "pi_beamTargetZ_", beam
     call RuntimeParameters_get (parameterString, targetZ)

     write (parameterString,'(a,i0)') "pi_beamTargetRadius_", beam
     call RuntimeParameters_get (parameterString, targetRadius)

     write (parameterString,'(a,i0)') "pi_beamApertureAngle_", beam
     call RuntimeParameters_get (parameterString, apertureAngle)

     write (parameterString,'(a,i0)') "pi_beamProtonEnergy_", beam
     call RuntimeParameters_get (parameterString, protonEnergy)

     write (parameterString,'(a,i0)') "pi_beamTime2Launch_", beam
     call RuntimeParameters_get (parameterString, time2Launch)

     write (parameterString,'(a,i0)') "pi_beamNumberOfProtons_", beam
     call RuntimeParameters_get (parameterString, numberOfProtons)

     write (parameterString,'(a,i0)') "pi_beamDetector_", beam
     call RuntimeParameters_get (parameterString, detector)

     write (parameterString,'(a,i0)') "pi_beamNoBoundaryCondition_", beam
     call RuntimeParameters_get (parameterString, noBoundaryCondition)
!
!
!     ...Catch any bad data at this point.
!
!  
     if (capsuleRadius <= 0.0) then
         write (*,*) ' Beam # , capsule radius = ',beam, capsuleRadius
         call Driver_abortFlash ("pi_setupBeams: Capsule radius of beam must be > 0!")
     end if

     if (capsuleGrainLevel < 0) then
         write (*,*) ' Beam # , capsule grain level = ',beam, capsuleGrainLevel
         call Driver_abortFlash ("pi_setupBeams: Capsule grain level of beam must be >= 0!")
     end if

     if (numberOfProtons < 1) then
         write (*,*) ' Beam # , number of protons = ',beam, numberOfProtons
         call Driver_abortFlash ("pi_setupBeams: Bad number of protons specified!")
     end if

     if (protonEnergy <= 0.0) then
         write (*,*) ' Beam # , proton energy (MeV) = ',beam, protonEnergy
         call Driver_abortFlash ("pi_setupBeams: The proton energy must be > 0!")
     end if

     if (detector < 1 .or. detector > pi_numberOfDetectors) then
         write (*,*) ' Beam # , target detector # = ',beam, detector
         call Driver_abortFlash ("pi_setupBeams: No valid detector specified for beam!")
     end if

     if (apertureAngle < 0.0 .or. apertureAngle >= 180.0) then
         write (*,*) ' Beam # , aperture angle = ',beam, apertureAngle
         call Driver_abortFlash ("pi_setupBeams: Aperture angle must be between 0 and < 180 deg!")
     end if
!
!
!     ...Perform some unit conversions.
!
!
     apertureAngle = apertureAngle * pi_degrees2rad        ! convert to radians
!
!
!     ...Store beam data obtained so far into appropriate places.
!
!  
     pi_beams (beam) % apertureAngle             = apertureAngle                   ! in radians
     pi_beams (beam) % capsuleGrainLevel         = capsuleGrainLevel
     pi_beams (beam) % capsuleRadius             = capsuleRadius                   ! in cm
     pi_beams (beam) % capsuleX                  = capsuleX                        ! in cm
     pi_beams (beam) % capsuleY                  = capsuleY                        ! in cm
     pi_beams (beam) % capsuleZ                  = capsuleZ                        ! in cm
     pi_beams (beam) % detector                  = detector
     pi_beams (beam) % noBoundaryCondition       = noBoundaryCondition             ! logical
     pi_beams (beam) % numberOfProtons           = numberOfProtons
     pi_beams (beam) % protonEnergy              = protonEnergy                    ! in MeV
     pi_beams (beam) % targetRadius              = targetRadius                    ! in cm
     pi_beams (beam) % targetX                   = targetX                         ! in cm
     pi_beams (beam) % targetY                   = targetY                         ! in cm
     pi_beams (beam) % targetZ                   = targetZ                         ! in cm
     pi_beams (beam) % time2Launch               = time2Launch                     ! in sec

  enddo
!
!
!     ...Extract more beam info (if needed) and check the gathered beams data.
!        These routines call geometry specific subroutines.
!
!  
  call pi_beamsInfo  ()
  call pi_beamsCheck ()
!
!
!     ...Set beams status indicator.
!
!
  pi_beamsAreSetup = .true.
!
!
!     ...Ready!
!
!
  return
end subroutine pi_setupBeams
