!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_setupBeams
!!
!! NAME
!!
!!  ed_setupBeams
!!
!! SYNOPSIS
!!
!!  call ed_setupBeams ()
!!
!! DESCRIPTION
!!
!!  Sets up and characterizes the beams to be used in the simulation.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  All needed info is read in as runtime parameters.
!!
!!***

subroutine ed_setupBeams ()

  use Driver_interface,            ONLY : Driver_abortFlash

  use ed_interface,                ONLY : ed_beamsCheck,      &
                                          ed_beamsInfo,       &
                                          ed_initializeBeams
  
  use EnergyDeposition_data,       ONLY : ed_beams,            &
                                          ed_beamsAreSetup,    &
                                          ed_degrees2rad,      &
                                          ed_gridGeometry,     &
                                          ed_laser3Din2D,      &
                                          ed_microns2cm,       &
                                          ed_numberOfBeams,    &
                                          ed_numberOfPulses,   &
                                          ed_pulsesAreSetup,   &
                                          ed_speedOfLight

  use Logfile_interface,           ONLY : Logfile_stampMessage

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  character :: semiAxisMajorTorsionAxis
  
  character (len = BEAM_STRING_LENGTH) :: crossSectionFunctionType
  character (len = BEAM_STRING_LENGTH) :: gridType
  character (len =  MAX_STRING_LENGTH) :: parameterString

  logical :: ambiguousGeometry
  logical :: beams1D, beams2D, beams3D
  logical :: ignoreBoundaryCondition

  integer :: beam
  integer :: gridnAngularTics
  integer :: gridnRadialTics
  integer :: gridnSemiAxisMajorTics
  integer :: gridnSemiAxisMinorTics
  integer :: numberOfRays
  integer :: pulseNumber

  real    :: gaussianExponent
  real    :: gaussianRadiusMajor
  real    :: gaussianRadiusMinor
  real    :: gaussianCenterMajor
  real    :: gaussianCenterMinor
  real    :: gridDeltaSemiAxisMajor
  real    :: gridDeltaSemiAxisMinor
  real    :: initialRaySpeed
  real    :: lensSemiAxisMajor
  real    :: lensX
  real    :: lensY
  real    :: lensZ
  real    :: semiAxisMajorTorsionAngle
  real    :: targetSemiAxisMajor
  real    :: targetSemiAxisMinor
  real    :: targetX
  real    :: targetY
  real    :: targetZ
  real    :: wavelength
!
!
!     ...Check, if the pulse(s) have been set up. If not, we cannot complete the beam(s)
!        data and we must stop.
!
!
  if (.not.ed_pulsesAreSetup) then
      call Driver_abortFlash ("ed_setupBeams: Pulses are not set up!")
  end if

  if (ed_numberOfBeams < 1) then
      call Driver_abortFlash ("ed_setupBeams: No beam defined!")
  end if

  if (ed_numberOfBeams > ED_MAXBEAMS) then
      call Logfile_stampMessage ("ed_setupBeams: ERROR")
      call Logfile_stampMessage ("# of beams > maximum # of beam runtime parameters!")
      call Logfile_stampMessage ("Not enough beam runtime parameters created.")
      call Logfile_stampMessage ("Increase the maximum number of beam runtime parameters")
      call Logfile_stampMessage ("using the ed_maxBeams=<number of beams> setup option.")
      call Driver_abortFlash    ("ed_setupBeams: Not enough beam runtime parameters. See Log File!")
  end if
!
!
!     ...Determine the geometrical nature of the beams.
!
!  
  beams1D =     (ed_gridGeometry == GRID_1DCARTESIAN)   &
           .or. (ed_gridGeometry == GRID_1DSPHERICAL)   &
           .or. (ed_gridGeometry == GRID_1DPOLAR)       &
           .or. (ed_gridGeometry == GRID_1DCYLINDRICAL  )

  beams2D =     (ed_gridGeometry == GRID_2DCARTESIAN)   &
           .or. (ed_gridGeometry == GRID_2DSPHERICAL)   &
           .or. (ed_gridGeometry == GRID_2DPOLAR)       &
           .or. (ed_gridGeometry == GRID_2DCYLINDRICAL .and. .not.ed_laser3Din2D)

  beams3D =     (ed_gridGeometry == GRID_3DCARTESIAN)   &
           .or. (ed_gridGeometry == GRID_3DCYLINDRICAL) &
           .or. (ed_gridGeometry == GRID_3DSPHERICAL)   &
           .or. (ed_gridGeometry == GRID_3DPOLAR)       &
           .or. (ed_gridGeometry == GRID_2DCYLINDRICAL .and. ed_laser3Din2D)

  ambiguousGeometry =     (beams1D .and. beams2D) &
                     .or. (beams1D .and. beams3D) &
                     .or. (beams2D .and. beams3D)

  if (ambiguousGeometry) then
      call Driver_abortFlash ("ed_setupBeams: Ambiguous geometrical nature of beams!")
  end if
!
!
!     ...Allocate the beams array and initialize all beams.
!
!  
  allocate (ed_beams (1:ed_numberOfBeams))
  
  call ed_initializeBeams ()
!
!
!     ...Read the beam runtime parameters.
!
!  
  do beam = 1, ed_numberOfBeams

     write (parameterString,'(a,i0)') "ed_lensX_", beam
     call RuntimeParameters_get (parameterString, lensX)

     write (parameterString,'(a,i0)') "ed_lensY_", beam
     call RuntimeParameters_get (parameterString, lensY)

     write (parameterString,'(a,i0)') "ed_lensZ_", beam
     call RuntimeParameters_get (parameterString, lensZ)

     write (parameterString,'(a,i0)') "ed_targetX_", beam
     call RuntimeParameters_get (parameterString, targetX)

     write (parameterString,'(a,i0)') "ed_targetY_", beam
     call RuntimeParameters_get (parameterString, targetY)

     write (parameterString,'(a,i0)') "ed_targetZ_", beam
     call RuntimeParameters_get (parameterString, targetZ)

     write (parameterString,'(a,i0)') "ed_targetSemiAxisMajor_", beam
     call RuntimeParameters_get (parameterString, targetSemiAxisMajor)

     write (parameterString,'(a,i0)') "ed_targetSemiAxisMinor_", beam
     call RuntimeParameters_get (parameterString, targetSemiAxisMinor)

     write (parameterString,'(a,i0)') "ed_lensSemiAxisMajor_", beam
     call RuntimeParameters_get (parameterString, lensSemiAxisMajor)

     write (parameterString,'(a,i0)') "ed_semiAxisMajorTorsionAngle_", beam
     call RuntimeParameters_get (parameterString, semiAxisMajorTorsionAngle)

     write (parameterString,'(a,i0)') "ed_semiAxisMajorTorsionAxis_", beam
     call RuntimeParameters_get (parameterString, semiAxisMajorTorsionAxis)

     write (parameterString,'(a,i0)') "ed_pulseNumber_", beam
     call RuntimeParameters_get (parameterString, pulseNumber)

     write (parameterString,'(a,i0)') "ed_wavelength_", beam
     call RuntimeParameters_get (parameterString, wavelength)

     write (parameterString,'(a,i0)') "ed_initialRaySpeed_", beam
     call RuntimeParameters_get (parameterString, initialRaySpeed)

     write (parameterString,'(a,i0)') "ed_ignoreBoundaryCondition_", beam
     call RuntimeParameters_get (parameterString, ignoreBoundaryCondition)

     write (parameterString,'(a,i0)') "ed_crossSectionFunctionType_", beam
     call RuntimeParameters_get (parameterString, crossSectionFunctionType)

     write (parameterString,'(a,i0)') "ed_gaussianExponent_", beam
     call RuntimeParameters_get (parameterString, gaussianExponent)

     write (parameterString,'(a,i0)') "ed_gaussianRadiusMajor_", beam
     call RuntimeParameters_get (parameterString, gaussianRadiusMajor)

     write (parameterString,'(a,i0)') "ed_gaussianRadiusMinor_", beam
     call RuntimeParameters_get (parameterString, gaussianRadiusMinor)

     write (parameterString,'(a,i0)') "ed_gaussianCenterMajor_", beam
     call RuntimeParameters_get (parameterString, gaussianCenterMajor)

     write (parameterString,'(a,i0)') "ed_gaussianCenterMinor_", beam
     call RuntimeParameters_get (parameterString, gaussianCenterMinor)

     write (parameterString,'(a,i0)') "ed_numberOfRays_", beam
     call RuntimeParameters_get (parameterString, numberOfRays)

     write (parameterString,'(a,i0)') "ed_gridType_", beam
     call RuntimeParameters_get (parameterString, gridType)

     write (parameterString,'(a,i0)') "ed_gridnRadialTics_", beam
     call RuntimeParameters_get (parameterString, gridnRadialTics)

     write (parameterString,'(a,i0)') "ed_gridnAngularTics_", beam
     call RuntimeParameters_get (parameterString, gridnAngularTics)

     write (parameterString,'(a,i0)') "ed_gridnSemiAxisMajorTics_", beam
     call RuntimeParameters_get (parameterString, gridnSemiAxisMajorTics)

     write (parameterString,'(a,i0)') "ed_gridnSemiAxisMinorTics_", beam
     call RuntimeParameters_get (parameterString, gridnSemiAxisMinorTics)

     write (parameterString,'(a,i0)') "ed_gridDeltaSemiAxisMajor_", beam
     call RuntimeParameters_get (parameterString, gridDeltaSemiAxisMajor)

     write (parameterString,'(a,i0)') "ed_gridDeltaSemiAxisMinor_", beam
     call RuntimeParameters_get (parameterString, gridDeltaSemiAxisMinor)
!
!
!     ...Catch any bad data at this point.
!
!  
     if (numberOfRays < 1) then
         write (*,*) ' Beam # , number of rays = ',beam, numberOfRays
         call Driver_abortFlash ("ed_setupBeams: Bad number of rays specified!")
     end if

     if ((pulseNumber < 1) .or. (pulseNumber > ed_numberOfPulses)) then
         write (*,*) ' Beam # , pulse number = ',beam, pulseNumber
         call Driver_abortFlash ("ed_setupBeams: invalid pulse number!")
     end if
!
!
!     ...Bad 1D beam data (with override, hence no abortion of program).
!
!  
     if (beams1D) then

         if (numberOfRays > 1) then
             write (*,*) ' Beam # , number of rays = ',beam, numberOfRays
             write (*,*) ' Only 1 ray is needed -> # of rays enforced to 1 '
             numberOfRays = 1
         end if

     end if
!
!
!     ...Bad 2D beam data.
!
!  
     if (beams2D) then

         if (targetSemiAxisMajor <= 0.0) then
             write (*,*) ' Beam # , target semiaxis length = ',beam, targetSemiAxisMajor
             call Driver_abortFlash ("ed_setupBeams: target semiaxis length <= 0!")
         end if

         if (lensSemiAxisMajor <= 0.0) then
             write (*,*) ' Beam # , lens semiaxis length = ',beam, lensSemiAxisMajor
             call Driver_abortFlash ("ed_setupBeams: lens semiaxis length <= 0!")
         end if

         if (      (gridType /= 'regular1D'    ) &
             .and. (gridType /= 'statistical1D') ) then
              write (*,*) ' Beam # , grid type = ',beam, gridType
              write (*,*) ' Should be: regular1D or statistical1D '
              call Driver_abortFlash ("ed_setupBeams: invalid 2D beam grid type!")
         end if

         if (      (crossSectionFunctionType /= 'uniform'          ) &
             .and. (crossSectionFunctionType /= 'gaussian1D'       ) &
             .and. (crossSectionFunctionType /= 'gaussianInverse1D') ) then
              write (*,*) ' Beam # , cross section function type = ',beam, crossSectionFunctionType
              write (*,*) ' Should be: uniform, gaussian1D or gaussianInverse1D '
              call Driver_abortFlash ("ed_setupBeams: invalid 2D beam cross section function type!")
         end if

     end if
!
!
!     ...Bad 3D beam data.
!
!  
     if (beams3D) then

         if (targetSemiAxisMajor <= 0.0) then
             write (*,*) ' Beam # , target major semiaxis length = ',beam, targetSemiAxisMajor
             call Driver_abortFlash ("ed_setupBeams: target major semiaxis length <= 0!")
         end if

         if (targetSemiAxisMinor <= 0.0) then
             write (*,*) ' Beam # , target minor semiaxis length = ',beam, targetSemiAxisMinor
             call Driver_abortFlash ("ed_setupBeams: target minor semiaxis length <= 0!")
         end if

         if (targetSemiAxisMajor < targetSemiAxisMinor) then
             call Driver_abortFlash ("ed_setupBeams: target semiaxes lengths out of order (major < minor)!")
         end if

         if (lensSemiAxisMajor <= 0.0) then
             write (*,*) ' Beam # , lens major semiaxis length = ',beam, lensSemiAxisMajor
             call Driver_abortFlash ("ed_setupBeams: lens major semiaxis length <= 0!")
         end if

         if ((semiAxisMajorTorsionAngle < 0.0) .or. (semiAxisMajorTorsionAngle > 360.0)) then
             write (*,*) ' Beam # , major semiaxis torsion angle = ',beam, semiAxisMajorTorsionAngle
             call Driver_abortFlash ("ed_setupBeams: invalid major semiaxis torsion angle!")
         end if

         if (      (gridType /= 'square2D'     ) &
             .and. (gridType /= 'delta2D'      ) &
             .and. (gridType /= 'radial2D'     ) &
             .and. (gridType /= 'rectangular2D') &
             .and. (gridType /= 'statistical2D') ) then
              write (*,*) ' Beam # , grid type = ',beam, gridType
              write (*,*) ' Should be: square2D, radial2D, delta2D, rectangular2D or statistical2D '
              call Driver_abortFlash ("ed_setupBeams: invalid 3D beam grid type!")
         end if

         if (      (semiAxisMajorTorsionAxis /= 'x') &
             .and. (semiAxisMajorTorsionAxis /= 'y') &
             .and. (semiAxisMajorTorsionAxis /= 'z') ) then
              write (*,*) ' Beam # , major semiaxis torsion axis = ',beam, semiAxisMajorTorsionAxis
              call Driver_abortFlash ("ed_setupBeams: invalid major semiaxis torsion axis!")
         end if

         if (      (crossSectionFunctionType /= 'uniform'          ) &
             .and. (crossSectionFunctionType /= 'gaussian2D'       ) &
             .and. (crossSectionFunctionType /= 'gaussianInverse2D') ) then
              write (*,*) ' Beam # , cross section function type = ',beam, crossSectionFunctionType
              write (*,*) ' Should be: uniform, gaussian2D or gaussianInverse2D '
              call Driver_abortFlash ("ed_setupBeams: invalid 3D beam cross section function type!")
         end if

     end if
!
!
!     ...Perform some unit conversions.
!
!
     wavelength                = wavelength                * ed_microns2cm         ! convert to cm
     semiAxisMajorTorsionAngle = semiAxisMajorTorsionAngle * ed_degrees2rad        ! convert to radians
!
!
!     ...Store beam data obtained so far into appropriate places.
!
!  
     ed_beams (beam) % crossSectionFunctionType  = crossSectionFunctionType        ! string
     ed_beams (beam) % gaussianExponent          = gaussianExponent
     ed_beams (beam) % gaussianRadiusMajor       = gaussianRadiusMajor             ! in cm
     ed_beams (beam) % gaussianRadiusMinor       = gaussianRadiusMinor             ! in cm
     ed_beams (beam) % gaussianCenterMajor       = gaussianCenterMajor             ! in cm
     ed_beams (beam) % gaussianCenterMinor       = gaussianCenterMinor             ! in cm
     ed_beams (beam) % frequency                 = ed_speedOfLight / wavelength    ! in Hz (s^-1)
     ed_beams (beam) % gridType                  = gridType                        ! string
     ed_beams (beam) % ignoreBoundaryCondition   = ignoreBoundaryCondition         ! logical
     ed_beams (beam) % initialRaySpeed           = initialRaySpeed                 ! in units of light speed
     ed_beams (beam) % lensSemiAxisMajor         = lensSemiAxisMajor               ! in cm
     ed_beams (beam) % lensX                     = lensX                           ! in cm
     ed_beams (beam) % lensY                     = lensY                           ! in cm
     ed_beams (beam) % lensZ                     = lensZ                           ! in cm
     ed_beams (beam) % numberOfRays              = numberOfRays
     ed_beams (beam) % pulseNumber               = pulseNumber
     ed_beams (beam) % semiAxisMajorTorsionAngle = semiAxisMajorTorsionAngle       ! in radians
     ed_beams (beam) % semiAxisMajorTorsionAxis  = semiAxisMajorTorsionAxis        ! one character
     ed_beams (beam) % targetSemiAxisMajor       = targetSemiAxisMajor             ! in cm
     ed_beams (beam) % targetSemiAxisMinor       = targetSemiAxisMinor             ! in cm
     ed_beams (beam) % targetX                   = targetX                         ! in cm
     ed_beams (beam) % targetY                   = targetY                         ! in cm
     ed_beams (beam) % targetZ                   = targetZ                         ! in cm
     ed_beams (beam) % wavelength                = wavelength                      ! in cm

     if (gridType == 'delta2D') then
         ed_beams (beam) % gridDelta1stDim = gridDeltaSemiAxisMajor
         ed_beams (beam) % gridDelta2ndDim = gridDeltaSemiAxisMinor
     else if (gridType == 'radial2D') then
         ed_beams (beam) % gridnTics1stDim = gridnRadialTics
         ed_beams (beam) % gridnTics2ndDim = gridnAngularTics
     else if (gridType == 'rectangular2D') then
         ed_beams (beam) % gridnTics1stDim = gridnSemiAxisMajorTics
         ed_beams (beam) % gridnTics2ndDim = gridnSemiAxisMinorTics
     end if

  enddo
!
!
!     ...Extract more beam info (if needed) and check the gathered beams data.
!        These routines call geometry specific subroutines.
!
!  
  call ed_beamsInfo  ()
  call ed_beamsCheck ()
!
!
!     ...Set beams status indicator.
!
!
  ed_beamsAreSetup = .true.
!
!
!     ...Ready!
!
!
  return
end subroutine ed_setupBeams
