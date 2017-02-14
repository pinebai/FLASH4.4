!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_printResults
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
!!   This routine prints the obtained ray exit coordinates and power (together with their
!!   percentage errors) to a .dat file. Only the master processor can make the printout, as it
!!   is the only processor having all ray exit data.
!!
!!***

subroutine sim_printResults ()

  use Simulation_data,             ONLY : sim_baseName,               &
                                          sim_geometry,               &
                                          sim_globalMe,               &
                                          sim_powerDecayFactor,       &
                                          sim_rayPexit,               &
                                          sim_rayXexit,               &
                                          sim_rayYexit,               &
                                          sim_rayZexit,               &
                                          sim_rayFexitPercentError1,  &
                                          sim_rayFexitPercentError12, &
                                          sim_rayPexitPercentError1,  &
                                          sim_rayPexitPercentError12, &
                                          sim_xw,                     &
                                          sim_yw,                     &
                                          sim_zw,                     &
                                          sim_lasersOrientation,      &
                                          GRID_2DCARTESIAN,           &
                                          GRID_2DCYLINDRICAL,         &
                                          GRID_3DCARTESIAN

  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer  :: fileUnit
  integer  :: ray
  integer  :: ut_getFreeFileUnit

  real     :: rayP, rayPanalytic, rayPerror, rayPerrorBar
  real     :: rayX, rayXanalytic, rayXerror, rayXerrorBar
  real     :: rayY, rayYanalytic, rayYerror, rayYerrorBar
  real     :: rayZ, rayZanalytic, rayZerror, rayZerrorBar
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
!     ...Print out the rays exit x-coordinates (if any).
!
!
  if (     sim_lasersOrientation == 'X'  &
      .or. sim_lasersOrientation == 'XY' &
      .or. sim_lasersOrientation == 'XZ' ) then

      write (fileUnit,'(/)')
      write (fileUnit,'(a)') "                         RAY X EXIT (errors in %)"
      write (fileUnit,'(/)')
      write (fileUnit,'(a)') " Ray         X exit         X exit analytic     X error    X error bar"
      write (fileUnit,'(a)') " ---------------------------------------------------------------------"

      if (sim_geometry == GRID_2DCARTESIAN .or. sim_geometry == GRID_2DCYLINDRICAL) then
          do ray = 1,2
             rayX         = sim_rayXexit (ray)
             rayXanalytic = sim_xw
             rayXerror    = abs (rayX - sim_xw) * 100.0 / sim_xw
             rayXerrorBar = sim_rayFexitPercentError1
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayX, rayXanalytic, rayXerror, rayXerrorBar
          end do
      else if (sim_geometry == GRID_3DCARTESIAN) then
          do ray = 1,4
             rayX         = sim_rayXexit (ray)
             rayXanalytic = sim_xw
             rayXerror    = abs (rayX - sim_xw) * 100.0 / sim_xw
             rayXerrorBar = sim_rayFexitPercentError1
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayX, rayXanalytic, rayXerror, rayXerrorBar
          end do
          do ray = 5,8
             rayX         = sim_rayXexit (ray)
             rayXanalytic = sim_xw
             rayXerror    = abs (rayX - sim_xw) * 100.0 / sim_xw
             rayXerrorBar = sim_rayFexitPercentError12
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayX, rayXanalytic, rayXerror, rayXerrorBar
          end do
      end if

  end if
!
!
!     ...Print out the rays exit y-coordinates (if any).
!
!
  if (     sim_lasersOrientation == 'Y'  &
      .or. sim_lasersOrientation == 'XY' &
      .or. sim_lasersOrientation == 'YZ' ) then

      write (fileUnit,'(/)')
      write (fileUnit,'(a)') "                         RAY Y EXIT (errors in %)"
      write (fileUnit,'(/)')
      write (fileUnit,'(a)') " Ray         Y exit         Y exit analytic     Y error    Y error bar"
      write (fileUnit,'(a)') " ---------------------------------------------------------------------"

      if (sim_geometry == GRID_2DCARTESIAN .or. sim_geometry == GRID_2DCYLINDRICAL) then
          do ray = 1,2
             rayY         = sim_rayYexit (ray)
             rayYanalytic = sim_yw
             rayYerror    = abs (rayY - sim_yw) * 100.0 / sim_yw
             rayYerrorBar = sim_rayFexitPercentError1
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayY, rayYanalytic, rayYerror, rayYerrorBar
          end do
      else if (sim_geometry == GRID_3DCARTESIAN) then
          do ray = 1,4
             rayY         = sim_rayYexit (ray)
             rayYanalytic = sim_yw
             rayYerror    = abs (rayY - sim_yw) * 100.0 / sim_yw
             rayYerrorBar = sim_rayFexitPercentError1
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayY, rayYanalytic, rayYerror, rayYerrorBar
          end do
          do ray = 5,8
             rayY         = sim_rayYexit (ray)
             rayYanalytic = sim_yw
             rayYerror    = abs (rayY - sim_yw) * 100.0 / sim_yw
             rayYerrorBar = sim_rayFexitPercentError12
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayY, rayYanalytic, rayYerror, rayYerrorBar
          end do
      end if

  end if
!
!
!     ...Print out the rays exit z-coordinates (if any).
!
!
  if (     sim_lasersOrientation == 'XZ' &
      .or. sim_lasersOrientation == 'YZ' ) then

      write (fileUnit,'(/)')
      write (fileUnit,'(a)') "                         RAY Z EXIT (errors in %)"
      write (fileUnit,'(/)')
      write (fileUnit,'(a)') " Ray         Z exit         Z exit analytic     Z error    Z error bar"
      write (fileUnit,'(a)') " ---------------------------------------------------------------------"

      if (sim_geometry == GRID_3DCARTESIAN) then
          do ray = 1,4
             rayZ         = sim_rayZexit (ray)
             rayZanalytic = sim_zw
             rayZerror    = abs (rayZ - sim_zw) * 100.0 / sim_zw
             rayZerrorBar = sim_rayFexitPercentError1
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayZ, rayZanalytic, rayZerror, rayZerrorBar
          end do
          do ray = 5,8
             rayZ         = sim_rayZexit (ray)
             rayZanalytic = sim_zw
             rayZerror    = abs (rayZ - sim_zw) * 100.0 / sim_zw
             rayZerrorBar = sim_rayFexitPercentError12
             write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayZ, rayZanalytic, rayZerror, rayZerrorBar
          end do
      end if

  end if
!
!
!     ...Print out the rays exit power.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         RAY P EXIT (errors in %)"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray         P exit         P exit analytic     P error    P error bar"
  write (fileUnit,'(a)') " ---------------------------------------------------------------------"

  if (sim_geometry == GRID_2DCARTESIAN .or. sim_geometry == GRID_2DCYLINDRICAL) then
      do ray = 1,2
         rayP         = sim_rayPexit (ray)
         rayPanalytic = sim_powerDecayFactor
         rayPerror    = abs (rayP - sim_powerDecayFactor) * 100.0 / sim_powerDecayFactor
         rayPerrorBar = sim_rayPexitPercentError1
         write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayP, rayPanalytic, rayPerror, rayPerrorBar
      end do
  else if (sim_geometry == GRID_3DCARTESIAN) then
      do ray = 1,4
         rayP         = sim_rayPexit (ray)
         rayPanalytic = sim_powerDecayFactor
         rayPerror    = abs (rayP - sim_powerDecayFactor) * 100.0 / sim_powerDecayFactor
         rayPerrorBar = sim_rayPexitPercentError1
         write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayP, rayPanalytic, rayPerror, rayPerrorBar
      end do
      do ray = 5,8
         rayP         = sim_rayPexit (ray)
         rayPanalytic = sim_powerDecayFactor
         rayPerror    = abs (rayP - sim_powerDecayFactor) * 100.0 / sim_powerDecayFactor
         rayPerrorBar = sim_rayPexitPercentError12
         write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayP, rayPanalytic, rayPerror, rayPerrorBar
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
