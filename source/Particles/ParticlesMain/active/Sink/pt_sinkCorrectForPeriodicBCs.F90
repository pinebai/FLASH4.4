!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkCorrectForPeriodicBCs
!!
!! NAME
!!
!!  pt_sinkCorrectForPeriodicBCs
!!
!! SYNOPSIS
!!
!!  call pt_sinkCorrectForPeriodicBCs(real(inout)  :: distx,
!!                                    real(inout)  :: disty,
!!                                    real(inout)  :: distz)
!!
!! DESCRIPTION
!!
!!  For a given distance (x,y,z) within the computational box,
!!  this function corrects the distances for periodic boundary conditions.
!!
!! ARGUMENTS
!!
!!   distx - x distance between two points in the computational domain
!!
!!   disty - y distance between two points in the computational domain
!!
!!   distz - z distance between two points in the computational domain
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2014
!!
!!***

subroutine pt_sinkCorrectForPeriodicBCs(distx, disty, distz)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  real, intent(INOUT) :: distx, disty, distz

  character(len=80), save :: grav_boundary_type
  logical, save           :: first_call = .true.
  real, save              :: xmin, xmax, ymin, ymax, zmin, zmax
  real, save              :: Lx, Ly, Lz, LxHalf, LyHalf, LzHalf

  if (first_call) then
     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
     if ((grav_boundary_type .ne. "isolated") .and. (grav_boundary_type .ne. "periodic")) then
        call Driver_abortFlash("Sink particles can only be used with periodic or isolated gravity type!")
     end if
     call RuntimeParameters_get('xmin', xmin)
     call RuntimeParameters_get('xmax', xmax)
     call RuntimeParameters_get('ymin', ymin)
     call RuntimeParameters_get('ymax', ymax)
     call RuntimeParameters_get('zmin', zmin)
     call RuntimeParameters_get('zmax', zmax)
     Lx = xmax-xmin
     Ly = ymax-ymin
     Lz = zmax-zmin
     LxHalf = 0.5*Lx
     LyHalf = 0.5*Ly
     LzHalf = 0.5*Lz
     first_call = .false.
  endif

  if (distx .lt. -LxHalf) distx = distx+Lx
  if (distx .gt. +LxHalf) distx = distx-Lx
  if (disty .lt. -LyHalf) disty = disty+Ly
  if (disty .gt. +LyHalf) disty = disty-Ly
  if (distz .lt. -LzHalf) distz = distz+Lz
  if (distz .gt. +LzHalf) distz = distz-Lz

  return

end subroutine pt_sinkCorrectForPeriodicBCs
