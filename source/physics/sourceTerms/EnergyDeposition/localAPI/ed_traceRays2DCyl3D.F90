!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_traceRays2DCyl3D
!!
!! NAME
!!
!!  ed_traceRays2DCyl3D
!!
!! SYNOPSIS
!!
!!  call ed_traceRays2DCyl3D (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Processes all 3D rays on the collection of 2D cylindrical blocks on the current processor.
!!  On exit, each ray in each block has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!! ARGUMENTS
!!
!!  timeStep : Current timestep value
!!
!! NOTES
!!
!!  The use of threading is possible for tracing the rays through each block. 
!!
!!***

subroutine ed_traceRays2DCyl3D (timeStep)

  implicit none

  real, intent (in) :: timeStep

  return
end subroutine ed_traceRays2DCyl3D
