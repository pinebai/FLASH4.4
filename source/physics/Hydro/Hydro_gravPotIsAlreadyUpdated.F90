!!****f* source/physics/Hydro/Hydro_gravPotIsAlreadyUpdated
!!
!! NAME
!!
!!  Hydro_gravPotIsAlreadyUpdated
!!
!! SYNOPSIS
!!
!!  logical result = Hydro_gravPotIsAlreadyUpdated()
!!
!! DESCRIPTION
!!
!!   This function returns an indication whether the most
!!   recent call to Hydro has, as a side effect, completely
!!   updated the gravitational potential in the UNK variable
!!   GPOT (and other state related to gravitational forces
!!   if applicable - see implementation of sink paticles).
!!
!!   Code that is responsible for updating GPOT to the most
!!   recent solution of the Poisson equation for gravity
!!   (typically by calling Gravity_potentialListOfBlocks
!!   at the end of each time advance iteration in
!!   Driver_evolveFlash) MAY use this function to avoid
!!   an unnecessary additional call to the Poisson solver.
!!
!!   An implementation of this interface should only return
!!   TRUE when it is certain that an additional call to
!!   Gravity_potentialListOfBlocks would not change
!!   results.
!!
!!   For any constant-in-time Gravity implementations,
!!   the return value does not matter.
!!   
!! ARGUMENTS
!!
!!   No arguments
!!
!! SEE ALSO
!!
!!  Driver_evolveFlash
!!  Hydro
!!  Gravity_potentialListOfBlocks
!!  Particles_advance
!!
!!***

logical function Hydro_gravPotIsAlreadyUpdated()
  implicit none

  Hydro_gravPotIsAlreadyUpdated = .FALSE.

end function Hydro_gravPotIsAlreadyUpdated
