!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_setupRays
!!
!! NAME
!!
!!  ed_setupRays
!!
!! SYNOPSIS
!!
!!  call ed_setupRays ()
!!
!! DESCRIPTION
!!
!!  Sets up the rays array. This simply means to allocate the needed rays array.
!!  Since the rays will be created afresh at each time step, there is no point in
!!  initializing the rays array at this stage. Initialization of the rays is done
!!  right before they are going to be created.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_setupRays ()

  implicit none

  return
end subroutine ed_setupRays
