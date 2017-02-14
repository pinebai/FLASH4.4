!!****if* source/Particles/localAPI/pt_findTagOffset
!!
!! NAME
!!    pt_findTagOffset
!!
!! SYNOPSIS
!!
!!    pt_findTagOffset()
!!
!! DESCRIPTION
!!    Add a unique tag each to particles generated during evolution.
!!    The algorithm first finds out the sum of number of particles 
!!    being added to all the processors to the left of MyPE. It then
!!    adds pt_startTagNumber (the largest tag number in use) to it.
!!    This number acts as the offset for the tag numbers being assigned
!!    to the newly added particles.
!!
!! NOTES
!!    This method of tag generation will work for up to 10^14 
!!    particles in a simulation
!!
!!
!!
!!
!!***

!!#define DEBUG_PARTICLES

subroutine pt_findTagOffset(newCount,tagOffset)

  implicit none

  integer, intent(IN) :: newCount
  integer, intent(OUT) :: tagOffset
  
  tagOffset=0

  return
end subroutine pt_findTagOffset
