##
## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for Generating an F90 subroutine
## For syntax of this file see "Readme.template"
##
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Particles_specifyMethods()

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  use Particles_data, ONLY : pt_typeInfo
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer :: n


  !! Maintain 3 vectors that hold information about the particles.
  !! The information is held in the same order in all vectors, e.g. 
  !! element 1 in all vectors refers to the same particle type.

  integer, dimension(1:NPART_TYPES) :: &
     particleTypeVect, particleInitVect, particleMapVect, particleAdvVect

  particleTypeVect(:) = NONEXISTENT
  particleInitVect(:) = NONEXISTENT
  particleMapVect(:)  = NONEXISTENT
  particleAdvVect(:)  = NONEXISTENT


  !! Include a failsafe mechanism in case the user includes 
  !! particle units, but does not specify a PARTICLETYPE line.
  if (NPART_TYPES == 0) then
    call Driver_abortFlash("[Particles_specifyMethods]: No PARTICLETYPE Config line")
  end if


  do n = 1, NPART_TYPES
    particleTypeVect(n) = PART_TYPES_BEGIN + (n - 1)
  end do 

  n = 1
  particleInitVect(n) = %(particleInitVect!\n  n = n + 1\n  particleInitVect(n) = |NONEXISTENT)s

  n = 1
  particleMapVect(n) = %(particleMapVect!\n  n = n + 1\n  particleMapVect(n) = |NONEXISTENT)s

  n = 1
  particleAdvVect(n) = %(particleAdvVect!\n  n = n + 1\n  particleAdvVect(n) = |NONEXISTENT)s


  !! Collate the information from our vectors into the main
  !! particle data structure.

  pt_typeInfo(PART_TYPE,1:NPART_TYPES) = particleTypeVect(1:NPART_TYPES)
  pt_typeInfo(PART_INITMETHOD,1:NPART_TYPES) = particleInitVect(1:NPART_TYPES)
  pt_typeInfo(PART_MAPMETHOD,1:NPART_TYPES) = particleMapVect(1:NPART_TYPES)
  pt_typeInfo(PART_ADVMETHOD,1:NPART_TYPES) = particleAdvVect(1:NPART_TYPES)

end subroutine Particles_specifyMethods
