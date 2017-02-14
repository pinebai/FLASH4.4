!!****if* source/Particles/localAPI/pt_mapStringParamToInt
!!
!! NAME
!!    pt_mapStringParamToInt
!!
!! SYNOPSIS
!!
!!    pt_mapStringParamToInt(integer(OUT):: mappedInt,
!!                           character(MAX_STRING_LENGTH)(IN)::paramString,
!!                           integer(IN) :: mapblock,
!!                           integer(IN) :: ind)
!!
!! DESCRIPTION
!!
!!  This routine is used for finding the integer value associated 
!!  with some string based runtime particles parameters that occur in groups
!!  and are therefore enumerated. 
!!  For example for runtimte parameters such as "particles_attributes_1", 
!!  "particles_attributes_2" etc, this routine takes paramString as 
!!   "particles_attributes_", ind provides the enumeration value, and the
!!  associated constant is returned in mappedInt
!!
!! ARGUMENTS
!! 
!!    mappedInt - the integer constant associated with the parameter
!!    paramString - the basename for the string based parameter
!!    mapblock - a reference to the distinguishing marker used by
!!               the setup script to know where to look for the 
!!               associated constant
!!    ind     - eunmeration index
!!
!!
!!***

subroutine pt_mapStringParamToInt(mappedInt,paramString, mapblock,ind  )

#include "constants.h"
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none
  integer,intent(OUT) :: mappedInt
  character(len=MAX_STRING_LENGTH),intent(IN) :: paramString
  integer,intent(IN) :: mapblock
  integer,intent(IN) :: ind

end subroutine pt_mapStringParamToInt
