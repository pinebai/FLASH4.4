!!****f* source/Simulation/Simulation_mapIntToStr
!!
!! NAME
!!  Simulation_mapIntToStr
!!
!!
!! SYNOPSIS
!!  Simulation_mapIntToStr(integer, intent(IN)             :: key, 
!!                         character, len=*, intent(INOUT) :: str,
!!                         integer, intent(IN)             :: map)
!!
!!
!! DESCRIPTION
!!
!!  This routine is created by the setup script, and should never be edited.
!!  This routine maps an index described in the Flash.h file to a string
!!  described in the Config file.  The integer can represent a variable,
!!  a species, a flux variable, or a particle property.  
!!
!!  For example, the Config file might say:
!!  VARIABLE velx
!!  At setup time, the Flash.h file would be created to look like
!!  #define VELX_VAR 6
!!  The result of "call Simulation_mapIntToStr(VELX_VAR,result,MAPBLOCK_UNK)" would
!!  be result="velx"    
!!
!! ARGUMENTS
!! 
!!  key   --  integer index
!!  str   --  returned string
!!  map   --  variable indicating the type of data structure within Flash.h.  Valid values are
!!            MAPBLOCK_UNK   for variables and species (NAME_VAR or NAME_SPEC in Flash.h)
!!            MAPBLOCK_FLUX  for flux variables (NAME_FLUX in Flash.h)          
!!            MAPBLOCK_PART  for particle properties (NAME_PART_PROP in Flash.h)
!!
!!***

subroutine Simulation_mapIntToStr(key, str, block)
implicit none 

#include "constants.h"

   integer, intent(in) :: key, block
   character(len=*), intent(inout) :: str
  
   str = "ERROR"

end subroutine Simulation_mapIntToStr


