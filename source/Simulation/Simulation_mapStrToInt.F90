!!****f* source/Simulation/Simulation_mapStrToInt
!!
!! NAME
!!  Simulation_mapStrToInt
!!
!!
!! SYNOPSIS
!!  Simulation_mapStrToInt(character, len=*, intent(IN) :: str,
!!                         integer, intent(OUT)         :: key, 
!!                         integer, intent(IN)          :: map)
!!
!!
!! DESCRIPTION
!!
!!  This routine is created by the setup script, and should never be edited.
!!  This routine maps a string described in the Config file to an integer
!!  index described in the Flash.h file.  The integer can represent a variable,
!!  a species, a flux variable, or a particle property.  
!!
!!  For example, the Config file might say:
!!  VARIABLE velx
!!  At setup time, the Flash.h file would be created to look like
!!  #define VELX_VAR 6
!!  The result of "call Simulation_mapStrToInt(keyout,'velx',MAPBLOCK_UNK)" would
!!  be keyout=6.    
!!
!! ARGUMENTS
!! 
!!  str   --  string input
!!  key   --  returned integer index
!!  map   --  variable indicating the type of data structure within Flash.h.  Valid values are
!!            MAPBLOCK_UNK   for variables and species (NAME_VAR or NAME_SPEC in Flash.h)
!!            MAPBLOCK_FLUX  for flux variables (NAME_FLUX in Flash.h)          
!!            MAPBLOCK_PART  for particle properties (NAME_PART_PROP in Flash.h)
!!
!!***

subroutine Simulation_mapStrToInt(str,key,map)
implicit none 

#include "constants.h"

   character(len=*), intent(in) :: str
   integer, intent(out) :: key 
   integer, intent(IN) :: map

   key = NONEXISTENT

end subroutine Simulation_mapStrToInt

