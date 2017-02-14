##
## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating an F90 subroutine
## For syntax of this file see "Readme.template"
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## partkeys -> list of particle properties
## varkeys -> list of variable names
## vartypes -> list of variable types
##
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!

subroutine Simulation_mapParticlesVar(part_key, var_key, var_type)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(in)  :: part_key
   integer, intent(out) :: var_key, var_type


   integer :: ctr, tmp_var_key, tmp_var_type

   integer, dimension(1:%(COUNT_partkeys|1)s), parameter :: part_keys = (/&
   & %(partkeys!,&\n     |NONEXISTENT)s&
   &/)

   integer, dimension(1:%(COUNT_varkeys|1)s), parameter :: var_keys = (/&
   & %(varkeys!,&\n      |NONEXISTENT)s&
   &/)

   integer, dimension(1:%(COUNT_vartypes|1)s), parameter :: var_types = (/&
   & %(vartypes!,&\n      |NONEXISTENT)s&
   &/)

   tmp_var_key = NONEXISTENT
   tmp_var_type = NONEXISTENT

   do ctr = LBOUND(part_keys,1), UBOUND(part_keys,1)
      if (part_keys(ctr) .eq. part_key) then
         tmp_var_key = var_keys(ctr)
         tmp_var_type = var_types(ctr)
      end if
   end do
   var_key = tmp_var_key
   var_type = tmp_var_type

end subroutine Simulation_mapParticlesVar


