##
## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating an F90 subroutine
## For syntax of this file see "Readme.template"
##
## Map Variable Name to Type
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## varnames -> list of strings (in Uppercase)
## vartypes -> list of strings (in Uppercase)
##
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_getVarnameType(varname,vartype)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(out) :: vartype
   integer, intent(in) :: varname
   integer :: ctr

   !! varnames is an array of possible input strings (with possibly extra spaces)
   !! vartypes is an array of corresponding answers

   integer, dimension(1:%(COUNT_varnames|1)s), parameter :: varnames = (/&
   & %(5varnames!, !, &\n     |VARTYPE_ERROR)s &
   &/)

   integer, dimension(1:%(COUNT_vartypes|1)s), parameter :: vartypes = (/&
   & %(5vartypes!, !,&\n     |VARTYPE_ERROR)s &
   &/)

   ! Species and mass scalars are kept in primitive form. - KW
   if (varname .ge. SPECIES_BEGIN .and. varname .le. MASS_SCALARS_END) then
      vartype = VARTYPE_PER_MASS
      return
   end if

   vartype = VARTYPE_ERROR  ! By default vartype is an ERROR
   do ctr = LBOUND(varnames,1), UBOUND(varnames,1)
      if ( varnames(ctr) .eq. varname ) vartype = vartypes(ctr)
   end do

end subroutine Simulation_getVarnameType

