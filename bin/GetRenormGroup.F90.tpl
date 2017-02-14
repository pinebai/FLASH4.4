##
## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for Generating an F90 subroutine
## For syntax of this file see "Readme.template"
##
## Map MassScalars to Group Name
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## mscalars -> list of integers
## groups -> list of integers
## group_names -> piece of text mapping group numbers to names 
##
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_getRenormGroup(mscalar,group)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(out) ::group 
   integer, intent(in) :: mscalar
   integer :: ctr

   !! values is an array of possible input strings (with possibly extra spaces)
   !! groups is an array of corresponding answers

   integer, dimension(1:%(COUNT_mscalars|1)s), parameter :: values = (/&
   & %(5mscalars !, !,&\n     |0)s &
   &/)

   integer, dimension(1:%(COUNT_mscalars|1)s), parameter :: groups = (/&
   & %(5groups !, !,&\n     |0)s &
   &/)

   !! The group numbers correspond to the following names used in the Config files
   !! %(group_names)s
   !! This is for informational purposes only

   group = 0  ! By default group 0 doesn't exist
   do ctr = LBOUND(values,1), UBOUND(values,1)
      if ( values(ctr) .eq. mscalar ) group = groups(ctr)
   end do

end subroutine Simulation_getRenormGroup

