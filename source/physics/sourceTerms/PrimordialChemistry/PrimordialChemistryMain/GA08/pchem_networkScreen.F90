!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_networkScreen
!!
!! NAME
!!
!! pchem_networkScreen
!!
!! SYNOPSIS
!!
!! call pchem_networkScreen(real intent(IN) :: y)
!!
!! DESCRIPTION
!!
!! routine pchem_networkScreen applies screening corrections to the raw reaction
!! rates
!!
!! this routine computers screening factors
!! and applies them to the raw reaction rates,
!! producing the final raction rates used by the 
!! right hand sides and jacobian matrix elements
!! input is y is ymass. Does not actually do anything besides rename ratdum-->ratraw
!!
!! ARGUMENTS
!!
!! y : ymass
!!
!!****

subroutine pchem_networkScreen(y)

   use PrimordialChemistry_data
   use pchem_data
   use pchem_interface !!I guess I need to write this too

   implicit none

#include "constants.h"
#include "Flash.h"

   !! declare
   real, intent(IN) :: y(NSPECIES)

   !! local declarations
   integer              :: i, jscr, screen_init
   integer, parameter    :: screen_on = 1
   real    sc1

   !! Since I am not sure if there is any screening involved. I will make
   !! this real simple and have it all be 1.0

   do i=1,nrat
        ratdum(i) = ratraw(i)
        scfac(i) = 1.0e0
   enddo
   return

end subroutine pchem_networkScreen