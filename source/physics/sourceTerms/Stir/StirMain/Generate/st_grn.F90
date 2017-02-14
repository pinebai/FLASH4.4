!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/st_grn
!!
!! NAME
!!
!!  st_grn
!!
!! SYNOPSIS
!!
!!  st_grn(real,intent (OUT)  :: grnval)
!!
!! DESCRIPTION
!!
!!  Subroutine draws a number randomly from a Gaussian distribution
!!    with the standard uniform distribution function "random_number"
!!    using the Box-Muller transformation in polar coordinates. The 
!!    resulting Gaussian has unit variance.
!!
!!  Reference : Numerical Recipes, section 7.2.
!!
!!
!!
!! ARGUMENTS
!!
!!   grnval : the number drawn from the distribution
!!
!!
!!
!!***


subroutine st_grn (grnval)
 
  use Stir_data, ONLY : st_reproducible, st_saveReproducible, st_randomSaveUnit
  use ut_randomInterface, ONLY : ut_randomNumber
  implicit none
  
#include "constants.h"

  real,intent (OUT) :: grnval 
  real              :: r1, r2, g1


  r1 = 0.; r2 = 0;
  if(.not.st_reproducible) then
     call ut_randomNumber (r1)
     call ut_randomNumber (r2)
  else
     read(st_randomSaveUnit,*)r1
     read(st_randomSaveUnit,*)r2
  end if
  if(st_saveReproducible) then
     write(st_randomSaveUnit,*)r1
     write(st_randomSaveUnit,*)r2
  end if

  g1 = sqrt (2. * log (1. / r1) ) * cos (2. * PI * r2)
  
  grnval = g1
  
  return 
  
end subroutine st_grn

