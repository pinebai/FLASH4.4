!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/charlette/gammaInt/fl_fsTFIFunc
!!
!! NAME
!!
!!  fl_fsTFIFunc
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIFunc(real(in) :: x)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   x : 
!!
!!
!!
!!***

#include "constants.h"

real function fl_fsTFIFunc(x)
   use fl_fsTFIData, ONLY : fl_fsTFIus, fl_fsTFIdedl, fl_fsTFIPrandtl
   implicit none
   real, intent(in) :: x

   real :: Re, Ck, a, h, arg1, arg2, Ceff, Ccdvp

   Ck = 1.5
   a = 18.0/55.0

   Re = fl_fsTFIus * fl_fsTFIdedl / fl_fsTFIPrandtl

   ! h function
   h = exp(-1.5*Ck*(x*PI)**(4.0/3.0)/Re)

   ! Ccdvp function
   arg1 = fl_fsTFIdedl / x
   arg2 = fl_fsTFIus*(x*PI)**(-1.0/3.0)*sqrt(a*Ck*h)
   Ccdvp = 0.5*(1.0 + erf(0.6*(log(arg1) - 1.0/sqrt(arg2))))

   ! C function
   Ceff = Ccdvp * 0.5 * (1.0 + erf(3.0*log10(2.0*arg2)))

   ! function f
   fl_fsTFIFunc = a * Ck * PI**(4.0/3.0) * Ceff**2 * x**(1.0/3.0) * h

end function fl_fsTFIFunc
