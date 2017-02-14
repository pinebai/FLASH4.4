!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_intSimpson
!!
!! NAME
!!
!!  grv_intSimpson
!!
!!
!! SYNOPSIS
!!
!!   real = grv_intSimpson(
!!                           real, external(in)    :: fce,
!!                           integer(in)           :: l,
!!                           real(in)              :: eta,
!!                           real(in)              :: dzeta,
!!                           real(in)              :: xmin,
!!                           real(in)              :: xmax,
!!                           integer(in)           :: n
!!        )
!!
!!
!! DESCRIPTION
!!
!!   grv_intSimpson integrates given function by Simpson method.
!!   The function has three parameters: two are of type real and
!!   the third is of type integer.
!!   
!!
!! ARGUMENTS
!!
!!  fce         - name of the integrated function
!!  l           - parameter of the function which is related to wavenumber
!!  eta, dzeta  - parameters of the function which are related to the Ewald method
!!                (see paragraph 2.2.3 in Wunsch et al. 2015)
!!  xmin, xmax  - endpoints of the interval where the function is integrated
!!  n           - number of bins used for numerical integration
!!
!!
!! RESULT
!!
!!  Value of the integral
!!
!! NOTES
!!
!!***

#include "Flash.h"

real function grv_intSimpson(fce,l,eta,dzeta,xmin,xmax,n)

#if defined(__INTEL_COMPILER)
  use IFPORT
#endif

  implicit none

  real,external :: fce
  real,intent(in) :: eta,dzeta,xmin,xmax
  integer,intent(in) :: l,n
  real x,dx
  real s,f0,f1,f2
  integer k

  s = 0.0
  x=xmin
  dx=0.5*(xmax-xmin)/real(n)

! integration by Simpson method
  do k=0,n
    f0=fce(x,l,eta,dzeta)
    f1=fce(x+dx,l,eta,dzeta)
    f2=fce(x+2.0*dx,l,eta,dzeta)
    s=s+(f0+f2+4.0*f1)
    x=x+2.0*dx
  enddo

  grv_intSimpson=s*dx/3.0

end

