!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_slopeLimiters
!!
!! NAME
!!
!!  ed_slopeLimiters
!!
!! SYNOPSIS
!!
!!  use ed_slopeLimiters
!!
!! ARGUMENTS
!!
!! DESCRIPTION
!!
!!  This module stores limiter functions that are specific to 
!!  the EnergyDeposition routines.
!!
!!***


Module ed_slopeLimiters

implicit none

contains

  function ed_checkMedian (a,b,c)
    implicit none
    real :: a,b,c,ed_checkMedian
    ed_checkMedian = a + ed_minmod (b-a,c-a)
  end function ed_checkMedian


  function ed_vanLeer (a,b)
    implicit none
    real :: a,b,ed_vanLeer
    if (a*b <=0.) then
        ed_vanLeer = 0.
    else
        ed_vanLeer = 2.*a*b/(a+b)
    endif
  end function ed_vanLeer


  function ed_mc (a,b)
    implicit none
    real :: a,b,ed_mc
    ed_mc = (sign(1.,a)+sign(1.,b))*min(.5*abs(a),.25*abs(a+b),.5*abs(b))
  end function ed_mc


  function ed_minmod (a,b)
    implicit none
    real :: a,b,ed_minmod
    ed_minmod = .5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
  end function ed_minmod


  function ed_signum (x)
    implicit none
    real :: x,ed_signum
    ed_signum = sign(.5,x)-sign(.5,-x)
  end function ed_signum


  function ed_firstDeriv (Umm,Um,U0,Up,Upp,order,upwindDir)
    implicit none
    real :: Umm,Um,U0,Up,Upp
    integer :: order,upwindDir
    real :: ed_firstDeriv

    if ( upwindDir > 0) then
       select case(order)
       case(1)
          ed_firstDeriv = U0-Um
       case(2)
          ed_firstDeriv = 0.5*(3.*U0-4.*Um+Umm)
       case(3)
          ed_firstDeriv = (2.*Up+3.*U0-6.*Um+Umm)/6.
       end select
    elseif (upwindDir < 0) then
       select case(order)
       case(1)
          ed_firstDeriv = Up-U0
       case(2)
          ed_firstDeriv = 0.5*(-Upp+4.*Up-3.*U0)
       case(3)
          ed_firstDeriv = (-Upp+6.*Up-3.*U0-2.*Um)/6.
       end select

    endif
  end function ed_firstDeriv

End Module ed_slopeLimiters
