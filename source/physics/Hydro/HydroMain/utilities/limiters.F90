!!****ih* source/physics/Hydro/HydroMain/utilities/limiters
!!
!! NAME
!!   limiters - modules for providing limiters
!!
!! SYNOPSIS
!!   use limiters
!!
!! DESCRIPTION
!!   This module provides limiters for use in shock-capturing schemes.
!!   Currently only a generalized minmod limiter is provided, but others
!!   could be added with little effort.
!!
!!***

module limiters

   implicit none

   interface minmod

!     Separate versions with and without theta as input:
      module procedure scalar_minmod2, scalar_minmod2t,           &
                       vector_minmod2, vector_minmod2t

!!!   Optional arguments versions:
!!    module procedure scalar_minmod2t, vector_minmod2t

   end interface

   contains

!!****if* /source/hydro/explicit/minmod
!!
!! NAME
!!   minmod - limiter for shock capturing schemes
!!
!! SYNOPSIS
!!   c   = minmod( a,    b    )
!!         minmod( real, real )
!!
!!   c() = minmod( a(),    b()    )
!!         minmod( real(), real() )
!!
!!   c   = minmod( a,    b,    theta )
!!         minmod( real, real, real  )
!!
!!   c() = minmod( a(),    b(),    theta )
!!         minmod( real(), real(), real  )
!!
!! DESCRIPTION
!!   The minmod function in its most basic form takes at least two numbers and
!!   returns the smallest if all are positive,
!!           the largest  if all are negative,
!!           zero         if they are not all positive or all negative
!!
!!   The first two forms, above, apply the minmod function to two numbers only.
!!   In the first form, two reals are input, and the result is a real.
!!   In the second form, two vectors of reals are input, and the minmod
!!     function is applied to each pair of elements. The result is a vector.
!!
!!   The generalized minmod function accepts two numbers, a and b, and
!!   a parameter, theta, and returns
!!      minmod(theta a, (a+b)/2, theta b)
!!
!!   theta is restricted to be between (inclusive) 1.0 and 2.0.
!!   theta = 1.0 for minmod limiter
!!   theta = 2.0 for monotonized centered limiter
!!   When used for shock-capturing schemes, larger implies less dissipation.
!!
!!   The third and fourth forms are the scalar and vector versions of the
!!     generalized minmod function. (Note the same theta is used for all
!!     pairs in the vector version.)
!! 
!! ARGUMENTS
!!   All arguments are intent(IN).
!!   a 
!!   b
!!   theta  The blending parameter, which adjusts the dissipation of the
!!          (most dissipation)  1.0 <= theta <= 2.0  (least)
!!          theta = 1.0 for minmod limiter
!!          theta = 2.0 for monotonized centered limiter
!!
!!***
!23456789*123456789*123456789*123456789*123456789*123456789*123456789*12
   function scalar_minmod2(a,b)
     implicit none
     real, intent(in)  :: a, b
!!     real, intent(out) :: scalar_minmod2
     real              :: scalar_minmod2
 
     scalar_minmod2 =                               &
       0.5*(  sign(1.0e0,a)                         &
            + sign(1.0e0,b))*min(abs(a),abs(b))
 
   end function scalar_minmod2
!_______________________________________________________________________

!23456789*123456789*123456789*123456789*123456789*123456789*123456789*12
   function scalar_minmod2t(a,b,theta)
     implicit none
     real, intent(in)   :: a, b, theta
!!     real, intent(out)  :: scalar_minmod2t
     real               :: scalar_minmod2t
     real               :: t, left, middle, right
     real               :: smmono, smmag
 
     t = min( 2.0e0, max(1.0e0,theta) )  ! restrict 1.0 < t < 2.0
     left = t*a
     right = t*b
     middle = 0.5e0*(a+b)
     smmono = 0.5e0*(  sign( 1.0e0, max(left,middle,right) )     &
                     + sign( 1.0e0, min(left,middle,right) ) )
     smmag = min( abs(left), abs(middle), abs(right) )
     scalar_minmod2t = smmag*smmono
 
   end function scalar_minmod2t
!_______________________________________________________________________

!23456789*123456789*123456789*123456789*123456789*123456789*123456789*12
   function vector_minmod2(a,b)
     implicit none
     real, intent(in)   :: a(:), b(:)
!!     real, intent(out)  :: vector_minmod2( max(size(a),size(b)) )
     real               :: vector_minmod2( max(size(a),size(b)) )
 
     vector_minmod2 =                               &
       0.5*(  sign(1.0e0,a)                         &
            + sign(1.0e0,b))*min(abs(a),abs(b))
 
   end function vector_minmod2
!_______________________________________________________________________

!23456789*123456789*123456789*123456789*123456789*123456789*123456789*12
   function vector_minmod2t(a,b,theta)
     implicit none
     real, intent(in)  :: a(:), b(:)
     real, intent(in)  :: theta
!!     real, intent(out) :: vector_minmod2t( max(size(a),size(b)) )
     real              :: vector_minmod2t( max(size(a),size(b)) )
     real              :: t
     real, dimension( max(size(a),size(b)) ) :: left, middle, right
     real, dimension( max(size(a),size(b)) ) :: smmono, smmag
 
     t = min( 2.0e0, max(1.0e0,theta) )  ! restrict 1.0 < t < 2.0
     left   = t*a
     right  = t*b
     middle = 0.5e0*(a+b)
     smmono = 0.5e0*(  sign( 1.0e0, max(left,middle,right) )     &
                     + sign( 1.0e0, min(left,middle,right) ) )
     smmag = min( abs(left), abs(middle), abs(right) )
     vector_minmod2t = smmag*smmono
 
   end function vector_minmod2t
!_______________________________________________________________________

!
!  Below are versions that use optional arguments.
!

!!$!23456789*123456789*123456789*123456789*123456789*123456789*123456789*12
!!$   function scalar_minmod2t(a,b,theta)
!!$     implicit none
!!$     real, intent(in)               :: a, b
!!$     real, intent(in), optional     :: theta
!!$!!     real, intent(out)              :: scalar_minmod2t
!!$     real                           :: scalar_minmod2t
!!$     real                           :: t, left, middle, right
!!$     real                           :: smsign, smmono, smmag
!!$ 
!!$     if(.not.present(theta)) then
!!$
!!$       scalar_minmod2t =                                 &
!!$         0.5*(  sign(1.0e0,a)                            &
!!$              + sign(1.0e0,b))*min(abs(a),abs(b))
!!$ 
!!$     else
!!$ 
!!$       t = min( 2.0e0, max(1.0e0,theta) )  ! restrict 1.0 < t < 2.0
!!$       left = t*a
!!$       right = t*b
!!$       middle = 0.5e0*(a+b)
!!$       smmono = 0.5e0*(  sign( 1.0e0, max(left,middle,right) )     &
!!$                       + sign( 1.0e0, min(left,middle,right) ) )
!!$       smmag = min( abs(left), abs(middle), abs(right) )
!!$       scalar_minmod2t = smmag*smmono
!!$ 
!!$     endif
!!$ 
!!$   end function scalar_minmod2t
!!$!_______________________________________________________________________

!!$!23456789*123456789*123456789*123456789*123456789*123456789*123456789*12
!!$   function vector_minmod2t(a,b,theta)
!!$     implicit none
!!$     real, intent(in)            :: a(:), b(:)
!!$     real, intent(in), optional  :: theta
!!$!!     real, intent(out) :: vector_minmod2t( max(size(a),size(b)) )
!!$     real              :: vector_minmod2t( max(size(a),size(b)) )
!!$     real              :: t
!!$     real, dimension( max(size(a),size(b)) ) :: left, middle, right
!!$     real, dimension( max(size(a),size(b)) ) :: smmono, smmag
!!$ 
!!$     if(.not.present(theta)) then
!!$ 
!!$       vector_minmod2t =                                 &
!!$         0.5*(  sign(1.0e0,a)                            &
!!$              + sign(1.0e0,b))*min(abs(a),abs(b))
!!$ 
!!$     else
!!$ 
!!$       t = min( 2.0e0, max(1.0e0,theta) )  ! restrict 1.0 < t < 2.0
!!$       left   = t*a
!!$       right  = t*b
!!$       middle = 0.5e0*(a+b)
!!$       smmono = 0.5e0*(  sign( 1.0e0, max(left,middle,right) )     &
!!$                       + sign( 1.0e0, min(left,middle,right) ) )
!!$       smmag = min( abs(left), abs(middle), abs(right) )
!!$       vector_minmod2t = smmag*smmono
!!$ 
!!$     endif
!!$ 
!!$   end function vector_minmod2t
!!$!_______________________________________________________________________


end module limiters
!_______________________________________________________________________

