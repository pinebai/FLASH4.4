!!****if* source/flashUtilities/interpolation/oneDim/ut_quadraticInterpol
!!
!!  NAME
!!    ut_quadraticInterpol
!!  SYNOPSIS 
!!    call ut_quadraticInterpol(x, y, xint, yint)
!!  FUNCTION 
!!    Given 3-element arrays (x,y) of places and values, and an interpolation point xint,
!!    returns a quadratic interpolant (or extrapolant -- no range checking done!)
!!    of the function to yint.
!!  INPUTS
!!      x  - real, sorted, array of cordinates
!!      y  - function values
!!      xint  - value to interpolate at
!!      yint  - function value at xint
!!  RESULT
!!      returns yint(xint; x,y)
!!
!!  HISTORY
!!
!!    This was in FLASH2 as setups/nova/non-numrecipies/polint.F90.
!!***
      subroutine ut_quadraticInterpol(x, y, xint, yint)

      implicit none
      real, intent(in)    :: xint           ! where function is to be eval'd
      real, intent(in)    :: y(*), x(*)     ! function values and locations
      real, intent(out)   :: yint           ! function evaluated at xint


      real :: dx12, dx23, dx13
      real :: a, b, c


      dx12 = x(1)-x(2)
      dx23 = x(2)-x(3)
      dx13 = x(1)-x(3)

      a = (xint-x(2))*(xint-x(3))/(dx12*dx13)
      b = (xint-x(1))*(xint-x(3))/(-dx12*dx23)
      c = (xint-x(1))*(xint-x(2))/(dx13*dx23)

      yint = y(1)*a + y(2)*b + y(3)*c
      return
      end subroutine ut_quadraticInterpol



!!****if* setups/nova/ut_quadraticCellAverageInterpol
!!
!!  NAME
!!    ut_quadraticCellAverageInterpol
!!  SYNOPSIS 
!!    call ut_quadraticInterpol(x, xl, y, xint, yint)
!!  FUNCTION 
!!    Given 3-element arrays (x,y) of locations and values, and an interpolation
!!    point xint, and the cell boundaries xl (the left-side cell wall),
!!    returns a quadratic interpolant (or extrapolant -- no range checking done!)
!!    of the function to yint such that the function doesn't interpolate y(), but
!!    the interpolant would integrate to cell-averaged quantities y.
!!  INPUTS
!!      x  - real, sorted, array of cordinates
!!      y  - function values
!!      xint  - value to interpolate at
!!      yint  - function value at xint
!!  RESULT
!!      returns yint(xint; x,y)
!!
!!  HISTORY
!!
!!    This was named polint_cav in FLASH2.
!!***
      subroutine ut_quadraticCellAverageInterpol(x, xl, y, xint, yint)

      implicit none
      real, intent(in)    :: xint           ! where function is to be eval'd
      real, intent(in)    :: y(*), x(*), xl(*) ! function values and locations
      real, intent(out)   :: yint           ! function evaluated at xint


      real :: a, b, c, Ai, Bi, Ci


      Ai = 1./3. * (xl(2)*xl(2)*xl(2)-xl(1)*xl(1)*xl(1))/(xl(2)-xl(1))
      Bi = 1./3. * (xl(3)*xl(3)*xl(3)-xl(2)*xl(2)*xl(2))/(xl(3)-xl(2))
      Ci = 1./3. * (xl(4)*xl(4)*xl(4)-xl(3)*xl(3)*xl(3))/(xl(4)-xl(3))

      a = ((x(3)-x(2))*(y(2)-y(1))) - ((x(2)-x(1))*(y(3)-y(2)))
      a = a/( (x(3)-x(2))*(Bi-Ai) - (x(2)-x(1))*(Ci-Bi) )

      b = ( (Ci-Bi)*(y(2)-y(1)) - (Bi-Ai)*(y(3)-y(2)) )
      b = b / ( (Ci-Bi)*(x(2)-x(1)) - (Bi-Ai)*(x(3)-x(2)) )

      c = y(1) - Ai*a - x(1)*b

      yint = a*xint*xint + b*xint + c
      return
      end subroutine ut_quadraticCellAverageInterpol


