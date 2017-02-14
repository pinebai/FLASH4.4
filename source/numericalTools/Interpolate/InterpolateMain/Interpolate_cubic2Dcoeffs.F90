!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic2Dcoeffs
!!
!! NAME
!!
!!  Interpolate_cubic2Dcoeffs
!!
!! SYNOPSIS
!!
!!  call Interpolate_cubic2Dcoeffs (integer, intent (in)    :: numberOfSquares,
!!                                  real,    intent (inout) :: a (1:16,*))
!!
!! DESCRIPTION
!!
!!  Calculates the bicubic expansion coefficients for a collection of squares. The bicubic
!!  expansion reads, for one square, in terms of rescaled [0,1] x,y coordinates:
!!
!!                                   3   3            i j
!!                        C (x,y) = sum sum  a (i,j) x y
!!                                  i=0 j=0
!!
!!
!!  and is uniquely determined by specifying the following 4 values at each of the 4 corners
!!  of the square (16 constraints for 16 unknown expansion coefficients):
!!
!!                             C (x,y)        (value of function at corner)
!!                              d/dx          (1st rescaled x-derivative)
!!                              d/dy          (1st rescaled y-derivative)
!!                             d2/dxdy        (2nd mixed rescaled xy-derivative)
!!
!!
!!  The square itself is defined by its 0 and 1 coordinates of its 4 corners:
!!
!!
!!                                  y
!!
!!                                  |
!!                                  |
!!                                  |
!!                                1  -----------
!!                                  |           |
!!                                  |           |
!!                                  |           |
!!                                  |           |
!!                                  |           |
!!                                0  ----------- ---- x
!!                                  0           1
!!
!!
!!  In order to obtain the function and global derivative values at a global position (X,Y)
!!  inside the square, one must first form the rescaled x,y coordinates:
!!
!!                                   x = (X - X0) / (X1 - X0)
!!                                   y = (Y - Y0) / (Y1 - Y0)
!!
!!  where X0,Y0 and X1,Y1 are the lower and upper global square x,y coordinates:
!!
!!
!!                                  Y
!!
!!                                  |
!!                                  |
!!                                  |
!!                               Y1  -----------
!!                                  |           |
!!                                  |           |
!!                                  |           |
!!                                  |           |
!!                                  |           |
!!                               Y0  ----------- ---- X
!!                                  X0         X1
!!
!!
!!  The function and global derivative (using the chain rule) values are then:
!!
!!
!!                              3   3            i j
!!                   C (X,Y) = sum sum  a (i,j) x y
!!                             i=0 j=0
!!
!!                              3   3                i-1 j
!!                      d/dX = sum sum  i * a (i,j) x   y  * dx/dX
!!                             i=1 j=0
!!
!!                              3   3                i j-1
!!                      d/dY = sum sum  j * a (i,j) x y    * dy/dY
!!                             i=0 j=1
!!
!!                              3   3                    i-1 j-1
!!                    d/dXdY = sum sum  i * j * a (i,j) x   y    * dx/dX * dy/dY
!!                             i=1 j=1
!!
!!
!!  where the rescaled to global coordinate differentials are:
!!
!!
!!                             dx/dX = 1 / (X1 - X0)
!!                             dy/dY = 1 / (Y1 - Y0)
!!
!!
!!  i.e. the inverses of the corresponding global square dimensions.
!!
!!
!!  Order of Input
!!  --------------
!!
!!  The order of the input values (function + derivatives) must be such, that each function/derivative
!!  must have its corner values clustered together in the order shown below. The function + derivatives
!!  order must follow the order mentioned above. We thus have the following ordering scheme:
!!  
!!
!!                                            Corners
!!
!!                                              x y
!!
!!                               C (x,y)    /   0 0
!!                                d/dx     /    1 0
!!                                d/dy    <     0 1
!!                               d2/dxdy   \    1 1
!!
!!  Order of Output
!!  ---------------
!!
!!  The order of the output coefficients a (i,j) is such that the j index has the highest ranking,
!!  followed by the i index. The overall location index of the a (i,j) inside the 16-dimensional
!!  vector is given by the following formula:
!!
!!
!!                  location index of (i,j)  =  1 + i + 4j
!!
!!
!! ARGUMENTS
!!
!!  numberOfSquares : the number of squares to be treated
!!  a (i,j)         : i-th function/derivative (input) and expansion coefficient (output) of j-th square
!!
!! NOTES
!!
!!  The array holding initially the function + derivative values and later the bicubic
!!  expansion coefficients is passed as an assumed size array (using *). This allows for
!!  compact looping over all cubes, even if in the calling routine this array is of different
!!  shape. The only drawback of this is that array operations on the assumed size array
!!  cannot be performed, i.e. every array element must be addressed by specific indices
!!  (like a (i,j) for example), which is the case here. Operations like 'a(1,:) = 1.0'
!!  or 'size (a,2)' cannot be done!
!!
!!  The code is written entirely without a single multiplication operation. Only
!!  the addition and substraction operation is used. This results in an optimum
!!  performance.
!!
!!  The code allows for threading to be used on the number of squares loop.
!!
!!***

subroutine Interpolate_cubic2Dcoeffs (numberOfSquares, a)

  implicit none

  integer, intent (in)    :: numberOfSquares
  real,    intent (inout) :: a (1:16,*)

  integer :: n

  real    :: p1, p2, p3, p4
  real    :: q1, q2, q3, q4
  real    :: r1, r2, r3, r4
  real    :: s1, s2, s3, s4

  real    :: x1,x2,x3,x4,x5,x6
!
!
!     ...Loop over all squares.
!
!
!$omp parallel do schedule (static)
  do n = 1, numberOfSquares
!
!
!     ...Collect intermediate values:
!
!           1st contribution from a  (1) to a  (4) (function   values at the 4 corners)
!           2nd contribution from a  (5) to a  (8) (  d/dx     values at the 4 corners)
!           3rd contribution from a  (8) to a (12) (  d/dy     values at the 4 corners)
!           4th contribution from a (12) to a (16) ( d2/dxdy   values at the 4 corners)
!
!
!    p1 = a ( 1,n)                       ! not needed, shown for informational purpose only
     q1 = a ( 1,n) - a ( 2,n)
     r1 = a ( 1,n) - a ( 3,n)
     s1 = a ( 4,n) - a ( 3,n)

     p2 = a ( 5,n)
     q2 = a ( 5,n) + a ( 6,n)
     r2 = a ( 5,n) - a ( 7,n)
     s2 = a ( 6,n) - a ( 8,n)

     p3 = a ( 9,n)
     q3 = a ( 9,n) + a (11,n)
     r3 = a ( 9,n) - a (10,n)
     s3 = a (11,n) - a (12,n)

     p4 = a (13,n)
     q4 = a (13,n) + a (14,n)
     r4 = a (13,n) + a (15,n)
     s4 = a (15,n) + a (16,n)
!
!
!     ...All function contributions to expansion coefficients.
!
!
     x1  = q1 + q1                       ! 2q          (2  terms)
     x2  = q1 + s1                       ! q + s       (2  terms)
     x3  = r1 + r1                       ! 2r          (2  terms)
     x4  = x2 + x2                       ! 2q + 2s     (4  terms)
     x5  = x4 + x2                       ! 3q + 3s     (6  terms)
     x6  = x5 + x5                       ! 6q + 6s     (12 terms)

!    a ( 1,n) = p1                       ! not needed, shown for informational purpose only
     a ( 3,n) = - x1 - q1
     a ( 4,n) = x1
     a ( 9,n) = - x3 - r1
     a (11,n) = x6 + x5
     a (12,n) = - x6
     a (13,n) = x3
     a (15,n) = - x6
     a (16,n) = x4 + x4
!
!
!     ...All d/dx contributions to expansion coefficients.
!
!
     x1  = r2 + r2                       ! r + r       (2  terms)
     x2  = r2 + s2                       ! r + s       (2  terms)
     x3  = x2 + r2                       ! 2r + s      (3  terms)
     x4  = x2 + x2                       ! 2r + 2s     (4  terms)
     x5  = x3 + x3                       ! 4r + 2s     (6  terms)

     a ( 2,n) = p2
     a ( 3,n) = a (3,n) - p2 - q2
     a ( 4,n) = a (4,n) + q2
     a (10,n) = - x1 - r2
     a (14,n) = x1
     a (11,n) = a (11,n) + x3 + x5
     a (12,n) = a (12,n) - x2 - x4
     a (15,n) = a (15,n) - x5
     a (16,n) = a (16,n) + x4
!
!
!     ...All d/dy contributions to expansion coefficients.
!
!
     x1  = r3 + r3                       ! r + r       (2  terms)
     x2  = r3 + s3                       ! r + s       (2  terms)
     x3  = x2 + r3                       ! 2r + s      (3  terms)
     x4  = x2 + x2                       ! 2r + 2s     (4  terms)
     x5  = x3 + x3                       ! 4r + 2s     (6  terms)

     a ( 5,n) = p3
     a ( 7,n) = - x1 - r3
     a ( 8,n) = x1
     a ( 9,n) = a ( 9,n) - p3 - q3
     a (11,n) = a (11,n) + x3 + x5
     a (12,n) = a (12,n) - x5
     a (13,n) = a (13,n) + q3
     a (15,n) = a (15,n) - x2 - x4
     a (16,n) = a (16,n) + x4
!
!
!     ...All d2/dxdy contributions to expansion coefficients.
!
!
     x1  = p4 + q4                       ! p + q       (2  terms)
     x2  = q4 + s4                       ! q + s       (2  terms)

     a ( 6,n) = p4
     a ( 7,n) = a ( 7,n) - x1
     a ( 8,n) = a ( 8,n) + q4
     a (10,n) = a (10,n) - p4 - r4
     a (11,n) = a (11,n) + x2 + x1 + r4
     a (12,n) = a (12,n) - x2 - q4
     a (14,n) = a (14,n) + r4
     a (15,n) = a (15,n) - x2 - r4
     a (16,n) = a (16,n) + x2
!
!
!     ...Next square.
!
!
  end do
!$omp end parallel do
!
!
!     ...Ready!
!
!
  return
end subroutine Interpolate_cubic2Dcoeffs
