!!****f* source/numericalTools/Interpolate/Interpolate_cubic2Dcoeffs
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

  a (1:16,1:numberOfSquares) = 0.0

  return
end subroutine Interpolate_cubic2Dcoeffs
