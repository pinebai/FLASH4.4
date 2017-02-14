!!****f* source/numericalTools/Interpolate/Interpolate_cubic3Dcoeffs
!!
!! NAME
!!
!!  Interpolate_cubic3Dcoeffs
!!
!! SYNOPSIS
!!
!!  call Interpolate_cubic3Dcoeffs (integer, intent (in)    :: numberOfCubes,
!!                                  real,    intent (inout) :: a (1:64,*))
!!
!! DESCRIPTION
!!
!!  Calculates the tricubic expansion coefficients for a collection of cubes. The tricubic
!!  expansion reads, for one cube, in terms of rescaled [0,1] x,y,z coordinates:
!!
!!                                   3   3   3              i j k
!!                      C (x,y,z) = sum sum sum  a (i,j,k) x y z
!!                                  i=0 j=0 k=0
!!
!!
!!  and is uniquely determined by specifying the following 8 values at each of the 8 corners
!!  of the cube (64 constraints for 64 unknown expansion coefficients):
!!
!!                             C (x,y,z)            (value of function at corner)
!!                                d/dx              (1st rescaled x-derivative)
!!                                d/dy              (1st rescaled y-derivative)
!!                                d/dz              (1st rescaled z-derivative)
!!                               d2/dxdy            (2nd mixed rescaled xy-derivative)
!!                               d2/dxdz            (2nd mixed rescaled xz-derivative)
!!                               d2/dydz            (2nd mixed rescaled yz-derivative)
!!                               d3/dxdydz          (3rd mixed rescaled xyz-derivative)
!!
!!
!!  The cube itself is defined by its 0 and 1 coordinates of its 8 corners:
!!
!!
!!                             z
!!                             
!!                             |
!!                             |    -----------
!!                             |  /|    y     /|
!!                               / |   /     / |
!!                              /  |  /     /  |
!!                           1  -----------    |
!!                             |   |       |   |
!!                             | 1  -----------
!!                             |  /        |  /
!!                             | /         | /
!!                             |/          |/
!!                           0  -----------    ---- x
!!                             0           1
!!
!!
!!  In order to obtain the function and global derivative values at a global position (X,Y,Z)
!!  inside the cube, one must first form the rescaled x,y,z coordinates:
!!
!!                                   x = (X - X0) / (X1 - X0)
!!                                   y = (Y - Y0) / (Y1 - Y0)
!!                                   z = (Z - Z0) / (Z1 - Z0)
!!
!!  where X0,Y0,Z0 and X1,Y1,Z1 are the lower and upper global cube x,y,z coordinates:
!!
!!
!!                                 Z
!!   
!!                                 |
!!                                 |    -----------
!!                                 |  /|    Y     /|
!!                                   / |   /     / |
!!                                  /  |  /     /  |
!!                                  -----------    |
!!                                /|   |       |   |
!!                           Y1 ---|--  -----------
!!                              /  |  /        |  /
!!                             Z1  | /         | /
!!                                 |/          |/
!!                          Y0 ---  -----------    ---- X
!!                               / |           |
!!                             Z0  X0          X1
!!
!!
!!  The function and global derivative (using the chain rule) values are then:
!!
!!
!!                         3   3   3              i j k
!!            C (X,Y,Z) = sum sum sum  a (i,j,k) x y z
!!                        i=0 j=0 k=0
!!
!!                         3   3   3                  i-1 j k
!!                 d/dX = sum sum sum  i * a (i,j,k) x   y z  * dx/dX
!!                        i=1 j=0 k=0
!!
!!                         3   3   3                  i j-1 k
!!                 d/dY = sum sum sum  j * a (i,j,k) x y   z  * dy/dY
!!                        i=0 j=1 k=0
!!
!!                         3   3   3                  i j k-1
!!                 d/dZ = sum sum sum  k * a (i,j,k) x y z    * dz/dZ
!!                        i=0 j=0 k=1
!!
!!                         3   3   3                      i-1 j-1 k
!!              d2/dXdY = sum sum sum  i * j * a (i,j,k) x   y   z  * dx/dX * dy/dY
!!                        i=1 j=1 k=0
!!
!!                         3   3   3                      i-1 j k-1
!!              d2/dXdZ = sum sum sum  i * k * a (i,j,k) x   y z    * dx/dX * dz/dZ
!!                        i=1 j=0 k=1
!!
!!                         3   3   3                      i j-1 k-1
!!              d2/dYdZ = sum sum sum  j * k * a (i,j,k) x y   z    * dy/dY * dz/dZ
!!                        i=0 j=1 k=1
!!
!!                         3   3   3                          i-1 j-1 k-1
!!            d3/dXdYdZ = sum sum sum  i * j * k * a (i,j,k) x   y   z    * dx/dX * dy/dY * dz/dZ
!!                        i=1 j=1 k=1
!!
!!
!!  where the rescaled to global coordinate differentials are:
!!
!!
!!                             dx/dX = 1 / (X1 - X0)
!!                             dy/dY = 1 / (Y1 - Y0)
!!                             dz/dZ = 1 / (Z1 - Z0)
!!
!!
!!  i.e. the inverses of the corresponding global cube dimensions.
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
!!
!!                                           Corners
!!
!!                                            x y z
!!
!!                                        /   0 0 0
!!                          C (x,y,z)    /    1 0 0
!!                            d/dx      /     0 1 0
!!                            d/dy     /      0 0 1
!!                            d/dz    <       1 1 0
!!                           d2/dxdy   \      1 0 1
!!                           d2/dxdz    \     0 1 1
!!                           d2/dydz     \    1 1 1
!!                           d3/dxdydz
!!
!!  Order of Output
!!  ---------------
!!
!!  The order of the output coefficients a (i,j,k) is such that the k index has the highest
!!  ranking, followed by the j index and the i index. The overall location index of the
!!  a (i,j,k) inside the 64-dimensional vector is given by the following formula:
!!
!!
!!                location index of (i,j,k)  =  1 + i + 4j + 16k
!!
!!
!! ARGUMENTS
!!
!!  numberOfCubes : the number of cubes to be treated
!!  a (i,j)       : i-th function/derivative (input) and expansion coefficient (output) of j-th cube
!!
!! NOTES
!!
!!  For more information please refer to the following paper: F. Lekien and J. Marsden,
!!  International Journal for Numerical Methods in Engineering 63, p.455-471 (2005).
!!
!!  The array holding initially the function + derivative values and later the tricubic
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
!!  The code allows for threading to be used on the number of cubes loop.
!!
!!***

subroutine Interpolate_cubic3Dcoeffs (numberOfCubes, a)

  implicit none

  integer, intent (in)    :: numberOfCubes        ! this is essential info, needs to be passed
  real,    intent (inout) :: a (1:64,*)           ! assumed size array

  a (1:64,1:numberOfCubes) = 0.0

  return
end subroutine Interpolate_cubic3Dcoeffs
