!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic3Dcoeffs
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

  integer :: n

  real    :: p1, p2, p3, p4, p5, p6, p7, p8
  real    :: q1, q2, q3, q4, q5, q6, q7, q8
  real    :: r1, r2, r3, r4, r5, r6, r7, r8
  real    :: s1, s2, s3, s4, s5, s6, s7, s8
  real    :: t1, t2, t3, t4, t5, t6, t7, t8
  real    :: u1, u2, u3, u4, u5, u6, u7, u8
  real    :: v1, v2, v3, v4, v5, v6, v7, v8
  real    :: w1, w2, w3, w4, w5, w6, w7, w8

  real    :: x1 ,x2 ,x3 ,x4 ,x5 ,x6 ,x7 ,x8 ,x9 ,x10, &
             x11,x12,x13,x14,x15,x16,x17,x18,x19,x20, &
             x21,x22
!
!
!     ...Loop over all cubes.
!
!
!$omp parallel do schedule (static)
  do n = 1, numberOfCubes
!
!
!     ...Collect intermediate values:
!
!           1st contribution from a  (1) to a  (8) (function   values at the 8 corners)
!           2nd contribution from a  (9) to a (16) (  d/dx     values at the 8 corners)
!           3rd contribution from a (17) to a (24) (  d/dy     values at the 8 corners)
!           4th contribution from a (25) to a (32) (  d/dz     values at the 8 corners)
!           5th contribution from a (33) to a (40) ( d2/dxdy   values at the 8 corners)
!           6th contribution from a (41) to a (48) ( d2/dxdz   values at the 8 corners)
!           7th contribution from a (49) to a (56) ( d2/dydz   values at the 8 corners)
!           8th contribution from a (57) to a (64) ( d3/dxdydz values at the 8 corners)
!
!
     p1 =   a ( 1,n)
     q1 =   a ( 1,n) - a ( 2,n)
     r1 =   a ( 1,n) - a ( 3,n)
     s1 =   a ( 1,n) - a ( 4,n)
     t1 =   a ( 6,n) - a ( 2,n)
     u1 =   a ( 5,n) - a ( 3,n)
     v1 =   a ( 7,n) - a ( 4,n)
     w1 = - a ( 8,n)

     p2 =   a ( 9,n)
     q2 =   a ( 9,n) + a (10,n)
     r2 =   a ( 9,n) - a (11,n)
     s2 =   a ( 9,n) - a (12,n)
     t2 =   a (10,n) - a (14,n)
     u2 = - a (11,n) - a (13,n)
     v2 = - a (12,n) + a (15,n)
     w2 =   a (16,n)

     p3 =   a (17,n)
     q3 =   a (17,n) + a (19,n)
     r3 =   a (17,n) - a (20,n)
     s3 =   a (17,n) - a (18,n)
     t3 =   a (19,n) - a (21,n)
     u3 = - a (20,n) - a (23,n)
     v3 = - a (18,n) + a (22,n)
     w3 =   a (24,n)

     p4 =   a (25,n)
     q4 =   a (25,n) + a (28,n)
     r4 =   a (25,n) - a (27,n)
     s4 =   a (25,n) - a (26,n)
     t4 =   a (28,n) - a (30,n)
     u4 = - a (27,n) - a (31,n)
     v4 = - a (26,n) + a (29,n)
     w4 =   a (32,n)

     p5 =   a (33,n)
     q5 =   a (33,n) + a (34,n)
     r5 =   a (33,n) + a (35,n)
     s5 =   a (33,n) - a (36,n)
     t5 =   a (34,n) - a (38,n)
     u5 =   a (35,n) + a (37,n)
     v5 = - a (36,n) - a (39,n)
     w5 =   a (37,n) - a (40,n)

     p6 =   a (41,n)
     q6 =   a (41,n) + a (44,n)
     r6 =   a (41,n) + a (42,n)
     s6 =   a (41,n) - a (43,n)
     t6 =   a (44,n) - a (47,n)
     u6 =   a (42,n) + a (46,n)
     v6 = - a (43,n) - a (45,n)
     w6 =   a (46,n) - a (48,n)

     p7 =   a (49,n)
     q7 =   a (49,n) + a (52,n)
     r7 =   a (49,n) + a (51,n)
     s7 =   a (49,n) - a (50,n)
     t7 =   a (52,n) - a (54,n)
     u7 =   a (51,n) + a (55,n)
     v7 = - a (50,n) - a (53,n)
     w7 =   a (55,n) - a (56,n)

     p8 =   a (57,n)
     q8 =   a (57,n) + a (58,n)
     r8 =   a (57,n) + a (59,n)
     s8 =   a (57,n) + a (60,n)
     t8 =   a (58,n) + a (62,n)
     u8 =   a (59,n) + a (61,n)
     v8 =   a (60,n) + a (63,n)
     w8 =   a (64,n)
!
!
!     ...All function contributions to expansion coefficients.
!
!
     x1  = q1 + q1                                  ! 2q                               (2  terms)
     x2  = q1 + u1                                  ! q + u                            (2  terms)
     x3  = r1 + r1                                  ! 2r                               (2  terms)
     x4  = r1 + v1                                  ! r + v                            (2  terms)
     x5  = s1 + s1                                  ! 2s                               (2  terms)
     x6  = s1 + t1                                  ! s + t                            (2  terms)
     x7  = x2 + x2                                  ! 2q + 2u                          (4  terms)
     x8  = x7 + x2                                  ! 3q + 3u                          (6  terms)
     x9  = x8 + x8                                  ! 6q + 6u                          (12 terms)
     x10 = x4 + x4                                  ! 2r + 2v                          (4  terms)
     x11 = x10 + x4                                 ! 3r + 3v                          (6  terms)
     x12 = x11 + x11                                ! 6r + 6v                          (12 terms)
     x13 = x6 + x6                                  ! 2s + 2t                          (4  terms)
     x14 = x13 + x6                                 ! 3s + 3t                          (6  terms)
     x15 = x14 + x14                                ! 6s + 6t                          (12 terms)
     x16 = p1 + t1 + u1 + v1 + w1                   ! p + t + u + v + w                (5  terms)
     x17 = x16 + x16                                ! 2p + 2t + 2u + 2v + 2w           (10 terms)
     x18 = x17 + x16                                ! 3p + 3t + 3u + 3v + 3w           (15 terms)
     x19 = x17 + x17                                ! 4p + 4t + 4u + 4v + 4w           (20 terms)
     x20 = x18 + x18 + x18                          ! 9p + 9t + 9u + 9v + 9w           (45 terms)
     x21 = x20 + x18                                ! 12p + 12t + 12u + 12v + 12w      (60 terms)
     x22 = x20 + x20                                ! 18p + 18t + 18u + 18v + 18w      (90 terms)

!    a ( 1,n) = p1                                  ! not needed, shown for informational purpose only
     a ( 3,n) = - x1 - q1
     a ( 4,n) = x1
     a ( 9,n) = - x3 - r1
     a (11,n) = x9 + x8
     a (12,n) = - x9
     a (13,n) = x3
     a (15,n) = - x9
     a (16,n) = x7 + x7
     a (33,n) = - x5 - s1
     a (35,n) = x15 + x14
     a (36,n) = - x15
     a (41,n) = x12 + x11
     a (43,n) = - x22 - x20
     a (44,n) = x22
     a (45,n) = - x12
     a (47,n) = x22
     a (48,n) = - x21
     a (49,n) = x5
     a (51,n) = - x15
     a (52,n) = x13 + x13
     a (57,n) = - x12
     a (59,n) = x22
     a (60,n) = - x21
     a (61,n) = x10 + x10
     a (63,n) = - x21
     a (64,n) = x19 + x19
!
!
!     ...All d/dx contributions to expansion coefficients.
!
!
     x1  = q2 + u2                                 ! q + u                            (2  terms)
     x2  = r2 + r2                                 ! r + r                            (2  terms)
     x3  = r2 + v2                                 ! r + v                            (2  terms)
     x4  = s2 + s2                                 ! s + s                            (2  terms)
     x5  = s2 + t2                                 ! s + t                            (2  terms)
     x6  = x1 + r2                                 ! q + u + r                        (3  terms)
     x7  = x5 + s2                                 ! 2s + t                           (3  terms)
     x8  = x1 + x1                                 ! 2q + 2u                          (4  terms)
     x9  = x5 + x5                                 ! 2s + 2t                          (4  terms)
     x10 = x3 + x3                                 ! 2r + 2v                          (4  terms)
     x11 = p2 + t2 + u2 + v2 + w2                  ! p + t + u + v + w                (5  terms)
     x12 = x6 + x6                                 ! 2q + 2u + 2r                     (6  terms)
     x13 = x7 + x7                                 ! 4s + 2t                          (6  terms)
     x14 = x10 + x3                                ! 3r + 3v                          (6  terms)
     x15 = x11 + x3                                ! p + r + t + u + 2v + w           (7  terms)
     x16 = x11 + x11                               ! 2p + 2t + 2u + 2v + 2w           (10 terms)
     x17 = x14 + x14                               ! 6r + 6v                          (12 terms)
     x18 = x15 + x15                               ! 2p + 2r + 2t + 2u + 4v + 2w      (14 terms)
     x19 = x16 + x11                               ! 3p + 3t + 3u + 3v + 3w           (15 terms)
     x20 = x18 + x15                               ! 3p + 3r + 3t + 3u + 6v + 3w      (21 terms)
     x21 = x19 + x19                               ! 6p + 6t + 6u + 6v + 6w           (30 terms)
     x22 = x20 + x20                               ! 6p + 6r + 6t + 6u + 12v + 6w     (42 terms)

     a ( 2,n) = p2
     a ( 3,n) = a (3,n) - p2 - q2
     a ( 4,n) = a (4,n) + q2
     a (10,n) = - x2 - r2
     a (11,n) = a (11,n) + x6 + x12
     a (12,n) = a (12,n) - x1 - x8
     a (14,n) = x2
     a (15,n) = a (15,n) - x12
     a (16,n) = a (16,n) + x8
     a (34,n) = - x4 - s2
     a (35,n) = a (35,n) + x7 + x13
     a (36,n) = a (36,n) - x5 - x9
     a (42,n) = x17 + x14
     a (43,n) = a (43,n) - x22 - x20
     a (44,n) = a (44,n) + x21 + x19
     a (46,n) = - x17
     a (47,n) = a (47,n) + x22
     a (48,n) = a (48,n) - x21
     a (50,n) = x4
     a (51,n) = a (51,n) - x13
     a (52,n) = a (52,n) + x9
     a (58,n) = - x17
     a (59,n) = a (59,n) + x22
     a (60,n) = a (60,n) - x21
     a (62,n) = x10 + x10
     a (63,n) = a (63,n) - x18 - x18
     a (64,n) = a (64,n) + x16 + x16
!
!
!     ...All d/dy contributions to expansion coefficients.
!
!
     x1  = q3 + u3                                 ! q + u                            (2  terms)
     x2  = r3 + r3                                 ! r + r                            (2  terms)
     x3  = r3 + v3                                 ! r + v                            (2  terms)
     x4  = s3 + s3                                 ! s + s                            (2  terms)
     x5  = s3 + t3                                 ! s + t                            (2  terms)
     x6  = x1 + r3                                 ! q + u + r                        (3  terms)
     x7  = x5 + s3                                 ! 2s + t                           (3  terms)
     x8  = x1 + x1                                 ! 2q + 2u                          (4  terms)
     x9  = x5 + x5                                 ! 2s + 2t                          (4  terms)
     x10 = x3 + x3                                 ! 2r + 2v                          (4  terms)
     x11 = p3 + t3 + u3 + v3 + w3                  ! p + t + u + v + w                (5  terms)
     x12 = x6 + x6                                 ! 2q + 2u + 2r                     (6  terms)
     x13 = x7 + x7                                 ! 4s + 2t                          (6  terms)
     x14 = x10 + x3                                ! 3r + 3v                          (6  terms)
     x15 = x11 + x3                                ! p + r + t + u + 2v + w           (7  terms)
     x16 = x11 + x11                               ! 2p + 2t + 2u + 2v + 2w           (10 terms)
     x17 = x14 + x14                               ! 6r + 6v                          (12 terms)
     x18 = x15 + x15                               ! 2p + 2r + 2t + 2u + 4v + 2w      (14 terms)
     x19 = x16 + x11                               ! 3p + 3t + 3u + 3v + 3w           (15 terms)
     x20 = x18 + x15                               ! 3p + 3r + 3t + 3u + 6v + 3w      (21 terms)
     x21 = x19 + x19                               ! 6p + 6t + 6u + 6v + 6w           (30 terms)
     x22 = x20 + x20                               ! 6p + 6r + 6t + 6u + 12v + 6w     (42 terms)

     a ( 5,n) = p3
     a ( 7,n) = - x4 - s3
     a ( 8,n) = x4
     a ( 9,n) = a ( 9,n) - p3 - q3
     a (11,n) = a (11,n) + x7 + x13
     a (12,n) = a (12,n) - x13
     a (13,n) = a (13,n) + q3
     a (15,n) = a (15,n) - x5 - x9
     a (16,n) = a (16,n) + x9
     a (37,n) = - x2 - r3
     a (39,n) = x17 + x14
     a (40,n) = - x17
     a (41,n) = a (41,n) + x6 + x12
     a (43,n) = a (43,n) - x22 - x20
     a (44,n) = a (44,n) + x22
     a (45,n) = a (45,n) - x1 - x8
     a (47,n) = a (47,n) + x21 + x19
     a (48,n) = a (48,n) - x21
     a (53,n) = x2
     a (55,n) = - x17
     a (56,n) = x10 + x10
     a (57,n) = a (57,n) - x12
     a (59,n) = a (59,n) + x22
     a (60,n) = a (60,n) - x18 - x18
     a (61,n) = a (61,n) + x8
     a (63,n) = a (63,n) - x21
     a (64,n) = a (64,n) + x16 + x16
!
!
!     ...All d/dz contributions to expansion coefficients.
!
!
     x1  = q4 + u4                                 ! q + u                            (2  terms)
     x2  = r4 + r4                                 ! r + r                            (2  terms)
     x3  = r4 + v4                                 ! r + v                            (2  terms)
     x4  = s4 + s4                                 ! s + s                            (2  terms)
     x5  = s4 + t4                                 ! s + t                            (2  terms)
     x6  = x1 + r4                                 ! q + u + r                        (3  terms)
     x7  = x5 + s4                                 ! 2s + t                           (3  terms)
     x8  = x1 + x1                                 ! 2q + 2u                          (4  terms)
     x9  = x5 + x5                                 ! 2s + 2t                          (4  terms)
     x10 = x3 + x3                                 ! 2r + 2v                          (4  terms)
     x11 = p4 + t4 + u4 + v4 + w4                  ! p + t + u + v + w                (5  terms)
     x12 = x6 + x6                                 ! 2q + 2u + 2r                     (6  terms)
     x13 = x7 + x7                                 ! 4s + 2t                          (6  terms)
     x14 = x10 + x3                                ! 3r + 3v                          (6  terms)
     x15 = x11 + x3                                ! p + r + t + u + 2v + w           (7  terms)
     x16 = x11 + x11                               ! 2p + 2t + 2u + 2v + 2w           (10 terms)
     x17 = x14 + x14                               ! 6r + 6v                          (12 terms)
     x18 = x15 + x15                               ! 2p + 2r + 2t + 2u + 4v + 2w      (14 terms)
     x19 = x16 + x11                               ! 3p + 3t + 3u + 3v + 3w           (15 terms)
     x20 = x18 + x15                               ! 3p + 3r + 3t + 3u + 6v + 3w      (21 terms)
     x21 = x19 + x19                               ! 6p + 6t + 6u + 6v + 6w           (30 terms)
     x22 = x20 + x20                               ! 6p + 6r + 6t + 6u + 12v + 6w     (42 terms)

     a (17,n) = p4
     a (19,n) = - x4 - s4
     a (20,n) = x4
     a (25,n) = - x2 - r4
     a (27,n) = x17 + x14
     a (28,n) = - x17
     a (29,n) = x2
     a (31,n) = - x17
     a (32,n) = x10 + x10
     a (33,n) = a (33,n) - p4 - q4
     a (35,n) = a (35,n) + x7 + x13
     a (36,n) = a (36,n) - x13
     a (41,n) = a (41,n) + x6 + x12
     a (43,n) = a (43,n) - x22 - x20
     a (44,n) = a (44,n) + x22
     a (45,n) = a (45,n) - x12
     a (47,n) = a (47,n) + x22
     a (48,n) = a (48,n) - x18 - x18
     a (49,n) = a (49,n) + q4
     a (51,n) = a (51,n) - x5 - x9
     a (52,n) = a (52,n) + x9
     a (57,n) = a (57,n) - x1 - x8
     a (59,n) = a (59,n) + x21 + x19
     a (60,n) = a (60,n) - x21
     a (61,n) = a (61,n) + x8
     a (63,n) = a (63,n) - x21
     a (64,n) = a (64,n) + x16 + x16
!
!
!     ...All d2/dxdy contributions to expansion coefficients.
!
!
     x1  = p5 + q5                               ! p + q                            (2  terms)
     x2  = q5 + u5                               ! q + u                            (2  terms)
     x3  = r5 + v5                               ! r + v                            (2  terms)
     x4  = s5 + s5                               ! s + s                            (2  terms)
     x5  = s5 + t5                               ! s + t                            (2  terms)
     x6  = t5 + w5                               ! t + w                            (2  terms)
     x7  = x3 + s5                               ! r + s + v                        (3  terms)
     x8  = x5 + s5                               ! 2s + t                           (3  terms)
     x9  = x5 + x5                               ! 2s + 2t                          (4  terms)
     x10 = x3 + x3                               ! 2r + 2v                          (4  terms)
     x11 = x3 + x6                               ! r + t + v + w                    (4  terms)
     x12 = x7 + x7                               ! 2r + 2s + 2v                     (6  terms)
     x13 = x8 + x8                               ! 4s + 2t                          (6  terms)
     x14 = x11 + x3                              ! 2r + t + 2v + w                  (6  terms)
     x15 = x11 + x5                              ! r + s + 2t + v + w               (6  terms)
     x16 = x11 + x11                             ! 2r + 2t + 2v + 2w                (8  terms)
     x17 = x14 + x8                              ! 2r + 2s + 2t + 2v + w            (9  terms)
     x18 = x14 + x14                             ! 4r + 2t + 4v + 2w                (12 terms)
     x19 = x15 + x15                             ! 2r + 2s + 4t + 2v + 2w           (12 terms)
     x20 = x17 + x17                             ! 4r + 4s + 4t + 4v + 2w           (18 terms)

     a ( 6,n) = p5
     a ( 7,n) = a ( 7,n) - x1
     a ( 8,n) = a ( 8,n) + q5
     a (10,n) = a (10,n) - p5 - r5
     a (11,n) = a (11,n) + x2 + x1 + r5
     a (12,n) = a (12,n) - x2 - q5
     a (14,n) = a (14,n) + r5
     a (15,n) = a (15,n) - x2 - r5
     a (16,n) = a (16,n) + x2
     a (38,n) = - x4 - s5
     a (39,n) = a (39,n) + x8 + x13
     a (40,n) = a (40,n) - x5 - x9
     a (42,n) = a (42,n) + x7 + x12
     a (43,n) = a (43,n) - x17 - x20
     a (44,n) = a (44,n) + x15 + x19
     a (46,n) = a (46,n) - x3 - x10
     a (47,n) = a (47,n) + x14 + x18
     a (48,n) = a (48,n) - x11 - x16
     a (54,n) = x4
     a (55,n) = a (55,n) - x13
     a (56,n) = a (56,n) + x9
     a (58,n) = a (58,n) - x12
     a (59,n) = a (59,n) + x20
     a (60,n) = a (60,n) - x19
     a (62,n) = a (62,n) + x10
     a (63,n) = a (63,n) - x18
     a (64,n) = a (64,n) + x16
!
!
!     ...All d2/dxdz contributions to expansion coefficients.
!
!
     x1  = p6 + q6                               ! p + q                            (2  terms)
     x2  = q6 + u6                               ! q + u                            (2  terms)
     x3  = r6 + v6                               ! r + v                            (2  terms)
     x4  = s6 + s6                               ! s + s                            (2  terms)
     x5  = s6 + t6                               ! s + t                            (2  terms)
     x6  = t6 + w6                               ! t + w                            (2  terms)
     x7  = x3 + s6                               ! r + s + v                        (3  terms)
     x8  = x5 + s6                               ! 2s + t                           (3  terms)
     x9  = x5 + x5                               ! 2s + 2t                          (4  terms)
     x10 = x3 + x3                               ! 2r + 2v                          (4  terms)
     x11 = x3 + x6                               ! r + t + v + w                    (4  terms)
     x12 = x7 + x7                               ! 2r + 2s + 2v                     (6  terms)
     x13 = x8 + x8                               ! 4s + 2t                          (6  terms)
     x14 = x11 + x3                              ! 2r + t + 2v + w                  (6  terms)
     x15 = x11 + x5                              ! r + s + 2t + v + w               (6  terms)
     x16 = x11 + x11                             ! 2r + 2t + 2v + 2w                (8  terms)
     x17 = x14 + x8                              ! 2r + 2s + 2t + 2v + w            (9  terms)
     x18 = x14 + x14                             ! 4r + 2t + 4v + 2w                (12 terms)
     x19 = x15 + x15                             ! 2r + 2s + 4t + 2v + 2w           (12 terms)
     x20 = x17 + x17                             ! 4r + 4s + 4t + 4v + 2w           (18 terms)

     a (18,n) = p6
     a (19,n) = a (19,n) - p6 - r6
     a (20,n) = a (20,n) + r6
     a (26,n) = - x4 - s6
     a (27,n) = a (27,n) + x7 + x12
     a (28,n) = a (28,n) - x3 - x10
     a (30,n) = x4
     a (31,n) = a (31,n) - x12
     a (32,n) = a (32,n) + x10
     a (34,n) = a (34,n) - x1
     a (35,n) = a (35,n) + x2 + x1 + r6
     a (36,n) = a (36,n) - x2 - r6
     a (42,n) = a (42,n) + x8 + x13
     a (43,n) = a (43,n) - x17 - x20
     a (44,n) = a (44,n) + x14 + x18
     a (46,n) = a (46,n) - x13
     a (47,n) = a (47,n) + x20
     a (48,n) = a (48,n) - x18
     a (50,n) = a (50,n) + q6
     a (51,n) = a (51,n) - x2 - q6
     a (52,n) = a (52,n) + x2
     a (58,n) = a (58,n) - x5 - x9
     a (59,n) = a (59,n) + x15 + x19
     a (60,n) = a (60,n) - x11 - x16
     a (62,n) = a (62,n) + x9
     a (63,n) = a (63,n) - x19
     a (64,n) = a (64,n) + x16
!
!
!     ...All d2/dydz contributions to expansion coefficients.
!
!
     x1  = p7 + q7                               ! p + q                            (2  terms)
     x2  = q7 + u7                               ! q + u                            (2  terms)
     x3  = r7 + v7                               ! r + v                            (2  terms)
     x4  = s7 + s7                               ! s + s                            (2  terms)
     x5  = s7 + t7                               ! s + t                            (2  terms)
     x6  = t7 + w7                               ! t + w                            (2  terms)
     x7  = x3 + s7                               ! r + s + v                        (3  terms)
     x8  = x5 + s7                               ! 2s + t                           (3  terms)
     x9  = x5 + x5                               ! 2s + 2t                          (4  terms)
     x10 = x3 + x3                               ! 2r + 2v                          (4  terms)
     x11 = x3 + x6                               ! r + t + v + w                    (4  terms)
     x12 = x7 + x7                               ! 2r + 2s + 2v                     (6  terms)
     x13 = x8 + x8                               ! 4s + 2t                          (6  terms)
     x14 = x11 + x3                              ! 2r + t + 2v + w                  (6  terms)
     x15 = x11 + x5                              ! r + s + 2t + v + w               (6  terms)
     x16 = x11 + x11                             ! 2r + 2t + 2v + 2w                (8  terms)
     x17 = x14 + x8                              ! 2r + 2s + 2t + 2v + w            (9  terms)
     x18 = x14 + x14                             ! 4r + 2t + 4v + 2w                (12 terms)
     x19 = x15 + x15                             ! 2r + 2s + 4t + 2v + 2w           (12 terms)
     x20 = x17 + x17                             ! 4r + 4s + 4t + 4v + 2w           (18 terms)

     a (21,n) = p7
     a (23,n) = - x4 - s7
     a (24,n) = x4
     a (25,n) = a (25,n) - p7 - r7
     a (27,n) = a (27,n) + x7 + x12
     a (28,n) = a (28,n) - x12
     a (29,n) = a (29,n) + r7
     a (31,n) = a (31,n) - x3 - x10
     a (32,n) = a (32,n) + x10
     a (37,n) = a (37,n) - x1
     a (39,n) = a (39,n) + x8 + x13
     a (40,n) = a (40,n) - x13
     a (41,n) = a (41,n) + x2 + x1 + r7
     a (43,n) = a (43,n) - x17 - x20
     a (44,n) = a (44,n) + x20
     a (45,n) = a (45,n) - x2 - r7
     a (47,n) = a (47,n) + x14 + x18
     a (48,n) = a (48,n) - x18
     a (53,n) = a (53,n) + q7
     a (55,n) = a (55,n) - x5 - x9
     a (56,n) = a (56,n) + x9
     a (57,n) = a (57,n) - x2 - q7
     a (59,n) = a (59,n) + x15 + x19
     a (60,n) = a (60,n) - x19
     a (61,n) = a (61,n) + x2
     a (63,n) = a (63,n) - x11 - x16
     a (64,n) = a (64,n) + x16
!
!
!     ...All d3/dxdydz contributions to expansion coefficients.
!
!
     x1  = q8 + r8                                 ! q + r                            (2  terms)
     x2  = q8 + s8                                 ! q + s                            (2  terms)
     x3  = q8 + u8                                 ! q + u                            (2  terms)
     x4  = r8 + s8                                 ! r + s                            (2  terms)
     x5  = r8 + v8                                 ! r + v                            (2  terms)
     x6  = s8 + t8                                 ! s + t                            (2  terms)
     x7  = x3 + r8                                 ! q + r + u                        (3  terms)
     x8  = x6 + q8                                 ! q + s + t                        (3  terms)
     x9  = x5 + s8                                 ! r + s + v                        (3  terms)
     x10 = p8 + t8 + u8 + v8 + w8                  ! p + t + u + v + w                (5  terms)
     x11 = x10 + x3                                ! p + q + t + 2u + v + w           (7  terms)
     x12 = x10 + x5                                ! p + r + t + u + 2v + w           (7  terms)
     x13 = x10 + x6                                ! p + s + 2t + u + v + w           (7  terms)
     x14 = x10 + x10                               ! 2p + 2t + 2u + 2v + 2w           (10 terms)

     a (22,n) = p8
     a (23,n) = a (23,n) - p8 - q8
     a (24,n) = a (24,n) + q8
     a (26,n) = a (26,n) - p8 - r8
     a (27,n) = a (27,n) + x3 + x1 + p8
     a (28,n) = a (28,n) - x3 - q8
     a (30,n) = a (30,n) + r8
     a (31,n) = a (31,n) - x7
     a (32,n) = a (32,n) + x3
     a (38,n) = a (38,n) - p8 - s8
     a (39,n) = a (39,n) + x6 + x2 + p8
     a (40,n) = a (40,n) - x8
     a (42,n) = a (42,n) + x5 + x4 + p8
     a (43,n) = a (43,n) - x14 - x1 - x2 - x4 + w8
     a (44,n) = a (44,n) + x11 + x8
     a (46,n) = a (46,n) - x5 - r8
     a (47,n) = a (47,n) + x12 + x7
     a (48,n) = a (48,n) - x11
     a (54,n) = a (54,n) + s8
     a (55,n) = a (55,n) - x6 - s8
     a (56,n) = a (56,n) + x6
     a (58,n) = a (58,n) - x9
     a (59,n) = a (59,n) + x13 + x9
     a (60,n) = a (60,n) - x13
     a (62,n) = a (62,n) + x5
     a (63,n) = a (63,n) - x12
     a (64,n) = a (64,n) + x10
!
!
!     ...Next cube.
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
end subroutine Interpolate_cubic3Dcoeffs
