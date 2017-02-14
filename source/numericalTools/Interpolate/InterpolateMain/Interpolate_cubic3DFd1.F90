!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic3DFd1
!!
!! NAME
!!
!!  Interpolate_cubic3DFd1
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic3DFd1 (real, intent (in) :: a (1:64),
!!                          real, intent (in) :: x,
!!                          real, intent (in) :: y,
!!                          real, intent (in) :: z)
!!
!! DESCRIPTION
!!
!!  Calculates the function value and all rescaled 1st derivative values for a triple
!!  [x,y,z] of rescaled [0,1] coordinates and the 64 tricubic expansion coefficients.
!!  The tricubic expansion reads, for one cube, in terms of rescaled [0,1] x,y,z coordinates:
!!
!!                                   3   3   3              i j k
!!                      F (x,y,z) = sum sum sum  a (i,j,k) x y z
!!                                  i=0 j=0 k=0
!!
!!  The order of the supplied expansion coefficients a (i,j,k) must be such, that the
!!  k index has the highest ranking, followed by the j index and the i index. The overall
!!  location index of the a (i,j,k) inside the 64-dimensional vector is given by the
!!  following formula:
!!
!!                location index of (i,j,k)  =  1 + i + 4j + 16k
!!
!!  The rescaled derivatives are given by the general formula:
!!
!!         r+s+t    r  s  t    3   3   3                                 i-r j-s k-t
!!        d     / dx dy dz  = sum sum sum (i) * (j) * (k) * a (i,j,k) * x   y   z
!!                            i=r j=s k=t    r     s     t
!!
!!  where the Pochhammer symbols are defined as:
!!
!!                  (i)  = i * (i-1) * (i-2) * ... * (i-r+1)
!!                     r
!!
!!  The rescaled derivatives are therefore also sums of appropriate expansion coefficients
!!  times monomial products. From the general formula we see, that the highest non-zero
!!  derivative is of 9-th order.
!!
!!  Since this function is (potentially) called many times from external applications,
!!  efficiency is key here and intermediate common summation terms are reused as much as
!!  possible. The strategy is partial summation and reduction at each index summation stage.
!!  The individual x-,y- and z-coordinate cubic polynomial sections are always evaluated
!!  using the Horner scheme to minimize accumulation of computation rounding errors.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th tricubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!  y     : rescaled [0,1] y coordinate
!!  z     : rescaled [0,1] z coordinate
!!
!! NOTES
!!
!!  1) The function is defined as a real array of size 4:
!!
!!           Interpolate_cubic3DFd1 (1) = the function value
!!           Interpolate_cubic3DFd1 (2) = the rescaled d/dx value
!!           Interpolate_cubic3DFd1 (3) = the rescaled d/dy value
!!           Interpolate_cubic3DFd1 (4) = the rescaled d/dz value
!!
!!  2) The code checks, if the supplied triple [x,y,z] is rescaled.
!!
!!***

function Interpolate_cubic3DFd1 (a,x,y,z)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:64)
  real, intent (in) :: x,y,z

  logical :: notRescaled

  real    :: x1p5, x2
  real    :: y1p5, y2
  real    :: z1p5, z2

  real    :: Interpolate_cubic3DFd1 (1:4)  ! declares the function as an array

  real    :: b (1:16)                      ! will hold partial i-index  summation values
  real    :: c (1:4)                       ! will hold partial ij-index summation values
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. (z < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0) .or. (z > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [Interpolate_cubic3DFd1]'        )
      call Logfile_stamp     (y, ' = rescaled y coordinate [Interpolate_cubic3DFd1]'        )
      call Logfile_stamp     (z, ' = rescaled z coordinate [Interpolate_cubic3DFd1]'        )
      call Driver_abortFlash ('[Interpolate_cubic3DFd1] ERROR: [x,y,z] triple not rescaled!')
  end if
!
!
!     ...Calculate needed multiples of the rescaled coordinates.
!
!
  x1p5 = x + (x / 2)
  y1p5 = y + (y / 2)
  z1p5 = z + (z / 2)
  x2   = x + x
  y2   = y + y
  z2   = z + z
!
!
!     ...Generate the function and derivative values. Intermediate vectors are
!        used to collect common summation contributions. The Horner scheme is
!        applied to calculate the individual polynomial contributions. The
!        multiplication count for each step is listed for information.
!
!                                                                                  ! used for   * count
!
  b (1:16) = ((a (4:64:4) * x + a (3:63:4)) * x + a (2:62:4)) * x + a (1:61:4)     ! f,dy,dz       48
  c (1:4)  = ((b (4:16:4) * y + b (3:15:4)) * y + b (2:14:4)) * y + b (1:13:4)     ! f,dz          12

  Interpolate_cubic3DFd1 (1) = ((c (4) * z    + c (3)) * z  + c (2)) * z + c (1)   ! -> f           3
  Interpolate_cubic3DFd1 (4) =  (c (4) * z1p5 + c (3)) * z2 + c (2)                ! -> dz          2
  
  c (1:4) = (b (4:16:4) * y1p5 + b (3:15:4)) * y2 + b (2:14:4)                     ! dy             8

  Interpolate_cubic3DFd1 (3)  = ((c (4) * z    + c (3)) * z  + c (2)) * z + c (1)  ! -> dy          3

  b (1:16) =  (a (4:64:4) * x1p5 + a (3:63:4)) * x2 + a (2:62:4)                   ! dx            32
  c (1:4)  = ((b (4:16:4) * y    + b (3:15:4)) * y  + b (2:14:4)) * y + b (1:13:4) ! dx            12

  Interpolate_cubic3DFd1 (2) = ((c (4) * z    + c (3)) * z  + c (2)) * z + c (1)   ! -> dx          3
!                                                                                                 ---
!                                                                                       total =   123
!
!     ...Ready!
!
!
  return
end function Interpolate_cubic3DFd1
