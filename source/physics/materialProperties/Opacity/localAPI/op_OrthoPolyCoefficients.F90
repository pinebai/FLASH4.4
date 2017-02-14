!!****if* source/physics/materialProperties/Opacity/localAPI/op_OrthoPolyCoefficients
!!
!! NAME
!!
!!  op_OrthoPolyCoefficients
!!
!! SYNOPSIS
!!
!!  call op_OrthoPolyCoefficients (integer (in)    :: n,
!!                                 integer (in)    :: m
!!                                 integer (in)    :: nRow,
!!                                 real    (in)    :: Mom,
!!                                 real    (in)    :: A,
!!                                 real    (in)    :: B,
!!                                 real    (inout) :: Row1,
!!                                 real    (inout) :: Row2,
!!                                 real    (out)   :: OrthoA,
!!                                 real    (out)   :: OrthoB)
!!
!! DESCRIPTION
!!
!!  Calculates the orthogonal polynomial 3-term recursion coefficients.
!!
!! ARGUMENTS
!!
!!  n        : the order of the polynomial
!!  m        : the number of modified moments
!!  nRow     : the declared dimension for the row vectors
!!  Mom      : the modified moments
!!  A        : the auxilliary polynomial A coefficients
!!  B        : the auxilliary polynomial B coefficients
!!  Row1     : row vector that will contain intermediate quantities
!!  Row2     : row vector that will contain intermediate quantities
!!  OrthoA   : the orthogonal polynomial A coefficients
!!  OrthoB   : the orthogonal polynomial B coefficients
!!
!!***
subroutine op_OrthoPolyCoefficients (n,m,nRow,Mom,A,B,Row1,Row2,OrthoA,OrthoB)

  implicit none

  integer, intent (in) :: n
  integer, intent (in) :: m
  integer, intent (in) :: nRow

  real,    intent (in),    dimension (1:m)    :: Mom
  real,    intent (in),    dimension (1:m)    :: A
  real,    intent (in),    dimension (1:m)    :: B
  real,    intent (inout), dimension (1:nRow) :: Row1
  real,    intent (inout), dimension (1:nRow) :: Row2
  real,    intent (out),   dimension (1:n)    :: OrthoA
  real,    intent (out),   dimension (1:n)    :: OrthoB

  Row1   = 0.0
  Row2   = 0.0
  OrthoA = 0.0
  OrthoB = 0.0

  return
end subroutine op_OrthoPolyCoefficients
