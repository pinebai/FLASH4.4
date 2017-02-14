!!****if* source/physics/materialProperties/Opacity/localAPI/op_generateRootsWeights
!!
!! NAME
!!
!!  op_generateRootsWeights
!!
!! SYNOPSIS
!!
!!  call op_generateRootsWeights (integer (in)    :: n,
!!                                integer (in)    :: nV,
!!                                real    (in)    :: MomZero,
!!                                real    (inout) :: V       (1:nV),
!!                                real    (inout) :: Diag    (1:n),
!!                                real    (inout) :: Offd    (1:n),
!!                                real    (out)   :: Roots   (1:n),
!!                                real    (out)   :: Weights (1:n))
!!
!! DESCRIPTION
!!
!!  Computes the quadrature roots and weights from a previously established symmetric
!!  terminal matrix by the Golub-Welsch algorithm (see G.H. Golub and J.H. Welsch,
!!  Math. of Computation 23, p. 221-230 and A1-A10, 1969), which is based on a result
!!  of Wilf (see H.S. Wilf, Mathematics for the Physical Sciences, New York: Wiley,
!!  Problem 9, p. 80). Wilf has shown that if V (k,i) is the k-th element of the i-th
!!  normalized eigenvector of the terminal matrix corresponding to the i-th eigenvalue
!!  e (i), then the roots (i.e. the zeros of the n-th orthogonal monic polynomial) and
!!  weights to be used for the quadrature are given by:
!!
!!                             Roots (i) = e (I)
!!                           Weights (i) = MomZero * (V (1,i)**2)
!!
!!  where MomZero is the value of the definite integral over the weight function W (x)
!!  alone. The routine performs hence a diagonalization of the tridiagonal symmetric
!!  terminal matrix keeping only the first component of the eigenvectors and sets the
!!  roots and weights equal to the above relations. The diagonalization code uses the
!!  implicit QL method. The original diagonals and offdiagonals of the terminal matrix
!!  are destroyed during the diagonalization process.
!!
!! ARGUMENTS
!!
!!  n       : the size of the terminal matrix
!!  nV      : the declared dimension for the first eigenvector component vector
!!  MomZero : the 0-th moment over the weight function
!!  V       : vector that will contain the first components of the eigenvectors
!!  Diag    : the n diagonals of the tridiagonal terminal matrix
!!  Offd    : the n-1 offdiagonals of the tridiagonal terminal matrix (located at 1:n-1)
!!  Roots   : the computed roots
!!  Weights : the computed weights
!!
!!***
subroutine op_generateRootsWeights (n, nV, MomZero, V, Diag, Offd, Roots, Weights)

  implicit none

  integer, intent (in)    :: n
  integer, intent (in)    :: nV
  real,    intent (in)    :: MomZero

  real,    intent (inout) :: V       (1:nV)
  real,    intent (inout) :: Diag    (1:n)
  real,    intent (inout) :: Offd    (1:n)
  real,    intent (out)   :: Roots   (1:n)
  real,    intent (out)   :: Weights (1:n)

  V       = 0.0
  Diag    = 0.0
  Offd    = 0.0
  Roots   = 0.0
  Weights = 0.0

  return
end subroutine op_generateRootsWeights
