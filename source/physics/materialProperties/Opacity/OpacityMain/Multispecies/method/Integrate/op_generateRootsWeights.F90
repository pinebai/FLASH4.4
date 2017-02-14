!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_generateRootsWeights
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

  use Driver_interface, ONLY : Driver_abortFlash
  use op_numericsData,  ONLY : zero,one,two

  implicit none

  integer, intent (in)    :: n
  integer, intent (in)    :: nV
  real,    intent (in)    :: MomZero

  real,    intent (inout) :: V       (1:nV)
  real,    intent (inout) :: Diag    (1:n)
  real,    intent (inout) :: Offd    (1:n)
  real,    intent (out)   :: Roots   (1:n)
  real,    intent (out)   :: Weights (1:n)

  integer :: i,j,m
  integer :: iter

  real    :: c,d,f,g,p,r,s
  real    :: test
  real    :: phytag

  integer, parameter :: maxIter = 30
!
!
!   ...Handle the special case of a 1x1 terminal matrix.
!
!
  select case (n)

    case (1)

      Roots   (1) = Diag (1)
      Weights (1) = MomZero
!
!
!   ...The general case of matrices with dimensions > 1x1. Check first the size
!      of the supplied first eigenvector component vector.
!
!
    case default

      if (nV < n) then
          call Driver_abortFlash ('[op_generateRootsWeights] ERROR: Eigenvector V too small')
      end if

      V (1)    = one
      V (2:n)  = zero
      Offd (n) = zero
!
!
!   ...QL iterates on a tridiagonal matrix.
!
!
      do j = 1,n

         iter = 0

 1000    do m = j,n
            if (m == n) exit
            test = abs (Diag (m)) + abs (Diag (m+1))
            if (abs (Offd (m)) <= epsilon (one) * test) exit
         end do

         if (m /= j) then

             if (iter == maxIter) then
                 call Driver_abortFlash ('[op_generateRootsWeights] ERROR: Too many iterations in QL method')
             end if

             iter = iter + 1

             g = (Diag (j+1) - Diag (j)) / (two * Offd (j))
             r = phytag (g,one)                                          ! r = sqrt (g * g + one)
             g = Diag (m) - Diag (j) + Offd (j) / (g + sign (r,g))
             s = one
             c = one
             p = zero

             do i = m-1,j,-1
                f = s * Offd (i)
                d = c * Offd (i)
                r = phytag (f,g)                                         ! r = sqrt (f * f + g * g)
                Offd (i+1) = r

                if (r == zero) then
                    Diag (i+1) = Diag (i+1) - p
                    Offd (m) = zero
                    goto 1000
                end if

                s = f / r
                c = g / r
                g = Diag (i+1) - p
                r = (Diag (i) - g) * s + two * c * d
                p = s * r
                Diag (i+1) = g + p
                g = c * r - d

                f = V (i+1)
                V (i+1) = s * V (i) + c * f
                V (i)   = c * V (i) - s * f

             end do

             Diag (j) = Diag (j) - p
             Offd (j) = g
             Offd (m) = zero

             goto 1000

         end if
      end do
!
!
!   ...Calculate the roots and weights.
!
!
      do i = 1,n
         Roots   (i) = Diag (i)
         Weights (i) = MomZero * V (i) * V (i)
      end do

  end select
!
!
!   ...Ready! 
!
!
  return
end subroutine op_generateRootsWeights

real function phytag (a,b)
implicit none
real, intent (in) :: a,b
real  :: absa,absb
real, parameter :: one  = 1.0
real, parameter :: zero = 0.0
absa = abs (a)
absb = abs (b)
if (absa > absb) then
    phytag = absa * sqrt (one + (absb/absa)**2)
else
    if (absb == zero) then
        phytag = zero
    else
        phytag = absb * sqrt (one + (absa/absb)**2)
    end if
end if
return
end
