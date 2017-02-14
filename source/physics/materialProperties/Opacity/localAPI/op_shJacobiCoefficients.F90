!!****if* source/physics/materialProperties/Opacity/localAPI/op_shJacobiCoefficients
!!
!! NAME
!!
!!  op_shJacobiCoefficients
!!
!! SYNOPSIS
!!
!!  call op_shJacobiCoefficients (integer (in)  :: n,
!!                                real    (in)  :: p,
!!                                real    (in)  :: q,
!!                                real    (out) :: JacA,
!!                                real    (out) :: JacB)
!!
!! DESCRIPTION
!!
!!  Generates the A and B coefficients defining the 3-term recursion relation for the
!!  shifted Jacobi (monic) polynomials:
!!
!!             G (p,q,x) =  (x - A ) G   (p,q,x)  -  B  G   (p,q,x)
!!              i                 i   i-1             i  i-2
!!
!!  We have:
!!
!!                  2(i-1)(i+p-1) + q(p-1)
!!          A   =  ------------------------              i = 1,2,...,n
!!           i         (2i+p-1)(2i+p-3)
!!
!!
!!                 (i-1)(i+p-2)(i+q-2)(i+p-q-1)
!!          B   =  ----------------------------          i = 2,3,...,n
!!           i      (2i+p-2)(2i+p-4)(2i+p-3)^2
!!
!! ARGUMENTS
!!
!!  n      : the maximum order of the polynomial wanted
!!  p      : the p value in the defining shifted Jacobi polynomial G(p,q,x)
!!  q      : the q value in the defining shifted Jacobi polynomial G(p,q,x)
!!  JacA   : the polynomial A coefficients
!!  JacB   : the polynomial B coefficients
!!
!!***
subroutine op_shJacobiCoefficients (n,p,q,JacA,JacB)

  implicit none

  integer, intent (in) :: n
  real,    intent (in) :: p,q

  real,    intent (out), dimension (1:n) :: JacA
  real,    intent (out), dimension (1:n) :: JacB

  JacA = 0.0
  JacB = 0.0

  return
end subroutine op_shJacobiCoefficients
