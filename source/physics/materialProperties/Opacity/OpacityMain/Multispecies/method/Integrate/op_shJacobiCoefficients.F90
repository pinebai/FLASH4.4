!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_shJacobiCoefficients
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

  use op_numericsData,  ONLY : one,two,three,four

  implicit none

  integer, intent (in)  :: n
  real,    intent (in)  :: p,q

  real,    intent (out) :: JacA (1:n)
  real,    intent (out) :: JacB (1:n)

  integer :: i
  real    :: x,xx
!
!
!   ...Calculate the coefficients.
!
!
  JacA (1) = q / (p + one)

  if (n == 1) return

  do i = 2,n
     x = real (i)
     xx = x + x
     JacA (i) = (       two*(x-one)*(x+p-one)+q*(p-one) ) / ( (xx+p-one)*(xx+p-three)                )
     JacB (i) = ( (x-one)*(x+p-two)*(x+q-two)*(x+p-q-1) ) / ( (xx+p-two)*(xx+p-four)*(xx+p-three)**2 )
  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_shJacobiCoefficients
