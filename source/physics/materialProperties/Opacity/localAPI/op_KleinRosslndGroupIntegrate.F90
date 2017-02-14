!!****if* source/physics/materialProperties/Opacity/localAPI/op_KleinRosslndGroupIntegrate
!!
!! NAME
!!
!!  op_KleinRosslndGroupIntegrate
!!
!! SYNOPSIS
!!
!!  call op_KleinRosslndGroupIntegrate (real (in)  :: KNintegralFactor,
!!                                      real (in)  :: B,
!!                                      real (in)  :: R,
!!                                      real (in)  :: S,
!!                                      real (out) :: rescaleBase10Exponent,
!!                                      real (out) :: KleinRosslndIntegral,
!!                                      real (out) :: RosslndIntegral)
!!
!! DESCRIPTION
!!
!!  This routine evaluates the following two kinds of integrals using the Romberg integration technique:
!!
!!
!!                                                           S
!!                                                          /
!!                                                         |       -1
!!   Klein-Rosseland Integral =  (1 / KNintegralFactor) *  |  F[Bx]   * x^4 e^x / (e^x - 1)^2 dx
!!                                                         |
!!                                                        /
!!                                                       R
!!
!!      
!!                                  S
!!                                 /
!!                                |
!!         Rosseland Integral =   |   x^4 e^x / (e^x - 1)^2 dx
!!                                |
!!                               /
!!                              R
!!
!!
!!  where the function F[Bx] in the Klein-Planck Integral is defined as:
!!
!!
!!                     2+2[Bx]-[Bx]^2       1            2(1+[Bx])^2        (1+3[Bx])
!!             F[Bx] = -------------- ln ---------  +  ---------------  -  -----------
!!                        2[Bx]^3        1 + 2[Bx]     (1+2[Bx])[Bx]^2     (1+2[Bx])^2
!!
!!
!!  Several points need to be observed when evaluating these integrals via the Romberg discrete
!!  integration scheme. The problems come from the exponential function e^x present in both
!!  integrals, which for large x leads to overflow on the machine. Since the range of x can
!!  be very substantial, special measures have to be taken to get reasonable answers.
!!
!!  As the lower integration limits, say R, get larger, the exponential function in both integrals
!!  get closer to the point of computational overflow. At a certain R value it is then advantageous
!!  to introduce the rescaling factor e^R during the romberg scheme, such that we have representable
!!  exponential values for a wide range in x. The integrals are then multiplied by e^R and the
!!  resulting factor (e^R e^x) / (e^x - 1)^2 in the integrand is approximated as (e^R e^x) / (e^x e^x)
!!  = e^(R-x). Hencetemporary rescaling can only be applied when it is safe to ignore the -1 in
!!  (e^x - 1), meaning that this truncation does not lead to accuracy errors. It is hence dependent
!!  on the number of mantissa digits.
!!
!!  The upper integration limit, say S, might also be way beyond representability of e^S. The code has
!!  thus to decide upon a different upper integration limit U < S, beyond which all integrands are
!!  considered to be equal to zero. U is determined from the largest possible ebase exponent. If the
!!  condition U < S is found, S is replaced by U in the above integrals.
!!
!!  Evaluation of the above F[Bx] function is problematic for small values < 0.1 of Bx, requiring
!!  multiple precision accuracy for intermediate values to achieve a decent accuracy in the final
!!  result. It is thus best to approximate F[Bx] via a rational function. Such a function is:
!!
!!        F'[Bx] = P / Q
!!
!!  where P is a polynomial of order 4 in [Bx]:
!!
!!        P = (((y4*[Bx] + y3)*[Bx] + y2)*[Bx] + y1)*[Bx] + y0
!!
!!  and Q is a polynomial of order 5 in [Bx]:
!!
!!        Q = ((((z5*[Bx] + z4)*[Bx] + z3)*[Bx] + z2)*[Bx] + z1)*[Bx] + z0
!!
!!  For the values of {y0,y1,y2,y3,y4} and {z0,z1,z2,z3,z4,z5} given below as real parameters,
!!  the maximum relative error (F'-F)/F is < 1.E-9 for the range [Bx]:0-1.
!!
!!
!! ARGUMENTS
!!
!!  KNintegralFactor      : the Klein Nishina integral factor
!!  B                     : the multiplying factor for the integration variable x
!!  R                     : lower integration limit
!!  S                     : upper integration limit
!!  rescaleBase10Exponent : the value of E in the integral rescaling factor 10^E
!!  KleinRosslndIntegral  : the value of the Klein-Rosseland integral (rescaled by 10^E)
!!  RosslndIntegral       : the value of the Rosseland integral (rescaled by 10^E)
!!
!!***

subroutine op_KleinRosslndGroupIntegrate (KNintegralFactor,              &
                                          B,                             &
                                          R,S,                           &
                                                  rescaleBase10Exponent, &
                                                  KleinRosslndIntegral,  &
                                                  RosslndIntegral        )

  implicit none

  real, intent (in)  :: KNintegralFactor
  real, intent (in)  :: B
  real, intent (in)  :: R,S
  real, intent (out) :: rescaleBase10Exponent
  real, intent (out) :: KleinRosslndIntegral
  real, intent (out) :: RosslndIntegral

  rescaleBase10Exponent = 0.0
  KleinRosslndIntegral  = 0.0
  RosslndIntegral       = 0.0

  return
end subroutine op_KleinRosslndGroupIntegrate
