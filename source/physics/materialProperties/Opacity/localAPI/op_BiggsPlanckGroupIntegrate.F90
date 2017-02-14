!!****if* source/physics/materialProperties/Opacity/localAPI/op_BiggsPlanckGroupIntegrate
!!
!! NAME
!!
!!  op_BiggsPlanckGroupIntegrate
!!
!! SYNOPSIS
!!
!!  call op_BiggsPlanckGroupIntegrate (integer (in)  :: nSec,
!!                                     real    (in)  :: A1Sec     (1:nSec),
!!                                     real    (in)  :: A2Sec     (1:nSec),
!!                                     real    (in)  :: A3Sec     (1:nSec),
!!                                     real    (in)  :: A4Sec     (1:nSec),
!!                                     real    (in)  :: intLimits (1:nSec+1),
!!                                     real    (out) :: rescaleBase10Exponent,
!!                                     real    (out) :: BiggsPlanckIntegral,
!!                                     real    (out) :: PlanckIntegral)
!!
!! DESCRIPTION
!!
!!  This routine evaluates the following two kinds of integrals using the Romberg integration technique:
!!
!!
!!                               intLimits (n+1)
!!                                   /
!!                           nSec   |
!!   Biggs-Planck Integral =  Sum   |  (A1[n]/x + A2[n]/x^2 + A3[n]/x^3 + A4[n]/x^4) * x^3 / (e^x - 1) dx
!!                             n    |
!!                                 /
!!                             intLimits (n)
!!
!!
!!      
!!                           intLimits (nSec+1)
!!                               /
!!                              |
!!         Planck Integral =    |   x^3 / (e^x - 1) dx
!!                              |
!!                             /
!!                        intLimits (1)
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
!!  resulting factor e^R/(e^x - 1) in the integrand is approximated as e^R/e^x = e^(R-x). Hence
!!  temporary rescaling can only be applied when it is safe to ignore the -1 in (e^x - 1), meaning
!!  that this truncation does not lead to accuracy errors. It is hence dependent on the number of
!!  mantissa digits.
!!
!!  The upper integration limit, say S, might also be way beyond representability of e^S. The code has
!!  thus to decide upon a different upper integration limit U < S, beyond which all integrands are
!!  considered to be equal to zero. U is determined from the largest possible ebase exponent. If the
!!  condition U < S is found, S is replaced by U in the above integrals.
!!
!!  Note the difference between evaluation of both integrals: the Biggs-Planck integral must be
!!  evaluated piece by piece and the Planck integral is evaluated in one shot. The possible need for
!!  rescaling at certain integration limits leads to the following strategy in evaluating the
!!  Biggs-Planck integral:
!!
!!               i) Start from the lowest integration limit R = intLimits (1)
!!
!!              ii) Set the rescaling factor equal to 1/e^R (note that this factor
!!                  is the same for both integrals so when forming the mean value
!!                  they will cancel out and there is no need to store it). In case
!!                  no rescaling is needed this factor is 1. Evaluate the 1st section
!!                  of the Biggs-Planck integral and call it BP1. Then the true
!!                  value of the 1st section is:
!!
!!                             1st BP section = (1/e^R) * BP1
!!
!!             iii) Take the next lowest limit S = intLimits (2) for the Biggs-Planck
!!                  integral. If rescaling is necessary, the factor will be 1/e^S.
!!                  After evaluation of the second rescaled Biggs-Planck integral section
!!                  BP2, the true value of the second section is:
!!
!!                             2nd BP section = (1/e^S) * BP2
!!
!!                  and the total sum of both sections is:
!!
!!                       1st + 2nd BP section (not rescaled) = (1/e^R) * BP1  +  (1/e^S) * BP2
!!
!!                  However, we need to have this sum rescaled to the very first rescaling
!!                  factor, since this is the basic rescaling. Thus we have:
!!
!!                       1st + 2nd BP section (rescaled with 1/e^R) = BP1  +  (e^R/e^S) * BP2
!!                                                                  = BP1  +    e^(R-S) * BP2
!!
!!              iv) Since R < S, if (R - S) exceeds the lowest ebase exponent we can simply
!!                  ignore the evaluation of BP2 and thus all other higher sections. If
!!                  e^(R-S) is non-zero, we evaluate BP2 and add e^(R-S)*BP2 to BP1. Then
!!                  we proceed to the next (3rd) section and repeat the analysis.
!!
!!
!! ARGUMENTS
!!
!!  nSec                  : total number of integral sections
!!  A1Sec                 : 1st Biggs expansion coefficients for the x^(-1) function
!!  A2Sec                 : 2nd Biggs expansion coefficients for the x^(-2) function
!!  A3Sec                 : 3rd Biggs expansion coefficients for the x^(-3) function
!!  A4Sec                 : 4th Biggs expansion coefficients for the x^(-4) function
!!  intLimits             : integration limits for each section (lower in n, upper in n+1)
!!  rescaleBase10Exponent : the value of E in the integral rescaling factor 10^E
!!  BiggsPlanckIntegral   : the value of the Biggs-Planck integral (rescaled by 10^E)
!!  PlanckIntegral        : the value of the Planck integral (rescaled by 10^E)
!!
!!***

subroutine op_BiggsPlanckGroupIntegrate (nSec,                                &
                                         A1Sec,A2Sec,A3Sec,A4Sec,             &
                                         intLimits,                           &
                                                       rescaleBase10Exponent, &
                                                       BiggsPlanckIntegral,   &
                                                       PlanckIntegral         )

  implicit none

  integer, intent (in)  :: nSec
  real,    intent (in)  :: A1Sec     (1:nSec)
  real,    intent (in)  :: A2Sec     (1:nSec)
  real,    intent (in)  :: A3Sec     (1:nSec)
  real,    intent (in)  :: A4Sec     (1:nSec)
  real,    intent (in)  :: intLimits (1:nSec+1)

  real,    intent (out) :: rescaleBase10Exponent
  real,    intent (out) :: BiggsPlanckIntegral
  real,    intent (out) :: PlanckIntegral

  rescaleBase10Exponent = 0.0
  BiggsPlanckIntegral   = 0.0
  PlanckIntegral        = 0.0

  return
end subroutine op_BiggsPlanckGroupIntegrate
