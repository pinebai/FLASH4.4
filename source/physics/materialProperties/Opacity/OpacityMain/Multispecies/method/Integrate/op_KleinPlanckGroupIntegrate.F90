!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_KleinPlanckGroupIntegrate
!!
!! NAME
!!
!!  op_KleinPlanckGroupIntegrate
!!
!! SYNOPSIS
!!
!!  call op_KleinPlanckGroupIntegrate (real    (in)  :: KNintegralFactor,
!!                                     real    (in)  :: B,
!!                                     real    (in)  :: R,
!!                                     real    (in)  :: S,
!!                                     real    (out) :: rescaleBase10Exponent,
!!                                     real    (out) :: KleinPlanckIntegral,
!!                                     real    (out) :: PlanckIntegral)
!!
!! DESCRIPTION
!!
!!  This routine evaluates the following two kinds of integrals using the Romberg integration technique:
!!
!!
!!                                                  S
!!                                                 /
!!                                                |
!!   Klein-Planck Integral =  KNintegralFactor *  |  F[Bx] * x^3 / (e^x - 1) dx
!!                                                |
!!                                               /
!!                                              R
!!
!!      
!!                               S
!!                              /
!!                             |
!!         Planck Integral =   |   x^3 / (e^x - 1) dx
!!                             |
!!                            /
!!                           R
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
!!  As the lower integration limit R gets larger, the exponential function in both integrals
!!  get closer to the point of computational overflow. At a certain R value it is then advantageous
!!  to introduce the rescaling factor e^R during the romberg scheme, such that we have representable
!!  exponential values for a wide range in x. The integrals are then multiplied by e^R and the
!!  resulting factor e^R/(e^x - 1) in the integrand is approximated as e^R/e^x = e^(R-x). Hence
!!  temporary rescaling can only be applied when it is safe to ignore the -1 in (e^x - 1), meaning
!!  that this truncation does not lead to accuracy errors. It is hence dependent on the number of
!!  mantissa digits.
!!
!!  The upper integration limit S might also be way beyond representability of e^S. The code has
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
!!  KleinPlanckIntegral   : the value of the Klein-Planck integral (rescaled by 10^E)
!!  PlanckIntegral        : the value of the Planck integral (rescaled by 10^E)
!!
!!***

subroutine op_KleinPlanckGroupIntegrate (KNintegralFactor,              &
                                         B,                             &
                                         R,S,                           &
                                                 rescaleBase10Exponent, &
                                                 KleinPlanckIntegral,   &
                                                 PlanckIntegral         )

  use Driver_interface, ONLY : Driver_abortFlash

  use op_numericsData,  ONLY : zero,half,one,two,three,four,ten,  &
                               op_eBaseLowestExponent,            &
                               op_eBaseLargestExponent,           &
                               op_eBaseIgnorePlusOneExponent


  use op_integrateData, ONLY : op_minRombergSteps, &
                               op_maxRombergSteps, &
                               op_RombergAccuracy, &
                               op_RombergIntegral, &
                               op_RombergRow

  implicit none

# include "Opacity.h"

  real,    intent (in)  :: KNintegralFactor
  real,    intent (in)  :: B
  real,    intent (in)  :: R,S

  real,    intent (out) :: rescaleBase10Exponent
  real,    intent (out) :: KleinPlanckIntegral
  real,    intent (out) :: PlanckIntegral

  logical :: converged
  logical :: rescaledIntegral
  logical :: totalExactRange
  logical :: totalRationalRange

  integer :: k,m
  integer :: nExactIntEvalSteps
  integer :: nRationalIntEvalSteps
  integer :: nTotalIntEvalSteps
  integer :: step

  real    :: accuracy
  real    :: a,c,x,P,Q,U
  real    :: Bx,Bx2,cx
  real    :: Binv, Binv2, twoBinv2, halfBinv3
  real    :: exphh,expRh
  real    :: FR3,FU3,Fx3
  real    :: h,hh,Rh
  real    :: Integral
  real    :: invExph,invExphh,invExpRh
  real    :: p4m,invp4m
  real    :: previous
  real    :: range
  real    :: scalingExponent
  real    :: sumValue

  real, parameter :: y0 = 1.333333333333333
  real, parameter :: y1 = 4.186213399858848
  real, parameter :: y2 = 5.052212207119021
  real, parameter :: y3 = 1.990894982015258
  real, parameter :: y4 = 1.608689417459069E-1

  real, parameter :: z0 = 1.000000000000000
  real, parameter :: z1 = 5.139659969840541
  real, parameter :: z2 = 8.868481513323712
  real, parameter :: z3 = 5.803866511398383
  real, parameter :: z4 = 1.284363647961355
  real, parameter :: z5 = 5.831913558200177E-2
!
!
!   ...Set the scaling exponent and the upper integration limit U.
!
!
  if (R < op_eBaseIgnorePlusOneExponent) then
      rescaledIntegral = .false.
      scalingExponent  = zero
  else
      rescaledIntegral = .true.
      scalingExponent  = R
  end if

  U = min (scalingExponent - op_eBaseLowestExponent , S )

  rescaleBase10Exponent = scalingExponent / log (ten)
!
!
!   ...Analyze integration range: rational integrand evaluation range (0 =< Bx =< 1)
!                                 exact integrand evaluation range (Bx > 1)
!
!
  range              =  U -  R
  totalExactRange    = (B * R >  one)
  totalRationalRange = (B * U <= one)

  if (.not.totalRationalRange) then
       Binv      = one / B
       Binv2     = Binv * Binv
       twoBinv2  = two  * Binv2
       halfBinv3 = half * Binv * Binv2
  end if
!
!
!   ...Start Romberg integration of Klein-Planck integral.
!
!
  converged = .false.

  x   = R
  Bx  = B * x

  if (Bx <= one) then
      P   = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
      Q   = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
      FR3 = (P/Q) * R * R * R
  else
      Bx2 = Bx + Bx
      a   = one + Bx
      c   = one / (one + Bx2)
      cx  = c * x
      FR3 = (three-(a-two)**2)*halfBinv3*log(c) + cx*(a*a*twoBinv2-cx*x*(a+Bx2))
  end if

  x   = U
  Bx  = B * x

  if (Bx <= one) then
      P   = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
      Q   = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
      FU3 = (P/Q) * U * U * U
  else
      Bx2 = Bx + Bx
      a   = one + Bx
      c   = one / (one + Bx2)
      cx  = c * x
      FU3 = (three-(a-two)**2)*halfBinv3*log(c) + cx*(a*a*twoBinv2-cx*x*(a+Bx2))
  end if

  if (rescaledIntegral) then
      sumValue = FR3 + exp (R-U) * FU3
  else
      sumValue = FR3/(exp(R)-one) + FU3/(exp(U)-one)
  end if

  op_RombergIntegral (0) = half * range * sumValue
  op_RombergRow  (0,OLD) = op_RombergIntegral (0)
!
!
!   ...Romberg iterations.
!
!
  h = range

  do step = 1,op_maxRombergSteps

     h  = h * half
     hh = h + h
     Rh = R + h

     nTotalIntEvalSteps = 2 ** (step-1) - 1

     if (totalRationalRange) then
         nRationalIntEvalSteps = nTotalIntEvalSteps
         nExactIntEvalSteps    = 0
     else if (totalExactRange) then
         nRationalIntEvalSteps = 0
         nExactIntEvalSteps    = nTotalIntEvalSteps
     else
         nRationalIntEvalSteps = max ( 0 , floor (half*(((Binv - R)/h) - one)) )
         nExactIntEvalSteps    = nTotalIntEvalSteps - nRationalIntEvalSteps
     end if
!
!
!   ...Evaluation of sum at all Romberg points.
!
!
     x   = Rh
     Bx  = B * x

     if (Bx <= one) then
         P   = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
         Q   = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
         Fx3 = (P/Q) * x * x * x
     else
         Bx2 = Bx + Bx
         a   = one + Bx
         c   = one / (one + Bx2)
         cx  = c * x
         Fx3 = (three-(a-two)**2)*halfBinv3*log(c) + cx*(a*a*twoBinv2-cx*x*(a+Bx2))
     end if

     if (rescaledIntegral) then

         invExph  = exp (-h)
         invExphh = invExph * invExph
         sumValue = Fx3 * invExph

         do k = 1,nRationalIntEvalSteps
            Rh       = Rh + hh
            x        = Rh
            Bx       = B * x
            P        = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
            Q        = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
            Fx3      = (P/Q) * x * x * x
            invExph  = invExph * invExphh
            sumValue = sumValue + Fx3 * invExph
         end do

         do k = 1,nExactIntEvalSteps
            Rh       = Rh + hh
            x        = Rh
            Bx       = B * x
            Bx2      = Bx + Bx
            a        = one + Bx
            c        = one / (one + Bx2)
            cx       = c * x
            Fx3      = (three-(a-two)**2)*halfBinv3*log(c) + cx*(a*a*twoBinv2-cx*x*(a+Bx2))
            invExph  = invExph * invExphh
            sumValue = sumValue + Fx3 * invExph
         end do

     else

         exphh    = exp (hh)
         expRh    = exp (Rh)
         sumValue = Fx3 / (expRh - one)

         do k = 1,nRationalIntEvalSteps
            Rh       = Rh + hh
            x        = Rh
            Bx       = B * x
            P        = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
            Q        = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
            Fx3      = (P/Q) * x * x * x
            expRh    = expRh * exphh
            sumValue = sumValue + Fx3 / (expRh - one)
         end do

         do k = 1,nExactIntEvalSteps
            Rh       = Rh + hh
            x        = Rh
            Bx       = B * x
            Bx2      = Bx + Bx
            a        = one + Bx
            c        = one / (one + Bx2)
            cx       = c * x
            Fx3      = (three-(a-two)**2)*halfBinv3*log(c) + cx*(a*a*twoBinv2-cx*x*(a+Bx2))
            expRh    = expRh * exphh
            sumValue = sumValue + Fx3 / (expRh - one)
         end do

     end if
!
!
!   ...Construction of new Romberg row.
!
!
     op_RombergRow (0,NEW) = half * op_RombergRow (0,OLD) + h * sumValue

     p4m = one
     do m = 1,step
        p4m    = p4m * four
        invp4m = one / (p4m - one)
        op_RombergRow (m,NEW) = (p4m * op_RombergRow (m-1,NEW) - op_RombergRow (m-1,OLD)) * invp4m
     end do

     op_RombergIntegral (step) = op_RombergRow (step,NEW)
!
!
!   ...check accuracy of the Romberg integral at the current step. Since the Klein-Nishina
!      function F[Bx] is well behaved, we only compare the relative accuracy with the
!      previous step integral. If this relative accuracy is below the specified accuracy
!      limit, the Romberg integral is considered to be converged.
!
!
     if (step > op_minRombergSteps) then

         converged  = .true.
         do m = step,step
            previous  = op_RombergIntegral (m-1)
            Integral  = op_RombergIntegral (m)
            accuracy  = abs ((Integral - previous) / Integral)
            converged = converged  .and. (accuracy  < op_RombergAccuracy)
         end do

     end if

     if (converged) then
!         write (*,*) ' Klein-Planck integral converged after ',step,' Romberg steps'
         exit
     else
         op_RombergRow (0:step,OLD) = op_RombergRow (0:step,NEW)
     end if

  end do
!
!
!   ...Check for convergence after maximum # of Romberg steps. Print a warning, if no convergence
!      was achieved.
!
!
  if (.not.converged) then
       write (*,*) ' Klein-Planck integral convergence problems! '
       write (*,*) ' Klein-Planck integral not converged after ',op_maxRombergSteps,' Romberg steps! '
       write (*,*) ' Highest Klein-Planck integral relative accuracy = ',accuracy
  end if

  KleinPlanckIntegral = KNintegralFactor * Integral
!
!
!   ...Start evaluation of the Planck integral. The info about rescaling, the new upper limit U
!      and the integration range U - R have already been established.
!
!
  converged = .false.
!
!
!   ...Start Romberg integration on Planck integral.
!
!
  if (rescaledIntegral) then
      sumValue = R*R*R + exp (R-U) * U*U*U
  else
      sumValue = R*R*R/(exp(R)-one) + U*U*U/(exp(U)-one)
  end if

  op_RombergIntegral (0) = half * range * sumValue
  op_RombergRow  (0,OLD) = op_RombergIntegral (0)
!
!
!   ...Romberg iterations.
!
!
  h = range

  do step = 1,op_maxRombergSteps

     h  = h * half
     hh = h + h
     Rh = R + h

     nTotalIntEvalSteps = 2 ** (step-1) - 1
!
!
!   ...Evaluation of sum at all Romberg points.
!
!
     if (rescaledIntegral) then

         invExph  = exp (-h)
         invExphh = invExph * invExph
         sumValue = invExph * Rh*Rh*Rh

         do k = 1,nTotalIntEvalSteps
            Rh       = Rh + hh
            invExph  = invExph * invExphh
            sumValue = sumValue + invExph * Rh*Rh*Rh
         end do

     else

         exphh    = exp (hh)
         expRh    = exp (Rh)
         sumValue = Rh*Rh*Rh / (expRh - one)

         do k = 1,nTotalIntEvalSteps
            Rh       = Rh + hh
            expRh    = expRh * exphh
            sumValue = sumValue + Rh*Rh*Rh / (expRh - one)
         end do

     end if
!
!
!   ...Construction of new Romberg row.
!
!
     op_RombergRow (0,NEW) = half * op_RombergRow (0,OLD) + h * sumValue

     p4m = one
     do m = 1,step
        p4m    = p4m * four
        invp4m = one / (p4m - one)
        op_RombergRow (m,NEW) = (p4m * op_RombergRow (m-1,NEW) - op_RombergRow (m-1,OLD)) * invp4m
     end do

     op_RombergIntegral (step) = op_RombergRow (step,NEW)
!
!
!   ...check accuracy of the Romberg integral at the current step. Since the Planck
!      integral is well behaved, we only compare the relative accuracy with the
!      previous step integral. If this relative accuracy is below the specified accuracy
!      limit, the Romberg integral is considered to be converged.
!
!
     if (step > op_minRombergSteps) then

         converged  = .true.
         do m = step,step
            previous  = op_RombergIntegral (m-1)
            Integral  = op_RombergIntegral (m)
            accuracy  = abs ((Integral - previous) / Integral)
            converged = converged  .and. (accuracy  < op_RombergAccuracy)
         end do

     end if

     if (converged) then
!         write (*,*) ' Planck integral converged after ',step,' Romberg steps'
         exit
     else
         op_RombergRow (0:step,OLD) = op_RombergRow (0:step,NEW)
     end if

  end do
!
!
!   ...Check for convergence of the Planck integral. Print a warning, if no convergence
!      was achieved.
!
!
  if (.not.converged) then
       write (*,*)   ' Planck integral not converged after ',op_maxRombergSteps,' Romberg steps! '
       write (*,*)   ' Highest Planck integral relative accuracy = ',accuracy
  end if

  PlanckIntegral = Integral
!
!
!   ...Ready! 
!
!
  return
end subroutine op_KleinPlanckGroupIntegrate
