!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_KleinRosslndGroupIntegrate
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

  real, intent (in)  :: KNintegralFactor
  real, intent (in)  :: B
  real, intent (in)  :: R,S
  real, intent (out) :: rescaleBase10Exponent
  real, intent (out) :: KleinRosslndIntegral
  real, intent (out) :: RosslndIntegral

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
  real    :: a,c,P,Q
  real    :: Binv
  real    :: Bx,Bx2,Bxinv,Bxinvsqr,Bxinvcube
  real    :: expR,expU,exphh,expRh
  real    :: F,FinvR4,FinvU4,Finvx4
  real    :: h,hh
  real    :: Integral
  real    :: invExph,invExphh,invExpRh
  real    :: p4m,invp4m
  real    :: previous
  real    :: R2,R4,Rh,Rh2,Rh4
  real    :: range
  real    :: scalingExponent
  real    :: sumValue
  real    :: U,U2,U4
  real    :: x,xsqr

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
       Binv = one / B
  end if
!
!
!   ...Start Romberg integration of Klein-Rosseland integral.
!
!
  converged = .false.

  x    = R
  xsqr = x * x
  Bx   = B * x

  if (Bx <= one) then
      P      = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
      Q      = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
      FinvR4 = (Q/P) * xsqr * xsqr
  else
      Bx2       = Bx + Bx
      Bxinv     = one / Bx
      Bxinvsqr  = Bxinv * Bxinv
      Bxinvcube = Bxinv * Bxinvsqr
      a         = one + Bx
      c         = one / (one + Bx2)
      F         = (three-(a-two)**2)*half*Bxinvcube*log(c) + c*(two*a*a*Bxinvsqr-c*(a+Bx2))
      FinvR4    = xsqr * xsqr / F
  end if

  x    = U
  xsqr = x * x
  Bx   = B * x

  if (Bx <= one) then
      P      = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
      Q      = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
      FinvU4 = (Q/P) * xsqr * xsqr
  else
      Bx2       = Bx + Bx
      Bxinv     = one / Bx
      Bxinvsqr  = Bxinv * Bxinv
      Bxinvcube = Bxinv * Bxinvsqr
      a         = one + Bx
      c         = one / (one + Bx2)
      F         = (three-(a-two)**2)*half*Bxinvcube*log(c) + c*(two*a*a*Bxinvsqr-c*(a+Bx2))
      FinvU4    = xsqr * xsqr / F
  end if

  if (rescaledIntegral) then
      sumValue = FinvR4 + exp (R-U) * FinvU4
  else
      expR = exp (R)
      expU = exp (U)
      sumValue = FinvR4 * (expR / (expR-one)**2) + FinvU4 * (expU / (expU-one)**2)
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

     h   = h * half
     hh  = h + h
     Rh  = R + h
     Rh2 = Rh * Rh
     Rh4 = Rh2 * Rh2

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
     x    = Rh
     xsqr = x * x
     Bx   = B * x

     if (Bx <= one) then
         P      = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
         Q      = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
         Finvx4 = (Q/P) * xsqr * xsqr
     else
         Bx2       = Bx + Bx
         Bxinv     = one / Bx
         Bxinvsqr  = Bxinv * Bxinv
         Bxinvcube = Bxinv * Bxinvsqr
         a         = one + Bx
         c         = one / (one + Bx2)
         F         = (three-(a-two)**2)*half*Bxinvcube*log(c) + c*(two*a*a*Bxinvsqr-c*(a+Bx2))
         Finvx4    = xsqr * xsqr / F
     end if

     if (rescaledIntegral) then

         invExph  = exp (-h)
         invExphh = invExph * invExph
         sumValue = Finvx4 * invExph

         do k = 1,nRationalIntEvalSteps
            Rh       = Rh + hh
            x        = Rh
            xsqr     = x * x
            Bx       = B * x
            P        = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
            Q        = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
            Finvx4   = (Q/P) * xsqr * xsqr
            invExph  = invExph * invExphh
            sumValue = sumValue + Finvx4 * invExph
         end do

         do k = 1,nExactIntEvalSteps
            Rh        = Rh + hh
            x         = Rh
            Bx        = B * x
            Bx2       = Bx + Bx
            Bxinv     = one / Bx
            Bxinvsqr  = Bxinv * Bxinv
            Bxinvcube = Bxinv * Bxinvsqr
            a         = one + Bx
            c         = one / (one + Bx2)
            F         = (three-(a-two)**2)*half*Bxinvcube*log(c) + c*(two*a*a*Bxinvsqr-c*(a+Bx2))
            Finvx4    = xsqr * xsqr / F
            invExph   = invExph * invExphh
            sumValue  = sumValue + Finvx4 * invExph
         end do

     else

         exphh    = exp (hh)
         expRh    = exp (Rh)
         sumValue = Finvx4 * (expRh / (expRh - one) ** 2)

         do k = 1,nRationalIntEvalSteps
            Rh       = Rh + hh
            x        = Rh
            xsqr     = x * x
            Bx       = B * x
            P        = (((y4*Bx + y3)*Bx + y2)*Bx + y1)*Bx + y0
            Q        = ((((z5*Bx + z4)*Bx + z3)*Bx + z2)*Bx + z1)*Bx + z0
            Finvx4   = (Q/P) * xsqr * xsqr
            expRh    = expRh * exphh
            sumValue = sumValue + Finvx4 * (expRh / (expRh - one) ** 2)
         end do

         do k = 1,nExactIntEvalSteps
            Rh        = Rh + hh
            x         = Rh
            Bx        = B * x
            Bx2       = Bx + Bx
            Bxinv     = one / Bx
            Bxinvsqr  = Bxinv * Bxinv
            Bxinvcube = Bxinv * Bxinvsqr
            a         = one + Bx
            c         = one / (one + Bx2)
            F         = (three-(a-two)**2)*half*Bxinvcube*log(c) + c*(two*a*a*Bxinvsqr-c*(a+Bx2))
            Finvx4    = xsqr * xsqr / F
            expRh     = expRh * exphh
            sumValue  = sumValue + Finvx4 * (expRh / (expRh - one) ** 2)
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
!   ...check accuracy of the Romberg integral at the current step. Since the inverse of the
!      Klein-Nishina function F[Bx] is well behaved, we only compare the relative accuracy with
!      the previous step integral. If this relative accuracy is below the specified accuracy
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
!         write (*,*) ' Klein-Rosseland integral converged after ',step,' Romberg steps'
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
       write (*,*) ' Klein-Rosseland integral convergence problems! '
       write (*,*) ' Klein-Rosseland integral not converged after ',op_maxRombergSteps,' Romberg steps! '
       write (*,*) ' Highest Klein-Rosseland integral relative accuracy = ',accuracy
  end if

  KleinRosslndIntegral = Integral / KNintegralFactor
!
!
!   ...Start evaluation of the Planck integral. The info about rescaling, the new upper limit U
!      and the integration range U - R have already been established.
!
!
  converged = .false.
!
!
!   ...Start Romberg integration on Rosseland integral.
!
!
  R2 = R * R
  R4 = R2 * R2
  U2 = U * U
  U4 = U2 * U2

  if (rescaledIntegral) then
      sumValue = R4 + exp (R-U) * U4
  else
      expR = exp (R)
      expU = exp (U)
      sumValue = R4 * (expR / (expR - one)**2) + U4 * (expU / (expU - one)**2)
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

     h   = h * half
     hh  = h + h
     Rh  = R + h
     Rh2 = Rh * Rh
     Rh4 = Rh2 * Rh2

     nTotalIntEvalSteps = 2 ** (step-1) - 1
!
!
!   ...Evaluation of sum at all Romberg points.
!
!
     if (rescaledIntegral) then

         invExph  = exp (-h)
         invExphh = invExph * invExph
         sumValue = invExph * Rh4

         do k = 1,nTotalIntEvalSteps
            Rh       = Rh + hh
            Rh2      = Rh * Rh
            Rh4      = Rh2 * Rh2
            invExph  = invExph * invExphh
            sumValue = sumValue + invExph * Rh4
         end do

     else

         exphh    = exp (hh)
         expRh    = exp (Rh)
         sumValue = Rh4 * (expRh / (expRh - one)**2)

         do k = 1,nTotalIntEvalSteps
            Rh       = Rh + hh
            Rh2      = Rh * Rh
            Rh4      = Rh2 * Rh2
            expRh    = expRh * exphh
            sumValue = sumValue + Rh4 * (expRh / (expRh - one)**2)
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
!   ...check accuracy of the Romberg integral at the current step. Since the Roseeland
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
!         write (*,*) ' Rosseland integral converged after ',step,' Romberg steps'
         exit
     else
         op_RombergRow (0:step,OLD) = op_RombergRow (0:step,NEW)
     end if

  end do
!
!
!   ...Check for convergence of the Rosseland integral. Print a warning, if no convergence
!      was achieved.
!
!
  if (.not.converged) then
       write (*,*)   ' Rosseland integral not converged after ',op_maxRombergSteps,' Romberg steps! '
       write (*,*)   ' Highest Rosseland integral relative accuracy = ',accuracy
  end if

  RosslndIntegral = Integral
!
!
!   ...Ready! 
!
!
  return
end subroutine op_KleinRosslndGroupIntegrate
