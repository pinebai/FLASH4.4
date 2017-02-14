!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_BiggsPlanckGroupIntegrate
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

  use Driver_interface, ONLY : Driver_abortFlash

  use op_numericsData,  ONLY : zero,half,one,four,ten,        &
                               op_eBaseLowestExponent,        &
                               op_eBaseLargestExponent,       &
                               op_eBaseIgnorePlusOneExponent

  use op_integrateData, ONLY : op_minRombergSteps, &
                               op_maxRombergSteps, &
                               op_RombergAccuracy, &
                               op_RombergIntegral, &
                               op_RombergRow

  implicit none

# include "Opacity.h"

  integer, intent (in)  :: nSec
  real,    intent (in)  :: A1Sec     (1:nSec)
  real,    intent (in)  :: A2Sec     (1:nSec)
  real,    intent (in)  :: A3Sec     (1:nSec)
  real,    intent (in)  :: A4Sec     (1:nSec)
  real,    intent (in)  :: intLimits (1:nSec+1)

  real,    intent (out) :: rescaleBase10Exponent
  real,    intent (out) :: BiggsPlanckIntegral
  real,    intent (out) :: PlanckIntegral

  logical :: converged
  logical :: rescaledIntegral

  integer :: k,m,n
  integer :: nIntEvalSteps
  integer :: step

  real    :: A1,A2,A3,A4
  real    :: accuracy
  real    :: basicExponent
  real    :: exphh,expRh
  real    :: h,hh,Rh
  real    :: Integral
  real    :: invExph,invExphh,invExpRh
  real    :: p4m,invp4m
  real    :: previous
  real    :: R,S,U
  real    :: range
  real    :: scalingExponent
  real    :: sumValue
!
!
!   ...Set the 1st (basic) rescaling exponent.
!
!
  R = intLimits (1)

  if (R < op_eBaseIgnorePlusOneExponent) then
      basicExponent = zero
  else
      basicExponent = R
  end if

  rescaleBase10Exponent = basicExponent / log (ten)
!
!
!   ...Start evaluation of the Biggs-Planck integral. Outer loop over all sections.
!      Skip section, if all expansion coefficients are equal to zero.
!      Skip rest of sections, if rescaling factor is too small.
!
!
  BiggsPlanckIntegral = zero

  do n = 1,nSec

     R = intLimits (n)
     S = intLimits (n+1)

     A1 = A1Sec (n)
     A2 = A2Sec (n)
     A3 = A3Sec (n)
     A4 = A4Sec (n)

     if (A1==zero .and. A2==zero .and. A3==zero .and. A4==zero) then
         cycle
     end if

     if (basicExponent - R < op_eBaseLowestExponent) then
         exit
     end if

     converged = .false.
!
!
!   ...Start evaluation of current n-th section. Set the scaling exponent and
!      the upper integration limit U.
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
!
!
!   ...Start Romberg integration on n-th section of Biggs-Planck integral.
!
!
     range = U - R

     if (rescaledIntegral) then
         sumValue = (R*(A1*R+A2)+A3+A4/R) + exp (R-U) * (U*(A1*U+A2)+A3+A4/U)
     else
         sumValue = (R*(A1*R+A2)+A3+A4/R)/(exp(R)-one) + (U*(A1*U+A2)+A3+A4/U)/(exp(U)-one)
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

        nIntEvalSteps = 2 ** (step-1) - 1
!
!
!   ...Evaluation of sum at all Romberg points.
!
!
        if (rescaledIntegral) then

            invExph  = exp (-h)
            invExphh = invExph * invExph
            sumValue = invExph * (Rh*(A1*Rh+A2)+A3+A4/Rh)

            do k = 1,nIntEvalSteps
               Rh       = Rh + hh
               invExph  = invExph * invExphh
               sumValue = sumValue + invExph * (Rh*(A1*Rh+A2)+A3+A4/Rh)
            end do

        else

            exphh    = exp (hh)
            expRh    = exp (Rh)
            sumValue = (Rh*(A1*Rh+A2)+A3+A4/Rh) / (expRh - one)

            do k = 1,nIntEvalSteps
               Rh       = Rh + hh
               expRh    = expRh * exphh
               sumValue = sumValue + (Rh*(A1*Rh+A2)+A3+A4/Rh) / (expRh - one)
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
!   ...check accuracy of the Romberg integral at the current step. Since the Biggs-Planck
!      integral has a rather complicated integrand, to assure convergence the set of the latest
!      four integrals are analyzed in terms of relative accuracy obtained at each of the
!      last four steps. If all three relative accuracies are below the specified accuracy
!      limit, the Romberg integral for the current section is considered to be converged.
!
!
         if (step > op_minRombergSteps) then

             converged  = .true.
             do m = step-2,step
                previous  = op_RombergIntegral (m-1)
                Integral  = op_RombergIntegral (m)
                accuracy  = abs ((Integral - previous) / Integral)
                converged = converged  .and. (accuracy  < op_RombergAccuracy)
             end do

         end if

         if (converged) then
!             write (*,*) ' Sec ',n,' converged after ',step,' Romberg steps'
             exit
         else
             op_RombergRow (0:step,OLD) = op_RombergRow (0:step,NEW)
         end if

     end do
!
!
!   ...Check for convergence of n-th section. Print a warning, if no convergence
!      was achieved.
!
!
!
!   ...Print section (activate for testing purposes).
!
!
!
     if (.not.converged) then
          write (*,*) n,'-th section of Biggs-Planck integral convergence problems! '
          write (*,*)   ' Biggs-Planck integral not converged after ',op_maxRombergSteps,' Romberg steps! '
          write (*,*)   ' Highest Biggs-Planck integral relative accuracy = ',accuracy
     end if
!
!
!   ...Add current properly rescaled Biggs-Planck integral section to overall integral.
!
!
     BiggsPlanckIntegral = BiggsPlanckIntegral + exp (basicExponent - scalingExponent) * Integral
    
  end do
!
!
!   ...Start evaluation of the Planck integral.
!
!
  R = intLimits (1)
  S = intLimits (nSec+1)

  converged = .false.

  if (R < op_eBaseIgnorePlusOneExponent) then
      rescaledIntegral = .false.
      U = min (- op_eBaseLowestExponent , S )
  else
      rescaledIntegral = .true.
      U = min (R - op_eBaseLowestExponent , S )
  end if
!
!
!   ...Start Romberg integration on Planck integral.
!
!
  range = U - R

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

     nIntEvalSteps = 2 ** (step-1) - 1
!
!
!   ...Evaluation of sum at all Romberg points.
!
!
     if (rescaledIntegral) then

         invExph  = exp (-h)
         invExphh = invExph * invExph
         sumValue = invExph * Rh*Rh*Rh

         do k = 1,nIntEvalSteps
            Rh       = Rh + hh
            invExph  = invExph * invExphh
            sumValue = sumValue + invExph * Rh*Rh*Rh
         end do

     else

         exphh    = exp (hh)
         expRh    = exp (Rh)
         sumValue = Rh*Rh*Rh / (expRh - one)

         do k = 1,nIntEvalSteps
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
end subroutine op_BiggsPlanckGroupIntegrate
