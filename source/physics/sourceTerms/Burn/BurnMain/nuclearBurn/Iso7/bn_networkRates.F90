!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_networkRates
!!
!! NAME
!!  
!!  bn_networkRates
!!
!!
!! SYNOPSIS
!! 
!!  bn_networkRates()
!!
!!  
!! DESCRIPTION
!!
!!  this routine generates nuclear reaction rates for the iso7 network.
!!       r1616 = o16+o16      r1212 = c12+c12      r3a   = 3a-c12
!!       rg3a  = c12-3a       roga  = o16(ga)      roag  = o16(ag)   
!!       rnega = ne20(ga)     rneag = ne20(ag)     rmgga = mg24(ga)     
!!       rmgag = mg24(ag)     rsiga = si28(ga)     r1216 = o16+c12
!!       rcaag = ca40(ag)     rtiga = ti44(ga)     rcag  = c12(ag)
!!  
!!
!!***

subroutine bn_networkRates

  use Burn_dataEOS, ONLY:  btemp, den
  use Burn_data
  use bn_dataIso7

#include "Flash.h"

   implicit none

 !
  !..declare  
  integer          i
  real tt9,t9r,t9,t912,t913,t923,t943,t953,t932,                    &
       &                 t92,t93,t972,t9r32,t9i,t9i13,t9i23,t9i32,        &
       &                 t9i12,t9ri,term,term1,term2,term3,rev,           &
       &                 r2abe,t9a,t9a13,t9a56,t9a23,gt9h,gca40,gti44,    &
       &                 rbeac,oneth,fivsix,sav(20)
  parameter        (oneth = 1.0e0/3.0e0, fivsix=5.0e0/6.0e0)

  !..zero the rates
  do i=1,nrat
     ratraw(i) = 0.0e0
  enddo

  !..some temperature factors 
  !..limit t9 to 10, except for the reverse ratios
  tt9   = btemp * 1.0e-9
  t9r   = tt9
  t9    = min(tt9,10.0e0)
  t912  = sqrt(t9)
  t913  = t9**oneth
  t923  = t913*t913
  t943  = t9*t913
  t953  = t9*t923
  t932  = t9*t912
  t92   = t9*t9
  t93   = t9*t92
  t972  = t912*t93
  t9r32 = t9r*sqrt(t9r)
  t9i   = 1.0e0/t9
  t9i13 = 1.0e0/t913
  t9i23 = 1.0e0/t923
  t9i32 = 1.0e0/t932
  t9i12 = 1.0e0/t912
  t9ri  = 1.0e0/t9r


  !..12c(ag)16o; note 1.7 times cf88 value
  term          = 1.04e8/(t92*(1.0e0 + 0.0489e0*t9i23)**2) *        &
       &                  exp(-32.120e0*t9i13 - t92/12.222016e0)          &
       &                + 1.76e8/(t92*(1.0e0 + 0.2654e0*t9i23)**2) *      &
       &                  exp(-32.120e0*t9i13)                            &
       &                + 1.25e3*t9i32*exp(-27.499*t9i)                   &
       &                + 1.43e-2*t92*t93*exp(-15.541*t9i)
  term          = 1.7e0 * term
  ratraw(ircag) = term * den 
  rev           = 5.13e10 * t9r32 * exp(-83.111*t9ri)
  ratraw(iroga) = rev * term


  !..triple alpha to c12 (has been divided by six)
  !..rate revised according to caughlan and fowler (nomoto ref.) 1988 q = -0.092
  r2abe = (7.40e+05 * t9i32)* exp(-1.0663*t9i) +                    &
       &        4.164e+09 * t9i23 * exp(-13.49*t9i13 - t92/0.009604) *    &
       &        (1.0e0 + 0.031*t913 + 8.009*t923 + 1.732*t9 +             &
       &        49.883*t943 + 27.426*t953)
  !..q = 7.367
  rbeac = (130.*t9i32) * exp(-3.3364*t9i) +                         &
       &        2.510e+07 * t9i23 * exp(-23.57*t9i13 - t92/0.055225) *    &
       &        (1.0e0 + 0.018*t913 + 5.249*t923 + 0.650*t9 +             &
       &        19.176*t943 + 6.034*t953)
  !..q = 7.275
  if (t9.gt.0.08) then
     term = 2.90e-16 * (r2abe*rbeac) +                                &
          &        0.1 * 1.35e-07 * t9i32 * exp(-24.811*t9i)
  else
     term = 2.90e-16*(r2abe*rbeac) *                                  &
          &        (0.01 + 0.2*(1.0e0 + 4.0e0*exp(-(0.025*t9i)**3.263)) /    &
          &        (1.0e0 + 4.0e0*exp(-(t9/0.025)**9.227))) +                &
          &        0.1 * 1.35e-07 * t9i32 * exp(-24.811*t9i)
  end if
  ratraw(ir3a)  = term * (den*den)/6.0e0
  rev           = 2.00e+20*exp(-84.424*t9ri)
  ratraw(irg3a) = rev*(t9r*t9r*t9r) * term


  !..c12 + c12 reaction; see cf88 references 47-49
  t9a            = t9/(1.0e0 + 0.0396*t9)
  t9a13          = t9a**oneth
  t9a56          = t9a**fivsix
  term           = 4.27e+26 * t9a56/t932 *                          &
       &                 exp(-84.165/t9a13-2.12e-03*t9*t9*t9)
  ratraw(ir1212) = 0.5e0 * den * term


  !..c12 + o16 reaction;  q = 16.755; valid for t9 .gt. 0.5
  !..y(nc12)*y(no16)*rc28 is the rate of formation of the si28 compound nucleus
  if (t9.ge.0.5) then
     t9a   = t9/(1.+0.055*t9)
     t9a13 = t9a**oneth
     t9a23 = t9a13*t9a13
     t9a56 = t9a**fivsix
     term  = 1.72e+31 * t9a56 * t9i32 * exp(-106.594/t9a13) /         &
          &         (exp(-0.18*t9a*t9a) + 1.06e-03*exp(2.562*t9a23))
     ratraw(ir1216) = den * term
  else
     ratraw(ir1216) = 0.0e0
  endif


  !..16o+16o rate; q = 16.542; 
  !..y16*y16*r32 is rate of formation of 32s compound nucleus
  term           = 7.10e36 * t9i23 *                                &
       &                 exp(-135.93 * t9i13 - 0.629*t923 -               &
       &                 0.445*t943 + 0.0103*t9*t9)
  ratraw(ir1616) = 0.5e0 * den * term                     


  !..16o(ag)20ne + inverse
  term1          = 9.37e9 * t9i23 *                                 &
       &                 exp(-39.757*t9i13- t92/2.515396)
  term2          = 62.1 * t9i32 * exp(-10.297*t9i) +                &
       &                 538.0e0 * t9i32 * exp(-12.226*t9i) +             &
       &                 13.0e0 * t92 * exp(-20.093*t9i)
  term           = term1 + term2
  ratraw(iroag)  = den * term
  rev            = 5.65e+10*t9r32*exp(-54.937*t9ri)
  ratraw(irnega) = rev * term


  !..20ne(ag)24mg + inverse
  term1          = 4.11e+11 * t9i23 *                               &
       &                 exp(-46.766*t9i13- t92/4.923961) *               &
       &                 (1.0e0 + 0.009*t913 + 0.882*t923 + 0.055*t9 +    &
       &                 0.749*t943 + 0.119*t953)
  term2          = 5.27e+03 * t9i32 * exp(-15.869*t9i) +            &
       &                 6.51e+03 * t912 * exp(-16.223*t9i)
  term3          = 0.1e0 * (42.1 * t9i32 * exp(-9.115*t9i) +        &
       &                 32.0 * t9i23 * exp(-9.383*t9i))
  term           = (term1+term2+term3)/                             &
       &                 (1.0e0 + 5.0e0*exp(-18.960*t9i))
  ratraw(irneag) = den * term
  rev            = 6.01e+10 * t9r32 * exp(-108.059*t9ri)
  ratraw(irmgga) = term * rev


  !..24mg(ag)28si + inverse
  term1          = (1.0e0 + 5.0e0*exp(-15.882*t9i))
  term           = (4.78e+01 * t9i32 * exp(-13.506*t9i) +           &
       &                 2.38e+03 * t9i32 * exp(-15.218*t9i) +            &
       &                 2.47e+02 * t932 * exp(-15.147*t9i) +             &
       &                 0.1 * (1.72e-09 * t9i32 * exp(-5.028*t9i) +      &
       &                 1.25e-03 * t9i32 * exp(-7.929*t9i) +             &
       &                 2.43e+01 * t9i * exp(-11.523*t9i)))/term1
  ratraw(irmgag) = den * term
  rev            = 6.27e+10 * t9r32 * exp(-115.862*t9ri)
  ratraw(irsiga) = rev * term 


  !..40ca(ag)44ti + inverse 
  !..take into account the temp dependence of the partition functions
  gca40=1.0e0 +exp((-4.150E+01 + 1.636E+00*t9 + 1.483E-01*t92)*t9ri)
  gti44=1.0e0 +exp((-1.111E+01 + 6.293E-01*t9 + 1.732E-01*t92)*t9ri)
  term           = 4.66e+24 * t9i23 * exp(-76.435 * t9i13 *         &
       &                 (1.0e0 + 1.650e-02*t9 + 5.973e-03*t92 -          &
       &                  3.889e-04*t93))
  rev            = 6.843e+10 * t9r32 * exp(-59.510*t9ri)
  ratraw(ircaag) = den*term
  ratraw(irtiga) = rev*term * gca40/gti44

  return
end subroutine  bn_networkRates


