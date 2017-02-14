!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_ecapnuc
!!
!! NAME
!!
!!   bn_ecapnuc
!!
!! SYNOPSIS
!!   bn_ecapnuc(real(IN)  ::etakep,
!!              real(IN)  ::temp,
!!              real(OUT) ::rpen,
!!              real(OUT) ::rnep,
!!              real(OUT) ::spen,
!!              real(OUT) ::snep)
!!
!! DESCRIPTION
!!
!!  routine bn_ecapnuc computes neutron and proton electron capture rates  
!!  
!!  given the electron degeneracy parameter etakep (chemical potential
!!  without the electron's rest mass divided by kt) and the temperature temp,
!!  this routine calculates rates for 
!!  electron capture on protons rpen (captures/sec/proton),
!!  positron capture on neutrons rnep (captures/sec/neutron), 
!!  and their associated neutrino energy loss rates 
!!  spen (ergs/sec/proton) and snep (ergs/sec/neutron)
!!  
!! ARGUMENTS
!!  etakep -- electron degeneracy parameter
!!  temp   -- temperature
!!  rpen   -- electron capture on protons rate
!!  rnep   -- positron capture on neutrons rate
!!  spen   -- neutrino energy loss rate per proton?
!!  snep   -- neutrino energy loss rate per neutron?
!!
!!***

subroutine bn_ecapnuc(etakep,temp,rpen,rnep,spen,snep)

  implicit none
 
  !!  declare  arguments
  real, intent(IN)     :: etakep, temp
  real, intent(OUT)    :: rpen,rnep,spen,snep

  !!  declare local variables
  integer, save ::           iflag
  real,save ::  t9,t5,qn,etaef,etael,zetan,eta,etael2, & 
       &                 etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g, & 
       &                 f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4, & 
       &                 fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1, & 
       &                 facv2,facv3,facv4,rjv1,rjv2,spenc,snepc, & 
       &                 exeta,zetan2,f0,etael5 
  real, parameter      :: qn1    = -2.0716446e-06, & 
       &                  ft     = 1083.9269e0, & 
       &                  twoln  = 0.6931472e0, & 
       &                  cmk5   = 1.3635675d-49, & 
       &                  cmk6   = 2.2993864d-59, & 
       &                  bk     = 1.38062e-16, & 
       &                  pi     = 3.1415927e0, & 
       &                  pi2    = pi * pi, & 
       &                  qn2    = 2.0716446e-06, & 
       &                  c2me   = 8.1872665e-07, & 
       &                  xmp    = 1.6726485e-24, & 
       &                  xmn    = 1.6749543e-24, & 
       &                  qndeca = 1.2533036e-06, & 
       &                  tmean  = 935.14e0

  !!  
  !!  tmean and qndeca are the mean lifetime and decay energy of the nuetro
  !!  xmp,xnp are masses of the p and n in grams.
  !!  c2me is the constant used to convert the neutrino energy
  !!  loss rate from mec2/s (as in the paper) to ergs/particle/sec.
  !!  
  !!  initialize
  rpen  = 0.0e0
  rnep  = 0.0e0
  spen  = 0.0e0
  snep  = 0.0e0
  t9    = temp * 1.0e-9
  iflag = 0
  qn    = qn1

  !!  chemical potential including the electron rest mass
  etaef = etakep + c2me/bk/temp

  !!  iflag=1 is for electrons,  iflag=2 is for positrons

502 iflag = iflag + 1
  if (iflag.eq.1) etael = qn2/bk/temp
  if (iflag.eq.2) etael = c2me/bk/temp
  if (iflag.eq.2) etaef = -etaef

  t5    = temp*temp*temp*temp*temp
  zetan = qn/bk/temp
  eta   = etaef - etael

  !!  protect from overflowing with large eta values
  if (eta .le. 6.8e+02) then
     exeta = exp(eta)
  else 
     exeta = 0.0e0
  end if
  etael2 = etael*etael
  etael3 = etael2*etael
  etael4 = etael3*etael
  etael5 = etael4*etael
  zetan2 = zetan*zetan
  if (eta .le. 6.8e+02) then
     f0 = log(1.0e0 + exeta)
  else
     f0 = eta
  end if

  !!  if eta le. 0., the following fermi integrals apply
  f1l = exeta
  f2l = 2.0e0   * f1l
  f3l = 6.0e0   * f1l
  f4l = 24.0e0  * f1l
  f5l = 120.0e0 * f1l

  !!  if eta gt. 0., the following fermi integrals apply:
  f1g = 0.0e0
  f2g = 0.0e0
  f3g = 0.0e0
  f4g = 0.0e0
  f5g = 0.0e0
  if (eta .gt. 0.0) then
     exmeta = exp(-eta)
     eta2   = eta*eta
     eta3   = eta2*eta
     eta4   = eta3*eta
     f1g = 0.5e0*eta2 + 2.0e0 - exmeta
     f2g = eta3/3.0e0 + 4.0e0*eta + 2.0e0*exmeta
     f3g = 0.25e0*eta4 + 0.5e0*pi2*eta2 + 12.0e0 - 6.0e0*exmeta
     f4g = 0.2e0*eta4*eta + 2.0e0*pi2/3.0e0*eta3 + 48.0e0*eta & 
          &       + 24.0e0*exmeta
     f5g = eta4*eta2/6.0e0 + 5.0e0/6.0e0*pi2*eta4  & 
          &       + 7.0e0/6.0e0*pi2*eta2  + 240.0e0 -120.e0*exmeta
  end if

  !!  factors which are multiplied by the fermi integrals
  fac3 = 2.0e0*zetan + 4.0e0*etael
  fac2 = 6.0e0*etael2 + 6.0e0*etael*zetan + zetan2
  fac1 = 4.0e0*etael3 + 6.0e0*etael2*zetan + 2.0e0*etael*zetan2
  fac0 = etael4 + 2.0e0*zetan*etael3 + etael2*zetan2

  !!  electron capture rates onto protons with no blocking
  rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
  rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

  !!  neutrino emission rate for electron capture:
  facv4 = 5.0e0*etael + 3.0e0*zetan
  facv3 = 10.0e0*etael2 + 12.0e0*etael*zetan + 3.0e0*zetan2
  facv2 = 10.0e0*etael3 + 18.0e0*etael2*zetan & 
       &        + 9.0e0*etael*zetan2 + zetan2*zetan
  facv1 = 5.0e0*etael4 + 12.0e0*etael3*zetan  & 
       &        + 9.0e0*etael2*zetan2 + 2.0e0*etael*zetan2*zetan
  facv0 = etael5 + 3.0e0*etael4*zetan & 
       &        + 3.0e0*etael3*zetan2 + etael2*zetan2*zetan
  rjv1  = f5l + facv4*f4l + facv3*f3l & 
       &        + facv2*f2l + facv1*f1l + facv0*f0
  rjv2  = f5g + facv4*f4g + facv3*f3g & 
       &        + facv2*f2g + facv1*f1g + facv0*f0

  !!  for electrons capture onto protons
  if (iflag.eq.2) go to 503
  if (eta.gt.0.) go to 505
  rpen  = twoln*cmk5*t5*rie1/ft
  spen  = twoln*cmk6*t5*temp*rjv1/ft
  spenc = twoln*cmk6*t5*temp*rjv1/ft*c2me
  go to 504

505 rpen = twoln*cmk5*t5*rie2/ft
  spen = twoln*cmk6*t5*temp*rjv2/ft
  spenc = twoln*cmk6*t5*temp*rjv2/ft*c2me

504 continue
  qn = qn2
  go to 502

  !!  for positrons capture onto neutrons
503 if (eta.gt.0.) go to 507
  rnep  = twoln*cmk5*t5*rie1/ft
  snep  = twoln*cmk6*t5*temp*rjv1/ft
  snepc = twoln*cmk6*t5*temp*rjv1/ft*c2me
  !      if (rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
  go to 506

507 rnep  = twoln*cmk5*t5*rie2/ft
  snep  = twoln*cmk6*t5*temp*rjv2/ft
  snepc = twoln*cmk6*t5*temp*rjv2/ft*c2me
  !      if (rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean

506 continue

  return
end subroutine bn_ecapnuc
