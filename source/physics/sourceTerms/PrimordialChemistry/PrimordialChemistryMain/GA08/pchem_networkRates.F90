!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_networkRates
!!
!! NAME
!!
!! pchem_networkRates
!!
!! SYNOPSIS
!!
!! pchem_networkRates()
!!
!!
!! DESCRIPTION
!! 
!!   routine networkRates generates the raw reaction rates for 
!!  
!!
!!
!!   see the pchem_initNetwork to see all the reactions used.
!!   
!!
!!***

subroutine pchem_networkRates
   
   use PrimordialChemistry_data
   use PrimordialChemistry_dataEOS
   use pchem_data
                                

   implicit none

#include "constants.h"
#include "Flash.h"

  integer       i
  real  lt1, lt2, lt3, lt4, lt5, lt6, lt7, lt8, lt9
  real   lnt1, lnt2, lnt3, lnt4, lnt5, lnt6, lnt7, lnt8, lnt9 
  real   te, mp, t3
  real  term, T, naterm, narate

!! Above is lt#, which is log10(Temp)^#, lnt# are ln(T)^#

!! zero the rates
  do i=1,nrat
     ratraw(i) = 0.0e0
  enddo
 
!! Say what lt#'s are

  lt1 = log10(ctemp)
  lt2 = log10(ctemp)**2
  lt3 = log10(ctemp)**3
  lt4 = log10(ctemp)**4
  lt5 = log10(ctemp)**5
  lt6 = log10(ctemp)**6
  lt7 = log10(ctemp)**7
  lt8 = log10(ctemp)**8

  
  te = ctemp*0.00008617  !This is kb*T, in units of electron-volts
  T = ctemp
  lnt1 = log(te)
  lnt2 = log(te)**2
  lnt3 = log(te)**3
  lnt4 = log(te)**4
  lnt5 = log(te)**5
  lnt6 = log(te)**6 
  lnt7 = log(te)**7
  lnt8 = log(te)**8
  lnt9 = log(te)**9
  t3 = ctemp / 300.0
  
 
  mp = 1.67262e-24 !!Mass of proton (used as nucleon) gms
  mp = 1/(6.02214e23) !!Changing this, want each 'rate' to be multiplied by Na, avogrado's number. easy way to do this

  narate = den / mp
  naterm = den / mp


!!R001: H + e --> HM + gam
if(ctemp .gt. 6000.0) then
   term = 10**(-16.4199 + 0.1998*lt2 - 0.005447*lt4 + 0.000040415*lt6)
else
    term = 10**(-17.8450 + 0.762*lt1 + 0.1523*lt2 - 0.03274*lt3)   
endif
   ratraw(iR001) = term * narate

!!R002: H + HM --> H2 + e
   term = 1.3e-9   !!*t3**(-0.1)
   ratraw(iR002) = term * narate

!!R003: H + HP --> H2P + gam
   term = 10**(-19.38 - 1.523*lt1 + 1.118*lt2 - 0.1269*lt3)
   ratraw(iR003) = term * narate

!!R004: H + H2P --> H2 + HP
   term = 6.4e-10
   ratraw(iR004) = term * narate

!!R005: HP + HM --> H + H
   term = 2.4e-6*ctemp**(-0.5)*(1.0+5.0e-5*ctemp)
   ratraw(iR005) = term * narate

!!R006 : H2P + e --> H + H
if(ctemp .gt. 617.0) then
   term = 1.32e-6*ctemp**(-0.76)
else
   term = 1.0e-8
endif
   ratraw(iR006) = term * narate

!!R007: H2 + HP --> H + H2P
if(ctemp .lt. 110.0) then
    term = 0.0e0
else if(ctemp .gt. 3.0e4) then
    term = (-3.3232183e-7 + 3.3735382e-7*log(3.0e4) - 1.4491368e-7*log(3.0e4)**2 &
          + 3.4172805e-8*log(3.0e4)**3- 4.7813720e-9*log(3.0e4)**4 + 3.9731542e-10*log(3.0e4)**5 &
          -1.8171411e-11*log(3.0e4)**6 + 3.5311932e-13*log(3.0e4)**7)*exp(-21237.15/3.0e4)
else
    term = (-3.3232183e-7 + 3.3735382e-7*log(ctemp) - 1.4491368e-7*log(ctemp)**2 &
          + 3.4172805e-8*log(ctemp)**3- 4.7813720e-9*log(ctemp)**4 + 3.9731542e-10*log(ctemp)**5 &
          -1.8171411e-11*log(ctemp)**6 + 3.5311932e-13*log(ctemp)**7)*exp(-21237.15/ctemp)
endif
   ratraw(iR007) = term * narate

!!R008: H2 + e --> H + H + e
!    term = nv*log10(4.49e-9*ctemp**(0.11)*exp(-101858.0/ctemp))
!    term = term + nlte*log10(1.91e-9*ctemp**(0.136)*exp(-53407.1/ctemp))
!    term = 10**(term)
if(pchem_rcCase .eq. 1) then 
    if(ctemp .lt. 1100.0) then
       term = 0.0e0
    else
       term = 4.49e-9*ctemp**(0.11)*exp(-101858.0/ctemp)
    endif
else
    if(ctemp .lt. 4050.0) then
       term = 0.0e0
    else
       term = 1.91e-9*ctemp**(0.136)*exp(-53407.1/ctemp) 
    endif
endif

    ratraw(iR008) = term * narate

!!R009: H2 + H --> H + H + H
!   term = nv*log10(6.67e-12*ctemp**(0.5)*exp(-(1.0+63593.0/ctemp)))
!   term = term + nlte*log10(3.52e-9*exp(-43900.0/ctemp))  !!Simon's Rate
!   term = 10**(term)
if(pchem_rcCase .eq. 1) then
   if(ctemp .lt. 700.0) then
      term = 0.0e0
   else
      term = 6.67e-12*ctemp**(0.5)*exp(-(1.0+63593.0/ctemp))
   endif
else
   if(ctemp .lt. 500.0) then
      term = 0.0e0
   else
      term = 3.52e-9*exp(-43900.0/ctemp)
   endif
endif
    ratraw(iR009) = term * narate

!!R010: H2 + H2 --> H + H + H2
!   term = (5.996e-30*ctemp**(4.1881)/(1.0+6.761e-6*ctemp)**(5.6881))*exp(-54657.4/ctemp)
!   term = nv*log10(term)
!   term = term + nlte*log10(1.3e-9*exp(-53300.0/ctemp)) !!Simon's Rate
!   term = 10**(term)
if(pchem_rcCase .eq. 1) then
   if(ctemp .lt. 750.0) then
      term = 0.0e0
   else
      term = (5.996e-30*ctemp**(4.1881)*exp(-54657.4/ctemp))/(1.0+6.761e-6*ctemp)
!!      term = 5.996e-30*ctemp**(4.1881)/(1.0+6.761e-6*ctemp)**(5.6881)*exp(-54657.4/ctemp)
   endif
else
   if(ctemp .lt. 600.0) then
      term = 0.0e0
   else
      term = 1.3e-9*exp(-53300.0/ctemp)
   endif
endif
   ratraw(iR010) = term * narate

!!R011: H2 + He --> H + H + He
!   term = nv*(-27.029+3.801*lt1-29487.0/ctemp)
!   term = term + nlte*(-2.729-1.75*lt1-23474.0/ctemp)   !!!6.6e-10*ctemp**(0.115)*exp(-52000.0/ctemp))
!   term = 10**(term)
if(pchem_rcCase .eq. 1) then
   if(ctemp .lt. 900.0) then
      term = 0.0e0
   else
      term = 10**(-27.029+3.801*lt1-29487.0/ctemp)
   endif
else
   if(ctemp .lt. 600.0) then
      term = 0.0e0
   else
      term = 10**(-2.729-1.75*lt1-23474.0/ctemp)
   endif
endif
   ratraw(iR011) = term * narate

!! R012: H + e --> HP + e + e
if(ctemp .lt. 2800.0) then
   term = 0.0e0
else
   term = exp(-3.271396786e1 + 1.35365560e1*lnt1 - 5.73932875e0*lnt2 &
           +1.56315498e0*lnt3 - 2.87705600e-1*lnt4 + 3.48255977e-2*lnt5 &
           -2.63197617e-3*lnt6 + 1.11954395e-4*lnt7 -2.03914985e-6*lnt8)
endif
   ratraw(iR012) = term * naterm

!!R013: HP + e --> H + gamma
if(pchem_ccCase .eq. 0) then
     term = 1.269e-13*(315614.0/ctemp)**(1.503)*(1.0+(604625.0/ctemp)**0.470)**(-1.923) !! Trying case A
else if (pchem_ccCase .eq. 1) then
     term = 2.753e-14*(315614.0/ctemp)**(1.500)*(1.0+(115188.0/ctemp)**0.407)**(-2.242)  !! This is case B
else
     term = 2.753e-14*(315614.0/ctemp)**(1.500)*(1.0+(115188.0/ctemp)**0.407)**(-2.242)  !! This is case B
endif
  ratraw(iR013) = term * narate

!!R014: HM + e --> H + e + e
if(ctemp .lt. 100.0) then
    term = 0.0e0
else
    term =exp(-1.801849334e1 + 2.36085220*lnt1 &
        -2.82744300e-1*lnt2 + 1.62331664e-2*lnt3 &
        -3.36501203e-2*lnt4 + 1.17832978e-2*lnt5 - 1.65619470e-3*lnt6 &
        +1.06827520e-4*lnt7 - 2.63128581e-6*lnt8)
endif
    ratraw(iR014) = term * narate

!!R015: HM + H --> H + H + e
if(te .gt. 0.1) then
    term = exp(-2.0372609e1 + 1.13944933e0*lnt1 - 1.4210135e-1*lnt2 &
           +8.4644554e-3*lnt3 - 1.4327641e-3*lnt4 + 2.0122503e-4*lnt5 &
           +8.6639632e-5*lnt6 - 2.5850097e-5*lnt7 + 2.4555012e-6*lnt8 &
           -8.0683825e-8*lnt9)
else
    term = 2.5634e-9*te**(1.78186)
endif
  ratraw(iR015) = term * narate

!!R016: HP + HM --> H2P + e
if(ctemp .gt. 8000.0) then
    term = 9.6e-7*(ctemp)**(-0.90)
else
    term = 6.9e-9*(ctemp)**(-0.35)
endif
    ratraw(iR016) = term * narate

!! R017: He + e --> HE+ + e + e
if(ctemp .lt. 2800.0) then
    term = 0.0e0
else
    term = exp(-4.409864886e1 + 2.391596563e1*lnt1 - 1.07532302e1*lnt2 &
         +3.05803875e0*lnt3 - 5.68511890e-1*lnt4 + 6.79539123e-2*lnt5 &
         -5.00905610e-3*lnt6 + 2.06723616e-4*lnt7 - 3.64916141e-6*lnt8)
endif
  ratraw(iR017) = term * naterm


!! R018: HEP + e --> HEPP + e + e
if(ctemp .lt. 5500.0) then
    term = 0.0e0
else
     term = exp(-6.87104099e1 + 4.393347633e1*lnt1 - 1.84806699e1*lnt2 &
         +4.70162649e0*lnt3 - 7.6924663e-1*lnt4 + 8.113042e-2*lnt5  &
         -5.32402063e-3*lnt6 + 1.97570531e-4*lnt7 - 3.16558106e-6*lnt8)
endif
  ratraw(iR018) = term * naterm


!!R019: Hep + e --> He + gamma
  term = 0.32*(1.0e-11*(ctemp)**(-0.5)*(11.19 - 1.676*lt1 - 0.2852*lt2 + 0.04433*lt3)) !! This is case B
  term = term + 0.68*(1.0e-11*(ctemp)**(-0.5)*(12.72 - 1.615*lt1 - 0.3162*lt2 + 0.0493*lt3)) !! This is case A,
  term = term + 1.9e-3*ctemp**(-1.50)*exp(-473421.0/ctemp)*(1.0+0.3*exp(-94684.0/ctemp)) 
  ratraw(iR019) = term * narate
  
!!R020: Hepp + e --> Hep + gam
if(pchem_ccCase .eq. 0 ) then
    term = 2.538e-13*(1262456.0/ctemp)**1.503*(1.0+(2418500.0/ctemp)**0.407)**(-1.923) !! CASE A
else if (pchem_ccCase .eq. 1) then
    term = 5.506e-14*(1262456.0/ctemp)**1.500*(1.0+(460752.0/ctemp)**0.407)**(-2.242) !! CASE B
else
    term = 5.506e-14*(1262456.0/ctemp)**1.500*(1.0+(460752.0/ctemp)**0.407)**(-2.242) !! CASE B
endif
  ratraw(iR020) = term * narate  

!!R021: H2P + HM --> H2 + H
   term = 1.4e-7*t3**(-0.5)
   ratraw(iR021) = term * narate

!!R022: H2P + HM --> H + H + H
   term = 1.4e-7*t3**(-0.5)
   ratraw(iR022) = term * narate

!!R023: H2 + e --> H + HM
if(ctemp .lt. 500.0) then
   term = 0.0e0
else
   term = 2.7e-8*ctemp**(-1.27)*exp(-43000.0/ctemp)
endif
   ratraw(iR023) = term * narate

!!R024: H2 + HeP --> He + H + HP
   term = 3.7e-14*exp(35.0/ctemp)
   ratraw(iR024) = term * narate

!!R025: H2 + HeP --> HE + H2P
   term = 7.2e-15
   ratraw(iR025) = term * narate

!!R026: H + HEP --> HE + HP + GAM
   term = 1.2e-15*t3**(0.25)
   ratraw(iR026) = term * narate

!!R027: He + HP --> H + HeP
if(ctemp .gt. 10000.0) then
   term = 4.0e-37*ctemp**(4.74)
else if(ctemp .lt. 1500.0) then
   term = 0.0e0
else
   term = 1.26e-9*ctemp**(-0.75)*exp(-127500.0/ctemp)
endif
   ratraw(iR027) = term * narate

!!R028: HeP + HM --> He + H
   term = 2.32e-7*t3**(-0.52)*exp(ctemp/22400.0)
   ratraw(iR028) = term * narate

!!R029: HM + He --> H + He + e
if(ctemp .lt. 240.0) then
   term = 0.0e0
else
   term = 4.1e-17*ctemp**2.0*exp(-19870.0/ctemp)
endif
   ratraw(iR029) = term * narate

!!R030: H + H + H --> H2 + H
!if(ctemp .gt. 300.0) then
!   term = 3.9e-30*ctemp**(-1.0)
!else
!   term = 1.14e-31*ctemp**(-0.38)
!endif
   term = 1.44e-26*ctemp**(-1.54)  
   ratraw(iR030) = term * narate * narate !!den 

!!R031: H + H + H2 --> H2 + H2
   ratraw(iR031) = ratraw(iR030)*0.125

!!R032: H + H + He --> H2 + He
   term = 6.9e-32*ctemp**(-0.4)
   ratraw(iR032) = term * narate * narate  !!* den

!!R033: DP + e --> D + gamma
  ratraw(iR033) = ratraw(iR013)

!!R034: D + HP --> H + DP
   if(ctemp .gt. 2.0e5) then
      term = 3.44e-10*ctemp**(0.35)
   else
      term = 2.0e-10*ctemp**(0.402)*exp(-37.1/ctemp)-3.31e-17*ctemp**(1.48)
   endif
   ratraw(iR034) = term * narate

!!R035: H + DP --> D + HP
   term = 2.06e-10*ctemp**(0.396)*exp(-33.0/ctemp)+2.03e-9*ctemp**(-0.332)
   ratraw(iR035) = term * narate

!!R036: H + D --> HD + gam
if(ctemp .gt. 200.0) then
    term = 1.0e-25*exp(507.207-370.889*log(ctemp) + 104.854*log(ctemp)**2.0 &
           -14.4192*log(ctemp)**3.0 + 0.971469*log(ctemp)**4.0 - &
            0.0258076*log(ctemp)**5.0)
elseif (ctemp .lt. 10.0) then
    term = 0.0
else
    term = 1.0e-25*(2.80202 - 6.63697*log(ctemp) + 4.75619*log(ctemp)**2.0 &
           -1.39325*log(ctemp)**3.0 + 0.178259*log(ctemp)**4.0 &
           -0.00817097*log(ctemp)**5.0)
endif
    ratraw(iR036) = term * narate

!!R037: H2 + D --> HD + H
if(ctemp .gt. 2000.0) then
   term = 3.17e-10*exp(-5207.0/ctemp)
else
   term = (-56.4737 + 5.888886*lt1 + 7.19692*lt2 + 2.25069*lt3 - 2.16903*lt4 + 0.317887*lt5)
   term = 10**term
endif
   ratraw(iR037) = term * narate

!!R038: HDP + H --> HD + HP
   ratraw(iR038) = ratraw(iR004)

!!R039: H2 + DP --> HD + HP
   term = 4.17e-10 + 8.46e-10*lt1 - 1.37e-10*lt2
   ratraw(iR039) = term * narate

!!R040: HD + H --> H2 + D
if(ctemp .gt. 200.0) then
   term = 5.25e-11*exp(-4430.0/ctemp + 173900.0/ctemp**2)
else if(ctemp .lt. 50.0) then
   term = 0.0e0
else
   term = 5.25e-11*exp(-4430.0/ctemp)
endif
   ratraw(iR040) = term * narate

!!R041: HD + HP --> H2 + D+
   term = 1.1e-9*exp(-488.0/ctemp)
   ratraw(iR041) = term * narate

!!R042: D + HP --> HDP + gam 
   term = 3.9e-19*t3**(1.8)*exp(20.0/ctemp)
   ratraw(iR042) = term * narate

!!R043: H + DP --> HDP + gam
   term = 3.9e-19*t3**(1.8)*exp(20.0/ctemp)
   ratraw(iR043) = term * narate

!!R044: HDP + e --> H + D
   term = 7.2e-8*ctemp**(-0.5)
   ratraw(iR044) = term * narate

!! R045: D + e --> DP + e + e
   ratraw(iR045) = ratraw(iR012)

!!R046: D + HEP --> HE + DP + gam
   term = 1.1e-15*t3**(0.25)
   ratraw(iR046) = term * narate

!!R047: He + DP --> D + HeP
if(ctemp .gt. 10000.0) then
   term = 5.9e-37*ctemp**(4.74)
elseif(ctemp .lt. 1500.0) then
   term = 0.0e0
else
   term = 1.85e-9*ctemp**(-0.75)*exp(-127500.0/ctemp)
endif
   ratraw(iR047) = term * narate

!!R048: H2P + D --> HDP + H
   term = 1.07e-9*t3**(0.062)*exp(-ctemp/41400.0)
   ratraw(iR048) = term * narate

!!R049: D + HDP --> HD + DP
   term = 6.4e-10
   ratraw(iR049) = term * narate

!!R050: HDP + H --> H2P + D
   term = 1.0e-9*exp(-154.0/ctemp)
   ratraw(iR050) = term * narate

!!R051: D + e --> DM + gam
   ratraw(iR051) = ratraw(iR001)

!!R052: H + DM --> D + HM
   term = 6.4e-9*t3**(0.41)
   ratraw(iR052) = term * narate

!!R053: D + HM --> D + DM
   term = 6.4e-9*t3**(0.41)
   ratraw(iR053) = term * narate

!!R054: D + HM --> HD + e
   !!term = 1.5e-9*t3**(-0.1)
   ratraw(iR054) = 0.5*ratraw(iR002)

!!R055: H + DM --> HD + e
   !!term = 1.5e-9*t3**(-0.1)
   ratraw(iR055) = 0.5*ratraw(iR002)

!!R056: D + DM --> D2 + e
   !!term = 1.6e-9*t3**(-0.1)
   ratraw(iR056) = ratraw(iR002)

!!R057: HD + e --> DM + H
if(ctemp .lt. 500.0) then
   term = 0.0e0
else
   term = 1.35e-9*ctemp**(-1.27)*exp(-43000.0/ctemp)
endif
   ratraw(iR057) = term * narate

!!R058: HD + e --> D + HM
   term = 1.35e-9*ctemp**(-1.27)*exp(-43000.0/ctemp)
   ratraw(iR058) = term * narate

!!R059: D2 + e --> D + DM
   term = 6.7e-11*ctemp**(-1.27)*exp(-43000.0/ctemp)
   ratraw(iR059) = term * narate

!!R060: HP + DM --> HDP + e
   term = 1.1e-9*t3**(-0.4)
   ratraw(iR060) = term * narate

!!R061: DP + HM --> HDP + e
   term = 1.1e-9*t3**(-0.4)
   ratraw(iR061) = term * narate

!!R062: DP + DM --> D2P + e
   term = 1.3e-9*t3**(-0.4)
   ratraw(iR062) = term * narate

!!R063: DM + e --> D + e +e
   ratraw(iR063) = ratraw(iR014)

!!R064: DM + H --> D + H +e
   ratraw(iR064) = ratraw(iR015)

!!R065: DM + He --> D + He + e
if(ctemp .lt. 250.0) then
   term = 0.0e0
else
   term = 1.5e-17*ctemp**2.0*exp(-19870.0/ctemp)
endif
   ratraw(iR065) = term * narate

!!R066: DP + HM --> D + H
   ratraw(iR066) = ratraw(iR005)

!!R067: HP + DM --> D + H
   ratraw(iR067) = ratraw(iR005)

!!R068: DP + DM --> D + D
   ratraw(iR068) = ratraw(iR005)

!!R069: H2P + DM --> H2 + D
   term = 1.7e-7*t3**(-0.5)
   ratraw(iR069) = term * narate

!!R070: H2P + DM --> H + H + D
   term = 1.7e-7*t3**(-0.5)
   ratraw(iR070) = term * narate

!!R071: HDP + HM --> HD + H
   term = 1.5e-7*t3**(-0.5)
   ratraw(iR071) = term * narate

!!R072: HDP + HM --> D + H + H
   term = 1.5e-7*t3**(-0.5)
   ratraw(iR072) = term * narate

!!R073: HDP + DM --> HD + D
   term = 1.9e-7*t3**(-0.5)
   ratraw(iR073) = term * narate

!!R074: HDP + DM --> D + H + D
   term = 1.9e-7*t3**(-0.5)
   ratraw(iR074) = term * narate

!!R075: D2P + HM --> D2 + H
   term = 1.5e-7*t3**(-0.5)
   ratraw(iR075) = term * narate

!!R076: D2P + HM --> D + D + H
   term = 1.5e-7*t3**(-0.5)
   ratraw(iR076) = term * narate

!!R077: D2P + DM --> D2 + D
   term = 2.0e-7*t3**(-0.5)
   ratraw(iR077) = term * narate

!!R078: D2P + DM --> D + D + D
   term = 2.0e-7*t3**(-0.5)
   ratraw(iR078) = term * narate

!!R079: HEP + DM --> He + D
   term = 3.03e-7*t3**(-0.52)*exp(ctemp/22400.0)
   ratraw(iR079) = term * narate

!!R080: D + DP --> D2P + gam
   term = 1.9e-19*t3**(1.8)*exp(20.0/ctemp)
   ratraw(iR080) = term * narate

!!R081: D + H2P --> H2 + DP
   term = 6.4e-10
   ratraw(iR081) = term * narate

!!R082: H2P + D --> HD + HP
   term = 1.0e-9
   ratraw(iR082) = term * narate

!!R083: HDP + H --> H2 + DP
   term = 1.0e-9
   ratraw(iR083) = term * narate

!!R084: HDP + D --> D2P + H 
   term = 1.0e-9
   ratraw(iR084) = term * narate

!!R085: HDP + D --> D2 + HP
   term = 1.0e-9
   ratraw(iR085) = term * narate

!!R086: D + D2P --> D2 + DP
   term = 6.4e-10
   ratraw(iR086) = term * narate

!!R087: H + D2P --> D2 + HP
   term = 6.4e-10
   ratraw(iR087) = term * narate

!!R088: D2P + H --> HDP + D
   term = 1.0e-9*exp(-472.0/ctemp)
   ratraw(iR088) = term * narate

!!R089: D2P + H --> HD + DP
   term = 1.0e-9
   ratraw(iR089) = term * narate

!!R090: H2 + DP --> D + H2P
   ratraw(iR090) = ratraw(iR007)


!!R091: H2 + DP --> HDP + H
if(ctemp .lt. 230.0) then
   term = 0.0e0
else
   term = (1.04e-9 + 9.52e-9*(ctemp/10000.0) - 1.81e-9*(ctemp/10000.0)**2)*exp(-21000.0/ctemp)
endif
   ratraw(iR091) = term * narate

!!R092: HD + HP --> H + HDP
   ratraw(iR092) = ratraw(iR007)

!!R093: HD + HP --> H2P + D
if(ctemp .lt. 230.0) then
   term = 0.0e0
else
   term = 1.0e-9*exp(-21600.0/ctemp)
endif
   ratraw(iR093) = term * narate

!!R094: HD + DP --> D + HDP
   ratraw(iR094) = ratraw(iR007)

!!R095: HD + DP --> D2 + HP
   term = 1.0e-9
   ratraw(iR095) = term * narate

!!R096: HD + DP --> D2P + H
   term = (3.54e-9 + 7.50e-10*(ctemp/10000.0) - 2.92e-10*(ctemp/10000.0)**2)*exp(-21100.0/ctemp)
   ratraw(iR096) = term * narate

!!R097: D2 + HP --> HD + DP
   term = 2.1e-9*exp(-491.0/ctemp)
   ratraw(iR097) = term * narate

!!R098: D2 + HP --> HDP + D
   term = (5.18e-11 + 3.05e-9*(ctemp/10000.0) - 5.42e-10*(ctemp/10000.0)**2)*exp(-20100.0/ctemp)
   ratraw(iR098) = term * narate

!!R099: D2 + HP --> H + D2P
   ratraw(iR099) = ratraw(iR007)

!!R100: D2 + DP --> D2P + D
   ratraw(iR100) = ratraw(iR007)

!!R101: HD + HEP --> HE + HDP
   term = 7.2e-15
   ratraw(iR101) = term * narate

!!R102: HD + HEP --> HE + HP + D
   term = 1.85e-14*exp(35.0/ctemp)
   ratraw(iR102) = term * narate

!!R103: HD + HEP --> HE + HP + D
   term = 1.85e-14*exp(35.0/ctemp)
   ratraw(iR103) = term * narate

!!R104: D2 + HeP --> HE + D2P
   term = 2.5e-14
   ratraw(iR104) = term * narate

!!R105: D2 + HeP --> He + DP + D
   term = 1.1e-13*t3**(-0.24)
   ratraw(iR105) = term * narate

!!R106: HD + D --> D2 + H
   term = 1.15e-11*exp(-3220.0/ctemp)
   ratraw(iR106) = term * narate

!!R107: D2 + H --> HD + D
if(ctemp .gt. 2200.0) then
   term = 2.67e-10*exp(-5945.0/ctemp)
else
   term = -86.1558 + 4.53978*lt1 + 33.5707*lt2 - 13.0449*lt3 + 1.22017*lt4 + 0.0482453*lt5
   term = 10**term
endif
   ratraw(iR107) = term * narate

!!R108: HD + H --> H + D +H
   ratraw(iR108) = ratraw(iR009)

!!R109: HD + H2 --> H + D + H2
   ratraw(iR109) = ratraw(iR010)

!!R110: HD + He --> H + D + He
   ratraw(iR110) = ratraw(iR011)

!!R111: HD + e --> H + D + e
!   term = 0.0e0
!   term = nv*log10(5.09e-9*ctemp**(0.128)*exp(-103258.0/ctemp))
!   term = term + nlte*log10(1.04e-9*ctemp**(0.218)*exp(-53070.7/ctemp))
!   term = 10**(term)
if(pchem_rcCase .eq. 1) then
    term = 5.09e-9*ctemp**(0.128)*exp(-103258.0/ctemp)
else
    term = 1.04e-9*ctemp**(0.218)*exp(-53070.7/ctemp)
endif
   ratraw(iR111) = term * narate

!!R112: D2 + H --> D + D + H
   ratraw(iR112) = ratraw(iR009)

!!R113: D2 + H2 --> D + D + H2
   ratraw(iR113) = ratraw(iR010)

!!R114: D2 + He --> D + D + He
   ratraw(iR114) = ratraw(iR011)

!!R115: D2 + e --> D + D + e
!   term = 0.0e0
!   term = nv*log10(8.24e-9*ctemp**(0.126)*exp(-105388.0/ctemp))
!   term = term + nlte*log10(2.75e-9*ctemp**(0.163)*exp(-53339.7/ctemp))
!   term = 10**(term)
if(pchem_rcCase .eq. 1) then
   term = 8.24e-9*ctemp**(0.126)*exp(-105388.0/ctemp)
else
   term = 2.75e-9*ctemp**(0.163)*exp(-53339.7/ctemp)
endif
   ratraw(iR115) = term * narate

!!R116: HM + gam --> H + e
   term = 1.36e-11*pchem_j21
   ratraw(iR116) = term 

!!R117: DM + gam --> D + e
   term = 1.36e-11*pchem_j21
   ratraw(iR117) = term 

!!R118: H2P + gam --> H + HP
   term = 4.11e-12*pchem_j21
   ratraw(iR118) = term

!!R119: HDP + gam --> H + DP
   term = 2.05e-12*pchem_j21
   ratraw(iR119) = term

!!R120: HDP + gam --> D + HP
   term = 2.05e-12*pchem_j21
   ratraw(iR120) = term

!!R121: D2P + gam --> D + DP
   term = 4.11e-12*pchem_j21
   ratraw(iR121) = term

!!R122: H2 + gam --> H + H
   term = 1.3e-12*pchem_fshh2*pchem_j21
   ratraw(iR122) = term

!!R123: HD + gam --> H + D
   term = 1.45e-12*pchem_fshhd*pchem_j21
   ratraw(iR123) = term

!!R124: D2 + gam --> D + D
   term = 1.3e-12*pchem_j21
   ratraw(iR124) = term


!!Try to look at rates
!!  open(10,file="rates", access='append')
!!  print *, 'In pchem_networkRates'
!!  print *, 'ctemp: ', ctemp, '  denrate: ', denrate  

!!  do i = 1,111
!!     write(10,"(e13.5)",advance="no") zz(i) 
!!  enddo
 
!!     write(10,"(e8.3)", advance="no") zz
!!     write(10,"(A)", advance="no") ' '
!!     write(10,"(f10.3)",advance="no") ctemp
!!    write(10,"(A)") ' '
!!  close(10)




  return

  end subroutine pchem_networkRates
