!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_network
!!
!! NAME
!! 
!! bn_network
!!
!! SYNOPSIS
!!
!! bn_network(real, intent(IN)                :: tt,
!!            real, intent(IN), dimension(:)  :: y,
!!            real, intent(OUT), dimension(:) :: dydt
!!
!! DESCRIPTION
!!
!!  This routine sets up the derivatives for the Iso7 network.
!!
!!  Iso7, the second network installed in flash 1.5
!!  like the original iso13 network, but adds
!!  (a,p)(p,g) sequences through fe52 to ni56 and neutrino losses
!!  Routine bn_burner drives the aprox13 network
!!     given time step tstep, temperature temp, density density, and
!!     composition xIn, this routine returns the burned composition xOut
!!     and the energy generation rate sdotRate.
!!
!! ARGUMENTS
!!
!!   tt  -- never used, so don't know what it does
!!   y --
!!   dydt --
!!
!!***
subroutine bn_network(tt,y,dydt)   

  use Burn_dataEOS, ONLY: btemp, den
  use Burn_data
  use bn_dataIso7

#include "Flash.h"
 
  implicit none

  !.. 
  !..this routine sets up the system of ode's for the iso7 nuclear reactions.
  !..
  !..isotopes: he4, c12, o16, ne20, mg25, si28, ni56
  !.. 
  !..declare
  real, intent(IN) :: tt
  real, intent(INOUT), dimension(*)  :: y
  real, intent(OUT), dimension(*) :: dydt

  integer i
  real    t9,yeff_ca40,yeff_ti44,                   &
       &                 t9i,t932,t9i32,rsi2ni,rni2si


  !..positive definite mass fractions
  do i=1,NSPECIES
     y(i) = min(1.0e0,max(y(i),1.0e-30))
  enddo


  !..some common factors and branching ratios
  t9    = btemp * 1.0e-9
  t9i   = 1.0e0/t9
  t932  = t9 * sqrt(t9)
  t9i32 = 1.0e0/t932


  !..rsi2ni is the rate for silicon to nickel
  !..rni2si is the rate for nickel to silicon
  rsi2ni = 0.0e0
  rni2si = 0.0e0
  !      if (t9 .gt. 20.0) then
  if (t9 .gt. 2.5 .and. y(ic12)+y(io16) .le. 4.0e-3) then
     yeff_ca40 = (t9i32**3) * exp(239.42*t9i-74.741)
     yeff_ti44 = (t932**3) * exp(-274.12*t9i+74.914)
     rsi2ni = yeff_ca40*den**3 * y(ihe4)**3 * ratdum(ircaag)*y(isi28)
     rni2si = min(1.0e20,yeff_ti44*ratdum(irtiga)/(den**3*y(ihe4)**3))
  end if


  !..set up the system of ode's : 
  !..4he reactions  
  dydt(ihe4) =  3.0 * y(ic12) * ratdum(irg3a)                       &
       &            - 3.0 * y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)    &
       &            + y(io16) * ratdum(iroga)                             &
       &            - y(ic12) * y(ihe4) * ratdum(ircag)                   &
       &            + y(ic12) * y(ic12) * ratdum(ir1212)                  &
       &            + 0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)          &
       &            + y(io16) * y(io16) * ratdum(ir1616)                  &
       &            - y(io16) * y(ihe4) * ratdum(iroag)                   &
       &            + y(ine20) * ratdum(irnega)                           &
       &            + y(img24) * ratdum(irmgga)                           &
       &            - y(ine20) * y(ihe4) * ratdum(irneag)                 &
       &            + y(isi28) * ratdum(irsiga)                           &
       &            - y(img24) * y(ihe4) * ratdum(irmgag)                 &
       &            - 7.0e0 * rsi2ni * y(ihe4)                            &
       &            + 7.0e0 * rni2si * y(ini56)


  !..12c reactions  
  dydt(ic12) =   y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)         &
       &             - y(ic12) * ratdum(irg3a)                            &
       &             + y(io16) * ratdum(iroga)                            &
       &             - y(ic12) * y(ihe4) * ratdum(ircag)                  &
       &             - 2.0e0 * y(ic12) * y(ic12) * ratdum(ir1212)         &
       &             - y(ic12) * y(io16) * ratdum(ir1216)

  !..16o reactions  
  dydt(io16) = -y(io16) * ratdum(iroga)                             &
       &            + y(ic12) * y(ihe4) * ratdum(ircag)                   &
       &            - y(ic12) * y(io16) * ratdum(ir1216)                  &
       &            - 2.0e0 * y(io16) * y(io16) * ratdum(ir1616)          &
       &            - y(io16) * y(ihe4) * ratdum(iroag)                   &
       &            + y(ine20) * ratdum(irnega)

  !..20ne reactions 
  dydt(ine20) =  y(ic12) * y(ic12) * ratdum(ir1212)                 &
       &             + y(io16) * y(ihe4) * ratdum(iroag)                  &
       &             - y(ine20) * ratdum(irnega)                          &
       &             + y(img24) * ratdum(irmgga)                          &
       &             - y(ine20) * y(ihe4) * ratdum(irneag)


  !..24mg reactions  
  dydt(img24) =  0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)         &
       &             - y(img24) * ratdum(irmgga)                          &
       &             + y(ine20) * y(ihe4) * ratdum(irneag)                &
       &             + y(isi28) * ratdum(irsiga)                          &
       &             - y(img24) * y(ihe4) * ratdum(irmgag)

  !..28si reactions  
  dydt(isi28) =  0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)         &
       &             + y(io16) * y(io16) * ratdum(ir1616)                 &
       &             - y(isi28) * ratdum(irsiga)                          &
       &             + y(img24) * y(ihe4) * ratdum(irmgag)                &
       &             - rsi2ni * y(ihe4)                                   &
       &             + rni2si * y(ini56)

  !..ni56 reactions  
  dydt(ini56) =  rsi2ni * y(ihe4)                                   &
       &             - rni2si*y(ini56)

  return
end subroutine bn_network

