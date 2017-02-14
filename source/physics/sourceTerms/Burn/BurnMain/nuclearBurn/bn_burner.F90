!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_burner
!!
!! NAME
!!  
!!  bn_burner
!!
!! SYNOPSIS
!! 
!!  bn_burner(
!!       real(in)  :: tstep,
!!       real(in)  :: temp,
!!       real(in)  :: density,
!!       real(in)  :: xIn(:),
!!       real(out) :: xOut(:),
!!       real(out) :: sdotRate)
!!
!!  
!! DESCRIPTION
!!
!!  Routine bn_burner drives the nuclear burning network  
!!     given time step tstep, temperature temp, density density, and 
!!     composition xIn, this routine returns the burned composition xOut
!!     and the energy generation rate sdotRate.
!!
!! ARGUMENTS
!!
!!  tstep:    time step 
!!  temp:     temperature
!!  density:  density
!!  xIn:      composition in
!!  xOut:     composition out
!!  sdotRate: energy generation rate
!!
!! PARAMETERS
!!
!!  bn_algebra  integer  specifies choice of bn_algebra ma28=1 or gift=2
!!  odeStepper  integer  specifies integration method bader-deuflhard=1 or rosenbrock=2
!!  useTable    integer  generate raw rates using table approximation=1 or direct method=2
!!
!! NOTES
!!
!!  Within the network setup process
!!  routine bn_network sets up the odes/derivatives 
!!  routine bn_networkDenseJakob sets up the dense jacobian 
!!  routine bn_networkSparseJakob sets up the sparse jacobian 
!!  routine bn_networkSparsePointers builds the nonzero locations for networkSparseJakob
!!  routine bn_networkRates generates the raw reaction rates
!!  routine bn_networkTable generates the raw rates using table interpolation
!!  routine bn_networkScreen applies screening corrections to the raw rates
!!  routine bn_networkWeak applies weak reaction forces (not found in every network, applied by 
!!         bn_networkScreen in networks that contain it)
!!
!!  Within the general network integration process, found in directory BurnIntegrate:
!!  routine bn_netIntegrate drives the integration of the odes
!!  this subroutine uses 8 routines called "steper" to drive the integration of 
!!  nuclear reaction networks with 
!!       2 choices of the time stepper and
!!       2 choices for the linear algebra.
!!  routine bn_baderMa28 drives a Bader-Deuflhard step with ma28 algebra
!!  routine bn_baderStepMa28 is a Bader-Deuflhard stepper with ma28 algebra
!!  routine bn_baderGift drives a Bader-Deuflhard step with gift algebra
!!  routine bn_baderStepGift is a Bader-Deuflhard stepper with ma28 algebra
!!  routine bn_rosenMa28 drives a Rosenbrock step with ma28 algebra
!!  routine bn_rosenGift drives a Rosenbrock step with gift algebra
!!  In addition, 
!!  routine bn_pzExtr does extrapolations for any of the Bader-Deuflhard bn_bader* routines
!!
!!  In this nuclearBurn directory, there are additional routines
!!  routine bn_azbar computes composition variables from the mass fractions; they are
!!              stored in Burn_dataEOS
!!  routine bn_ecapnuc computes neutron and proton electron capture rates
!!  routine bn_sneutx computes neutrino losses
!!  routine bn_ifermi12 does an inverse fermi integral for bn_sneutx
!!  routine bn_mazurek computes NI56 electron capture rates
!!  routine bn_mcord operates on the sparse matrix pointers
!!  routine bn_screen4 calculates screening factors for nuclear reaction rates in the
!!             weak, intermediate and strong regimes
!!
!!
!!***

subroutine bn_burner(tstep,temp,density,xIn,xOut,sdotRate)

  use Burn_dataEOS, ONLY:  btemp, den
  use Burn_data, ONLY: bn_algebra, bn_odeStepper, bn_useBurnTable, &
       & xmass, ymass, xoktot, xbadtot, bion, sneut, aion

  use bnIntegrate_interface, ONLY: bn_netIntegrate
  !  This are routine names to be passed as arguments.  Cannot be included
  !   in an EXTERNAL statement if you're going to use an interface
  use bnIntegrate_interface, ONLY: bn_baderMa28, bn_baderGift, &
       bn_rosenMa28, bn_rosenGift
  use bnNetwork_interface, ONLY: bn_network, bn_networkSparsePointers
  use bn_interface, ONLY: bn_azbar, bn_sneutx

  implicit none

#include "constants.h"
#include "Flash.h"

  ! See notes in bnNetwork_interface on why you can't use the interface for
  !  SparseJakob/DenseJakob
  external bn_networkSparseJakob, bn_networkDenseJakob
!  external bn_baderMa28, bn_baderGift, bn_rosenMa28, bn_rosenGift

  ! arguments
  real, intent(IN)                       :: tstep,temp,density
  real, intent(OUT)                      :: sdotRate
  real, intent(IN), dimension(NSPECIES)  :: xIn
  real, intent(OUT), dimension(NSPECIES) :: xOut

  !..local varaibles      
  integer         ::  i,k,nok,nbad,kount
  integer, parameter ::  tdim=10, iprint=0, nostore=0

  real            ::  stptry,stpmin,ys2(NSPECIES),                           &
       &                 ttime(tdim),elem(NSPECIES,tdim)

  real, parameter ::  avo  = 6.0221367e23,     ev2erg = 1.602e-12,           &
       &                  conv = ev2erg*1.0e6*avo, tol    = 1.0e-5,          &
       &                  beg = 0.0e0,             odescal = 1.0e-6

  !..set the the material and network variables
  btemp = temp
  den  = density

  do i=1,NSPECIES
     xmass(i) = xIn(i)
  enddo

  call bn_azbar()  !! generates ymass

  do i=1,NSPECIES
     ys2(i) = ymass(i)
  enddo


  !..get the reaction rates from a table or formula

  if (bn_useBurnTable) then
     call bn_networkTable
  else
     call bn_networkRates
  endif
  !! in most netoworks, the weak subroutine does not exist.
  !! so it is now called from Aprox19's bn_networkScreen
!!  call bn_networkWeak(ymass)
  call bn_networkScreen(ymass)


  !..set the time step variables for a single point burn
  stptry = tstep
  stpmin = stptry * 1.0e-20

  !..bader-deuflhard integration
  if (bn_odeStepper .eq. 1) then

     !..with ma28 bn_algebra
     if (bn_algebra .eq. 1) then
        call bn_netIntegrate(beg,stptry,stpmin,tstep,ys2,                        &
             &              tol,beg,nostore,                                    &
             &              ttime,elem,tdim,NSPECIES,tdim,NSPECIES,                 &
             &              nok,nbad,kount,odescal,iprint,                      &
             &              bn_network,bn_networkSparseJakob,bn_networkSparsePointers,bn_baderMa28)
        !!  last line is derivs,    jakob,   bjakob,  steper

        !..with gift bn_algebra
     else if (bn_algebra .eq. 2) then
        call bn_netIntegrate(beg,stptry,stpmin,tstep,ys2,                        &
             &              tol,beg,nostore,                                    &
             &              ttime,elem,tdim,NSPECIES,tdim,NSPECIES,                 &
             &              nok,nbad,kount,odescal,iprint,                      &
             &              bn_network,bn_networkDenseJakob,bn_networkSparsePointers,bn_baderGift)

     end if


     !..rosenbrock integration
  else if (bn_odeStepper .eq. 2) then

     !..with ma28 bn_algebra
     if (bn_algebra .eq. 1) then
        call bn_netIntegrate(beg,stptry,stpmin,tstep,ys2,                        &
             &              tol,beg,nostore,                                    &
             &              ttime,elem,tdim,NSPECIES,tdim,NSPECIES,                 &
             &              nok,nbad,kount,odescal,iprint,                      &
             &              bn_network,bn_networkSparseJakob,bn_networkSparsePointers,bn_rosenMa28)

        !..with gift bn_algebra
     else if (bn_algebra .eq. 2) then
        call bn_netIntegrate(beg,stptry,stpmin,tstep,ys2,                        &
             &              tol,beg,nostore,                                    &
             &              ttime,elem,tdim,NSPECIES,tdim,NSPECIES,                 &
             &              nok,nbad,kount,odescal,iprint,                      &
             &              bn_network,bn_networkDenseJakob,bn_networkSparsePointers,bn_rosenGift)

     end if
  end if


  xoktot  = xoktot + real(nok)
  xbadtot = xbadtot + real(nbad)


  !..the average energy generated over the time step
  sdotRate = 0.0e0
  do k=1,NSPECIES
     sdotRate = sdotRate + (ys2(k) - ymass(k)) * bion(k)
  enddo
  sdotRate = sdotRate * conv/tstep

  !..take into neutrino losses
  call bn_sneutx()
  sdotRate = sdotRate - sneut

  !..update the composition 
  do i=1,NSPECIES
     xOut(i)  = ys2(i) * aion(i)
  enddo
  return
end subroutine bn_burner




