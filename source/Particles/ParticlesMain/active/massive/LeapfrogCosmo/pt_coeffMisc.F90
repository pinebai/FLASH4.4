!!****if* source/Particles/ParticlesMain/active/massive/LeapfrogCosmo/pt_coeffMisc
!!
!! NAME
!!
!!  pt_coeffMisc
!!
!! SYNOPSIS
!!
!!  function BNCOEFF( real, INTENT(in):: DT,
!!                    real, INTENT(in):: DT_OLD,
!!                    real, INTENT(in):: ALPHA,
!!                    real, INTENT(in):: DOTALPHA)
!!
!!   returns real B/C/DCOEFF
!!
!! ARGUMENTS
!!
!!    DT:
!!    DT_OLD: 
!!    ALPHA:
!!    DOTALPHA:   
!!
!!  
!! PARAMETERS
!!
!! DESCRIPTION
!!
!!  Calculation of coefficients "A,B,C" for Particles_advance.  See notes
!!       within each routine
!!  
!!  This version is the second-order leapfrog advancement for the active 
!!  submodule, including cosmological redshift terms.
!!
!!  The method used is that worked out by Scott Dodelson for COSMOS
!!  (2000, ApJ, 536, 122).
!!
!!***

!===============================================================================

!               Functions for the second-order variable-timestep leapfrog
!               velocity update.

function BNCOEFF (DT, DT_OLD, ALPHA, DOTALPHA)
  !       This function is the "B_n" function from
  !       equation 29 in SecondOrder.tex

  implicit none

  real :: BNCOEFF
  real, INTENT(in)  :: DT, DT_OLD, ALPHA, DOTALPHA
  real, parameter   :: onetwelfth = 1./12.
  real, parameter   :: onesixth   = 1./6.

  BNCOEFF = 1. - 0.5*ALPHA*DT + DT**2*(ALPHA**2 - DOTALPHA)*onesixth
  BNCOEFF = BNCOEFF * (1. - 0.5*ALPHA*DT_OLD + &
       DT_OLD**2*(ALPHA**2+2*DOTALPHA)*onetwelfth)

  return
end function BNCOEFF

!===============================================================================

function CNCOEFF (DT, DT_OLD, ALPHA, DOTALPHA)
  !       This function is the "C_n" function from
  !       equation 31 in SecondOrder.tex

  implicit none

  real :: CNCOEFF
  real, INTENT(in)   :: DT, DT_OLD, ALPHA, DOTALPHA
  real, parameter    :: onesixth   = 1./6.
  real, parameter    :: onethird   = 1./3.

  CNCOEFF = 0.5*DT + DT_OLD*onethird + DT**2/(6*DT_OLD) - &
       ALPHA*DT*(DT+DT_OLD)*onesixth

  return

end function CNCOEFF

!===============================================================================

function DNCOEFF (DT, DT_OLD, ALPHA, DOTALPHA)
  !       This function is the "D_n" function from
  !       equation 33 in SecondOrder.tex

  implicit none

  real :: DNCOEFF
  real, INTENT(in)  :: DT, DT_OLD, ALPHA, DOTALPHA
  real, parameter   :: onetwelfth = 1./12.

  DNCOEFF = (DT_OLD**2 - DT**2)/(6*DT_OLD) - ALPHA*DT_OLD*(DT+DT_OLD)*onetwelfth

  return
end function DNCOEFF

!===============================================================================
