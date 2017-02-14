!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/Eos_getAbarZbar
!!
!! NAME
!!
!!  Eos_getAbarZbar
!!
!! SYNOPSIS
!!
!!  call Eos_getAbarZbar(OPTIONAL,real(IN)  :: solnVec(NUNK_VARS),
!!                       OPTIONAL,real(OUT) :: abar,
!!                       OPTIONAL,real(OUT) :: zbar,
!!                       OPTIONAL,real(OUT) :: sumY,
!!                       OPTIONAL,real(OUT) :: Ye,
!!                       OPTIONAL,real(IN)  :: massFrac(:) )
!!
!! DESCRIPTION
!!
!!  Special overriding version for avoiding the inclusion of the EOS unit. It sets the Abar
!!  and Zbar values equal to 1. This routine is needed to override the version in the Energy
!!  Deposition unit.
!!
!! ARGUMENTS
!!
!!   solnVec : optional - the solution vector for one cell
!!
!!   abar    : optional - if present, will be filled with the average atomic mass
!!   zbar    : optional - if present, will be filled with the average inonization level
!!
!!   sumY    : optional - if present, will be filled with the inverse of the
!!                        average atomic mass
!!   Ye      : optional - if present, will be filled with Ye, which is defined
!!                        according to   Ye = zbar / abar
!!              
!!   massFrac : this is an optional argument which may be used when there is more 
!!              than one species in the simulation, as an alternative to providing
!!              the complete solution vector in solnVec
!!   
!! NOTES
!!
!!  The argument list must match the argument list as declared in the Eos interface,
!!  even if those arguments are never needed here.
!!
!!***

subroutine Eos_getAbarZbar (solnVec,  abar,zbar,sumY,Ye,massFrac)

  implicit none

#include "Flash.h"

  real, OPTIONAL,dimension (NUNK_VARS),intent (IN)  :: solnVec
  real, OPTIONAL,                      intent (OUT) :: abar, zbar, Ye, sumY
  real, OPTIONAL,dimension (NSPECIES), intent (IN)  :: massFrac

  real :: dummy
!
!
!     ...Avoid compiler warnings about unused variables.
!
!
  if (present ( solnVec)) dummy = 0.0
  if (present (massFrac)) dummy = 0.0
!
!
!     ...Set the required values.
!
!
  if (present (abar)) abar = 1.0
  if (present (zbar)) zbar = 1.0
  if (present (  Ye)) Ye   = 1.0
  if (present (sumY)) sumY = 1.0
!
!
!     ...Ready!
!
!
  return
end subroutine Eos_getAbarZbar
