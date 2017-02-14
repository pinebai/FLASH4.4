!!****if* source/physics/Cosmology/unitTest/Cosmology_unitTest
!!
!! NAME
!!
!!  Cosmology_unitTest
!!
!! SYNOPSIS
!!
!!  Cosmology_unitTest(integer, intent(IN)  :: fileunit,
!!                     logical, intent(INOUT)  :: perfect)
!!
!! DESCRIPTION
!!
!!  This unit test checks the generation of cosmological redshifts and scaling 
!!  factors, by comparing them to an analytical solution, ensuring that the 
!!  Friedmann equation is solved reasonably accurately, and that the redshift
!!  generated and time reported from the redshift is also reasonably accurate.
!!
!!  Presently, the tests are run with a tolerance of 10e-6 vs the analytical
!!  solution.
!!
!!  The tests are run under these conditions:
!!  OmegaMatter = 1.0
!!  OmegaLambda = 0.0
!!  OmegaBaryon = 1.0
!!  OmegaRadiation = 0.0
!!  HubbleConstant(H0) = 1.62038e-18
!!
!!  The cosmological scaling factor is related to the time by the equation:
!!
!!   a(t) = (t/t0)**(2/3)
!!   where t0 = 2/(3 * H0)
!!   and is related to the cosmolgical redshift by the 
!!   equation z(t) = (1/a(t)) - 1
!! 
!!   The test case then iterates over time starting at the time given by a 
!!   redshift of 50 and generates a number of uniform dt steps.The test will 
!!   then check the agreement of scaling factor, redshift and time and see 
!!   if the scaling factor is solved correctly and ensure that nothing 
!!   pathological is going on in solving for the current redshift or 
!!   converting from the current redshift to time.
!!
!! ARGUMENTS
!!
!!
!!   fileunit : This is the file unit number passed in by Driver_evolveFlash
!!              This goest to the "unitTest_xxxx" output file
!!
!!   perfect : true if the test was completed within tolerances.  
!!
!! NOTES
!!  
!!  This test is presently run over 10000 steps running from z=50 to 0.
!!  Shrinking the timestep will increase the accuracy, unless the redshift 
!!  ever goes below zero.  The results of doing this in this code are 
!!  not well defined, and it will show a very large amount of error. 
!!
!!
!!***

subroutine Cosmology_unitTest (  fileUnit, perfect )
  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Cosmology_interface, ONLY : Cosmology_redshiftToTime, &
    Cosmology_solveFriedmannEqn, Cosmology_getRedshift
  use Cosmology_data, ONLY : csm_scaleFactor, csm_oldScaleFactor, &
                             csm_hubble, csm_meshMe

  implicit none
   
#include "constants.h" 

  integer, intent(IN) ::  fileUnit
  logical, intent(INOUT) :: perfect
  
  
  real :: tInitial, tFinal, dt, time
  real :: zInitial, zFinal
  real :: t0
  
  real :: testTime, testRedshift
  real :: analytical_scaleFactor, analytical_redshift
  real :: magErrScaleFactor, magErrRedshift, magErrTime

  real, parameter :: tol = 1.e-6
  integer, parameter :: count = 10000

  perfect = .true.
  t0 = 2.0 / (3.0 * csm_hubble)
  
  call RuntimeParameters_get("zInitial", zInitial)
  call RuntimeParameters_get("zFinal", zFinal)

  if (MASTER_PE == csm_meshMe) print*, "Z = ", zFinal, zInitial
  call Cosmology_redshiftToTime(zInitial, tInitial)
  call Cosmology_redshiftToTime(zFinal, tFinal)

  if(MASTER_PE == csm_meshMe) print *, "Times = ", tFinal, tInitial

  dt = (tFinal - tInitial) / real(count)
  if(MASTER_PE == csm_meshMe) print *, "dt = ", dt
  time = tInitial
  
  !open (UNIT=10, FILE="testOut")

  testRedshift = zInitial

  do while(time < tFinal)

     analytical_scaleFactor = (time / t0)**(2.0/3.0)
     analytical_redshift = 1.0/analytical_scaleFactor - 1.0

     call Cosmology_redshiftToTime(testRedshift, testTime)

     !check scaleFactor
     if(abs((analytical_scaleFactor - csm_scaleFactor)/analytical_scaleFactor) > MAX(TINY(1.), tol)) then
       magErrScaleFactor = abs((analytical_scaleFactor - csm_scaleFactor)/analytical_scaleFactor)
        if(MASTER_PE == csm_meshMe) print *, "[Cosmology_unitTest]: Mismatch in scaleFactor, at time ", time
        if(MASTER_PE == csm_meshMe) print *, "expected: ", analytical_scaleFactor
        if(MASTER_PE == csm_meshMe) print *, "found: ", csm_scaleFactor
        if(MASTER_PE == csm_meshMe) print *, "magnitude: ", magErrScaleFactor
        perfect = .false.
     end if
     !check redshift
     if(abs((analytical_redshift - testRedshift)/analytical_redshift) > MAX(TINY(1.), tol)) then
        magErrRedshift = abs((analytical_redshift - testRedshift)/analytical_redshift)
        if(MASTER_PE == csm_meshMe) print*, "[Cosmology_unitTest]: Mismatch in redshift, at time ", time
        if(MASTER_PE == csm_meshMe) print *, "expected: ", analytical_redshift
        if(MASTER_PE == csm_meshMe) print *, "found: ", testRedshift
        if(MASTER_PE == csm_meshMe) print *, "magnitude: ", magErrRedshift
        perfect = .false.
     end if
     !check time
     if(abs((time - testTime)/time) > MAX(TINY(1.), tol)) then
        magErrTime = abs((time - testTime) / time)
        if(MASTER_PE == csm_meshMe) print*, "[Cosmology_unitTest]: Mismatch in time, at time ", time
        if(MASTER_PE == csm_meshMe) print *, "expected: ", time
        if(MASTER_PE == csm_meshMe) print *, "found: ", testTime
        if(MASTER_PE == csm_meshMe) print *, "magnitude: ", magErrTime

        perfect = .false.
     end if
   
   !  write(10, '(6(ES23.15))') time, analytical_scaleFactor, csm_scaleFactor, analytical_redshift, testRedshift, testTime
     
     call Cosmology_solveFriedmannEqn(time+dt, dt)
     call Cosmology_getRedshift(testRedshift)
     
     time = time + dt
    
     !if(.not. perfect) return
  end do
    
  !close(10)
  
  return

end subroutine Cosmology_unitTest
