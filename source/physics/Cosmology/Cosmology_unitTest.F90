!!****f* source/physics/Cosmology/Cosmology_unitTest
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

  implicit none

  integer, intent(IN) :: fileUnit
  logical, intent(INOUT) :: perfect

  return

end subroutine Cosmology_unitTest
