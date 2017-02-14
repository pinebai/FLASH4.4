!!****f* source/Particles/Particles_unitTest
!!
!! NAME
!!
!!  Particles_unitTest
!!
!! SYNOPSIS
!!
!!  call Particles_unitTest(integer(in)     :: fileUnit,
!!                          logical(inout)  :: perfect  )
!!
!! DESCRIPTION
!!
!!  Entry point to testing code within the Particles unit.
!!
!!  Different implementations exist at different locations within the
!!  Particles code hierarchy, for testing different features of the unit.
!!  Test setups should provide appropriate code and configuaration directives
!!  in the part of the code that resides in the Simulation unit, in order
!!  to select the right implementation of Particles_unitTest, require any
!!  additional code units needed, provide initialization etc., and call
!!  Particles_unitTest.
!!
!!  See comments in specific implementations for more specific documentation.
!!
!! ARGUMENTS
!!
!!  fileUnit:  integer   number of file for output
!!  perfect:   logical   may be set to false if any test fails
!!
!! NOTES
!!
!!  Typically, Particles_unitTest will be called from a modified
!!  version of Driver_evolveFlash. That implementation of
!!  Driver_evolveFlash will typically
!!
!!   * initialize the actual argument associated with the dummy argument 'perfect'
!!     to TRUE, and
!!   * open a file for writing, whose unit number is passed as 'fileUnit';
!!     (Note that usually several files will be opened if the test is executed
!!     in parallel on several MPI tasks, each with a uniquely constructed filename
!!     that follows the pattern unitTest_NNNN.)
!!
!!   * call Particles_unitTest, possibly more than one time in a task depending
!!     on the test performed;
!!
!!   * consider the test successful (on the current MPI rank) iff perfect has not
!!     been set to FALSE by any invocation of Particles_unitTest, and
!!   * write a line containing the magic phrase "all results conformed with expected
!!     values" to the file given by fileUnit before closing it.
!!
!!  FlashTest considers a "unitTest" successful if it finds the magic phrase in
!!  each of the unitTest_NNNN files created by the test run.
!!
!! SEE ALSO
!!
!!  Driver_evolveFlash
!!
!!***

subroutine Particles_unitTest(fileUnit,perfect)

#include "Flash.h"
#include "constants.h"

  implicit none

  integer, intent(IN)           :: fileUnit ! Output to file
  logical, intent(INOUT)        :: perfect  ! Flag to indicate errors

  return
 
end subroutine Particles_unitTest
