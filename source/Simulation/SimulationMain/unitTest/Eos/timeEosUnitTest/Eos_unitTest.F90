!!****if* source/Simulation/SimulationMain/unitTest/Eos/timeEosUnitTest/Eos_unitTest
!! NAME
!!
!!  Eos_unitTest
!! 
!! SYNOPSIS
!!
!!  call Eos_unitTest(integer(IN) :: fileUnit,
!!                    logical(OUT) :: perfect
!!
!! DESCRIPTION
!!
!! This function is the unit test for the Eos unit. It is invoked in
!! the setup unitTest/Eos. The Config file for Eos unit test setup
!! requests a few extra variables in the main grid data structure for
!! Grid scope temporary storage. The Simulation_initBlock of the Eos
!! unit test initializes density in the right place for the DENS_VAR
!! variable (see Flash.h for DENS_VAR, TEMP_VAR etc definitions), and
!! temperature and pressure in the extra storage space CTMP_VAR
!! and CPRS_VAR. The physical quantities at this point are not in
!! thermal equilibrium. 
!!
!! The Eos_unit test starts by copying the initialized
!! temperature into the TEMP_VAR location and calling the
!! Eos_wrapped function with eosMode = MODE_DENS_TEMP, where
!! density and temperature are given and pressure and energy are
!! calculated. Now PRES_VAR and EINT_VAR contain values of pressure
!! and internal energy that are in thermal equilibrium, and the pressure values
!! are not necessarily what was stored in the extra storage space
!! during intialization. 
!! 
!! At this point in time three quantities; temperature,
!! pressure and energy are saved in the extra storage requested by
!! the unitTest/Eos setup, say OTMP_VAR, OPRS_VAR and OENT_VAR. Now
!! the Eos_unitTest function calls Eos_wrapped with eosMode =
!! MODE_DENS_PRES, followed by eosMode= MODE_DENS_EI.  If the
!! newly calculated values of temperature, pressure and energy are
!! the same as those saved in OTMP_VAR, OPRS_VAR and OENT_VAR, then
!! we can conclude that the Eos is working in MODE_DENS_PRES and
!! MODE_DENS_EI modes. However, we still can't say anything about the
!! MODE_DENS_TEMP mode. So we repeat the process by copying CPRS_VAR
!! into PRES_VAR and calling Eos_wrapped with MODE_DENS_PRES. We
!! again save the calculated values in the extra storage and make two
!! more Eos_wrapped calls with the remaining two modes. This time if
!! the new and old values of variables compare, we can conclude that
!! MODE_DENS_TEMP works too, and hence the unit test is successful.
!!
!! A final test calculates the optional derivates by setting the mask
!! argument to Eos equal to true.  The values of these arguments are not
!! tested in any way.  The unitTest simply makes sure that they can be calculated
!! without NaNs or the like.  This Helmholtz version writes out two
!! randomly selected variables to the checkpoint files for
!! visualization.  It's kludgy, but the user can change which
!! derivatives they want written.  Sadly, xflash or visit cannot
!! yet view grid scratch variables.
!!
!!
!!  ARGUMENTS 
!!   
!!   fileUnit : unit number for file opened by the unitTest/Eos setup
!!              in which to write results of the test
!!
!!   perfect : indicates test ran without error is true.
!!
!!  PARAMETERS
!!
!!  eintSwitch  a rarely used switch which ensures that internal energy calculations 
!!        maintain sufficient precision. Important only if energyTotal is dominated 
!!        by energyKinetic.
!!
!!***


subroutine Eos_unitTest(fileUnit, perfect)
  use Eos_interface, ONLY : Eos_wrapped
  use Eos_data, ONLY : eos_meshMe
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkPtr, Grid_getBlkIndexLimits, Grid_releaseBlkPtr
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

# include "Eos.h"
# include "constants.h"
# include "Flash.h"
  include "Flash_mpi.h"
  integer, intent(in) :: fileUnit
  logical, intent(out) :: perfect
  integer :: localBlkCount, blockID
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, pointer, dimension(:,:,:,:):: solnData
  integer :: ib,ie,jb,je,kb,ke
  real :: startTime, endTime
  integer :: repeat, ierr
  integer :: num_eos_calls

  call RuntimeParameters_get("num_eos_calls", num_eos_calls)
  call Grid_getLocalNumBlks(localBlkCount)

  if (localBlkCount >= 1) then
     blockID=1
     call Grid_getBlkPtr(blockId,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     !! In Simulation_initBlock,
     !! temperature is initialized in CTMP_VAR and pressure is
     !! initialized in CPRS_VAR. We don't change these variables
     !! we copy them into the usual variable name as needed.

     ib=blkLimits(LOW,IAXIS)
     ie=blkLimits(HIGH,IAXIS)

     jb=blkLimits(LOW,JAXIS)
     je=blkLimits(HIGH,JAXIS)

     kb=blkLimits(LOW,KAXIS)
     ke=blkLimits(HIGH,KAXIS)

     ! Testing density/temperature in; energy/pressure out
     solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)

     if (eos_meshMe == 0) then
        write(6,"(a,i10,a)") " We will now call Eos_wrapped in MODE_DENS_TEMP", &
             num_eos_calls, " times"
     end if

     startTime = MPI_WTime()
     do repeat = 1, num_eos_calls
        call Eos_wrapped(MODE_DENS_TEMP, blkLimits,blockID)
     end do
     endTime = MPI_WTime()

     if (eos_meshMe == 0) then
        write(6,"(a,i10,a,f12.4,a)") " Completed ", num_eos_calls, &
             " EOS calls in ", endTime - startTime, " seconds"
     end if
     call Grid_releaseBlkPtr(blockID,solnData)
  end if

  perfect = .true.
  call MPI_Finalize(ierr)
  stop
     
end subroutine Eos_unitTest
