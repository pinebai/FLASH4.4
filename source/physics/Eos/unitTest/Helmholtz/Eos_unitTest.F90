!!****if* source/physics/Eos/unitTest/Helmholtz/Eos_unitTest
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

!!REORDER(4): solnData

subroutine Eos_unitTest(fileUnit, perfect)

  use Eos_interface, ONLY : Eos_wrapped, Eos
  use Grid_interface, ONLY :  Grid_getLocalNumBlks, &
       Grid_getBlkPtr, Grid_getBlkIndexLimits, Grid_releaseBlkPtr, &
       Grid_getBlkType, Grid_putRowData
  use IO_interface, ONLY : IO_writeCheckpoint
  use Eos_data, ONLY : eos_meshMe, eos_meshNumProcs
  use eos_testData, ONLY: eos_testPresModeStr, &
                          eos_testEintModeStr, &
                          eos_testTempModeStr, &
                          eos_testPresMode, &
                          eos_testEintMode, &
                          eos_testTempMode
  implicit none

# include "Eos.h"
# include "constants.h"
# include "Flash.h"

  integer, intent(in) :: fileUnit
  logical, intent(out) :: perfect
  integer :: localBlkCount, blockID
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, parameter :: tolerance = 1e-9
  real :: presErr, tempErr, eintErr

  real, pointer, dimension(:,:,:,:):: solnData
  logical:: test1,test2,test3,test4 !for a block
  logical:: test1allB,test2allB,test3allB,test4allB !for all blocks

  integer :: vecLen, blockOffset,  pres, dens, temp, e, n, m
  integer :: isize, jsize, ksize, i,j,k, nStartsAtOne
  real, dimension(:), allocatable :: eosData
  real, dimension(:), allocatable :: massFrac
  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
  real, allocatable, dimension(:,:,:,:) :: derivedVariables
  real, dimension(:,:,:), allocatable :: deriv1, deriv2

  character(len=7),pointer:: ap
  character(len=7),target :: a
  integer,parameter :: maxPrintPE = 20
  integer :: nodeType
  integer :: ib,ie,jb,je,kb,ke
  integer, dimension(3) :: startingPos, dataSize, startRow
     real presErr1, presErr2

  if (eos_meshNumProcs==1) then
     a = ''
     ap => a
  else
20   format(I6,':')
     write(a,20) eos_meshMe
     a = trim(adjustl(a))
     ap => a
  end if

! info for checkpoint output


  mask = .true.

  call Grid_getLocalNumBlks(localBlkCount)

  do blockID=1,localBlkCount
     call Grid_getBlkType(blockId,nodeType)
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
     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_TEMP or similar: Density, temperature in; energy, pressure out; mode=', &
          eos_testTempMode,eos_testTempModeStr
     if (eos_meshMe<maxPrintPE) then
        print*,ap,'The initialized extreme values are '
        print*,ap,'Initialized Density min',minval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Density max',maxval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Temperature min',minval(abs(solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Temperature max',maxval(abs(solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Pressure min',minval(solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Initialized Pressure max',maxval(solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))
     end if 

     solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)
    call Eos_wrapped(eos_testTempMode, blkLimits,blockID)

    !! Summarize results of MODE_DENS_TEMP (or similar) call
    if (eos_meshMe<maxPrintPE) then
        print*,ap,'The resulting extreme values are '
!!$        print*,ap,'Resulting Density min',minval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
!!$        print*,ap,'Resulting Density max',maxval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting Temperature min',minval(abs(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting Temperature max',maxval(abs(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting Pressure min',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Resulting Pressure max',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Resulting internal energy min',minval(abs(solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting internal energy max',maxval(abs(solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     end if

     ! Save the equilibrium values in O variables and (hopefully) write them out
     solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
     solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
     solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)

     call Grid_releaseBlkPtr(blockID,solnData)
  end do

  call IO_writeCheckpoint()   !! This is checkpoint 001


  test1allB = .TRUE.
  do blockID=1,localBlkCount
     call Grid_getBlkType(blockId,nodeType)
     call Grid_getBlkPtr(blockId,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ib=blkLimits(LOW,IAXIS)
     ie=blkLimits(HIGH,IAXIS)

     jb=blkLimits(LOW,JAXIS)
     je=blkLimits(HIGH,JAXIS)

     kb=blkLimits(LOW,KAXIS)
     ke=blkLimits(HIGH,KAXIS)

!  Testing density/energy in, temperature/pressure out
     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI or similar: Density, energy in; temperature, pressure out; mode=', &
          eos_testEintMode,eos_testEintModeStr
         !  Zero output variables
         ! solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=0  ! don't zero TEMP or eos_helm cannot converge in MODE_DENS_EI
         solnData(PRES_VAR,:,:,:)=0 
     call Grid_releaseBlkPtr(blockID,solnData)
     call Eos_wrapped(eos_testEintMode,blkLimits,blockID)

     call Grid_getBlkPtr(blockId,solnData)

     if (eos_meshMe<maxPrintPE) then !! Summarize results of MODE_DENS_EI (or similar) call
        print*,ap,'  Temperature min ',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Temperature max ',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Pressure min ',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Pressure max ',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     end if

     !! Calculate error from MODE_DENS_EI (or similar) call.
     tempErr = maxval(   &
             abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)- solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/ &
                 solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))  )
     if (eos_meshMe<maxPrintPE) then
        print*,ap,'  The calculated error in temperature is ',tempErr
     end if

     presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/&
          solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in pressure is ',presErr

     call Grid_releaseBlkPtr(blockID,solnData)

     test1 = (tolerance > tempErr)
     test1 = test1.and.(tolerance > presErr)
     if(test1) then
         if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_EI or similar is fine'
     else
        if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI or similar is BAD!!!'
        test1allB = .FALSE.
     endif
  end do

  call IO_writeCheckpoint()  !! This is checkpoint 002



  test2allB = .TRUE.
  do blockID=1,localBlkCount
     call Grid_getBlkType(blockId,nodeType)
     call Grid_getBlkPtr(blockId,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ib=blkLimits(LOW,IAXIS)
     ie=blkLimits(HIGH,IAXIS)

     jb=blkLimits(LOW,JAXIS)
     je=blkLimits(HIGH,JAXIS)

     kb=blkLimits(LOW,KAXIS)
     ke=blkLimits(HIGH,KAXIS)

! Testing density/pressure in, energy/temp out
     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_PRES or similar: Density, pressure in; energy, temperature out; mode=', &
          eos_testPresMode,eos_testPresModeStr
         solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0
         ! solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=0  ! don't zero TEMP or eos_helm cannot converge in any mode
     call Grid_releaseBlkPtr(blockID,solnData)
     call Eos_wrapped(eos_testPresMode,blkLimits,blockID)
     call Grid_getBlkPtr(blockId,solnData)

     !! Summarize results of MODE_DENS_PRES (or similar) call;
     !! calculate error from MODE_DENS_PRES (or similar) call.
     tempErr = maxval(abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
     eintErr = maxval(abs((solnData(EINT_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OENT_VAR,ib:ie,jb:je,kb:ke))/solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'  Energy min is',minval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'  Energy max is',maxval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'  Temperature min is',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'  Temperature max is',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in energy is ',eintErr
     if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in temperature is ',tempErr

     test2 = (tolerance > tempErr)
     test2 = test2.and.(tolerance > eintErr)
     if(test2) then
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_PRES or similar is fine'
     else
        if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_PRES or similar is BAD!!!'
        test2allB = .FALSE.
     endif

     call Grid_releaseBlkPtr(blockID,solnData)
  end do

  call IO_writeCheckpoint()   !! This is checkpoint 003



  test3allB = .TRUE.
  test4allB = .TRUE.
  do blockID=1,localBlkCount
     call Grid_getBlkType(blockId,nodeType)
     call Grid_getBlkPtr(blockId,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ib=blkLimits(LOW,IAXIS)
     ie=blkLimits(HIGH,IAXIS)

     jb=blkLimits(LOW,JAXIS)
     je=blkLimits(HIGH,JAXIS)

     kb=blkLimits(LOW,KAXIS)
     ke=blkLimits(HIGH,KAXIS)

     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print*,ap,'And now to verify the other solutions'

     ! Generate the initial conditions again
     solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=solnData(CPRS_VAR,ib:ie,jb:je,kb:ke)
     solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)
 
     !!  Do calculations in reverse order
     ! Density and pressure in, energy and temperature out
         !solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=0   ! don't zero TEMP or eos_helm cannot converge
         solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0 
     call Eos_wrapped(MODE_DENS_PRES, blkLimits,blockID)
     ! Now we have a "true"  temperature and internal energy; save them for comparison
     solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
     solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
     solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)
     

     ! Density and energy in, temperature and pressure out
        !! zero output values to make sure they're being calculated
       solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=0.0
       !solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=0.0   ! don't zero TEMP or eos_helm cannot converge 
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
     presErr1 = maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     presErr2 = maxval(solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print *,ap,'maxval PRES_VAR OPRS_VAR',presErr1,presErr2
     presErr = maxval(abs(&
          (solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/&
           solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     ! NOTE this ALWAYS comes out to zero. suspect  something wrong here....
     tempErr = maxval(abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
     test3 = (tolerance > tempErr)
     test3 = test3.and.(tolerance > presErr)
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in pressure from EI is ',presErr
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in temperature from EI is ',tempErr
     if(test3) then
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_EI is fine '
     else
        if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI is BAD!!!'
        test3allB = .FALSE.
     endif

         solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0
         solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=0 
     call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
     presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     eintErr = maxval(abs((solnData(EINT_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OENT_VAR,ib:ie,jb:je,kb:ke))/solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     test4 = (tolerance > presErr)
     test4 = test4.and.(tolerance > eintErr)
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in pressure is ',presErr
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in energy is ',eintErr
     if(test4) then
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_TEMP is fine'
     else
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_TEMP is BAD!!!!'
        test4allB = .FALSE.
     endif
     solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
     solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
     solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)
     call Grid_releaseBlkPtr(blockID,solnData)

     !! Finally, do a test of the derived variables just for exercise.....
     if (eos_meshMe<maxPrintPE) print *,ap,' Now testing the derived variables'

     !  Allocate the necessary arrays for an entire block of data
     isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
     jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
     ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

!!$     call Eos_wrapped(MODE_DENS_PRES,blkLimits,blockID)
     vecLen=isize

     allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
     allocate(eosData(vecLen*EOS_NUM))
     allocate(massFrac(vecLen*NSPECIES))

     ! Initialize them
     derivedVariables = 0.0
     vecLen=isize
     pres = (EOS_PRES-1)*vecLen
     dens = (EOS_DENS-1)*vecLen
     temp = (EOS_TEMP-1)*vecLen

     call Grid_getBlkPtr(blockID,solnData)
     
     ! Space and dimensions for scratch variables
     dataSize(1) = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
     dataSize(2) = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
     dataSize(3) = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
     startingPos = 1
     allocate(deriv1(dataSize(1),dataSize(2),dataSize(3)))
     allocate(deriv2(dataSize(1),dataSize(2),dataSize(3)))

     !! Get DENS and PRES to fill up input, also massFraction
     solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
           do i = 1,vecLen
              massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
                   solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
           end do
           
           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
           
           call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)

           
           startRow(1) = blkLimits(LOW,IAXIS)
           startRow(2) = j
           startRow(3) = k
           do e=EOS_VARS+1,EOS_NUM
              m = (e-1)*vecLen
              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
              if (e==EOS_DEA) &
                 call Grid_putRowData(blockID,SCRATCH_CTR,DRV1_SCRATCH_CENTER_VAR,EXTERIOR,IAXIS, &
                      startRow,eosData(m+1:m+vecLen),vecLen)
              if (e==EOS_DPT) &
                 call Grid_putRowData(blockID,SCRATCH_CTR,DRV2_SCRATCH_CENTER_VAR,EXTERIOR,IAXIS, &
                      startRow,eosData(m+1:m+vecLen),vecLen)

           end do

          !!Stuff a few test derivatives into scratch storage so you can see what they look like
          !!  Feel free to change the variable inserted
           do i= 1, vecLen
              deriv1(i,j-NGUARD,k-NGUARD) = derivedVariables(i,j-NGUARD,k-NGUARD,EOS_DEA)
              deriv2(i,j-NGUARD,k-NGUARD) = derivedVariables(i,j-NGUARD,k-NGUARD,EOS_DPT)
           end do
        end do
     end do
     
     !! Stuff a few test derivatives into scratch storage so you can see what they look like
     !!  Feel free to change the variable inserted
     !! This is an alternate way of outputting data.  But since LBR can't see anything, who can
     !!  tell which is working
     !!call Grid_putBlkData(blockID,SCRATCH_CTR,DRV1_SCRATCH_CENTER_VAR,INTERIOR,startingPos, &
     !!           deriv1,dataSize)
     !!call Grid_putBlkData(blockID,SCRATCH_CTR,DRV2_SCRATCH_CENTER_VAR,INTERIOR,startingPos, &
     !!           deriv2,dataSize)

     call Grid_releaseBlkPtr(blockID,solnData)

     deallocate(deriv1)
     deallocate(deriv2)

     
     if (All(derivedVariables .eq. 0.0)) then
       if (eos_meshMe<maxPrintPE) print*,ap,"No derived variables were set!"
     end if 
     if (ANY(derivedVariables .ne. 0.0)) then
      if (eos_meshMe<maxPrintPE) print*,ap,"Some derived variables were set."
     end if 
     deallocate(eosData)
     deallocate(massFrac)
     deallocate(derivedVariables)
     
  end do

  !! Output to get the derived variables
  call IO_writeCheckpoint()

  if (eos_meshMe.EQ.MASTER_PE) print*,'out of the loop'
  perfect = test1allB.and.test2allB.and.test3allB.and.test4allB
  if(perfect) then
     if (eos_meshMe<maxPrintPE) print*,ap,'SUCCESS all tests were fine'
  else
     if (eos_meshMe<maxPrintPE) print*,ap,'FAILURE some tests failed'
  end if
  return
end subroutine Eos_unitTest




