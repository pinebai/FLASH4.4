!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/Opacity_init
!!
!! NAME
!!
!!  Opacity_init
!!
!! SYNOPSIS
!!
!!  call Opacity_init ()
!!
!! DESCRIPTION
!!
!!  Initializes data for the OPAL opacity model using run time
!!  parameters.
!!
!! ARGUMENTS
!!
!!***
subroutine Opacity_init ()

  use Opacity_data

  use op_opalData,                 ONLY : op_opalMaxLowT

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Simulation_interface,        ONLY : Simulation_mapStrToInt

  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getMype

  use op_interface,                ONLY : op_initNumerics,                &
                                          op_initConstant,                &
                                          op_initConstcm2g,               &
                                          op_initTabulated,               &
                                          op_initIntegrate,               &
                                          op_setEnergyGroupBoundaries,    &
                                          op_writeOpacity

  use Logfile_interface,           ONLY : Logfile_stamp
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  implicit none

# include "constants.h"
# include "Opacity.h"
! # include "OpacityOPAL.h"
# include "Flash.h"
  
  character(len=MAX_STRING_LENGTH) :: massFracVarStr

  call Driver_getMype(GLOBAL_COMM,op_globalMe)

!
!
!    ...Check, if the opacity unit is needed at all.
!
!
  call RuntimeParameters_get ("useOpacity",   op_useOpacity)

  if (.not.op_useOpacity) then
       return
  end if

  call RuntimeParameters_get("rt_useMGD", op_useMGD)

  !!if(.not. op_useMGD) then
  !!   call Driver_abortFlash("[Opacity_init] useOpacity is .true. but rt_useMGD is .false.")
  !!end if

  call Timers_start("Opacity_init")

  call RuntimeParameters_get("op_emitConst", op_emitConst)
  call RuntimeParameters_get("op_absorbConst", op_absorbConst)
!
!
!    ...Get the write and low temperature action parameters.
!
!
  call RuntimeParameters_get ("opacity_writeOpacityInfo",   op_writeOpacityInfo)
!
!
!    ...Initialize the numeric subunit to have all numerical constants available.
!
!
  call op_initNumerics ()
!
!
!    ...Sets all the data related to species and elements.
!
!
  op_opalNumHydrogenAbundances = OP_OPAL_NUM_H_ABUNDANCES
  call RuntimeParameters_get("op_opalNumHydrogenAbundances",op_opalNumHydrogenAbundances)
!
!
!    ...Set internal parameters:
!
!          i) # of radiation energy groups
!
!
  call RuntimeParameters_get("rt_mgdNumGroups", op_nEnergyGroups)

  if (op_nEnergyGroups < 1) then
      call Driver_abortFlash ('[Opacity_init] ERROR: no energy groups found')
  end if

  call RuntimeParameters_get("op_hydrogenMassFracVar", massFracVarStr)
  call Simulation_mapStrToInt(massFracVarStr  , op_hydrogenMassFracVar,MAPBLOCK_UNK)
  call RuntimeParameters_get("op_hydrogenMassFrac"   , op_hydrogenMassFrac)

  ! The following live in op_opalData:
  call RuntimeParameters_get("op_opalMaxLowT"        , op_opalMaxLowT)
!
!
!    ...Allocate all needed arrays.
!
!
  allocate (op_absorptionKind           (1:op_opalNumHydrogenAbundances))
  allocate (op_emissionKind             (1:op_opalNumHydrogenAbundances))
  allocate (op_transportKind            (1:op_opalNumHydrogenAbundances))
  allocate (op_ionNumberDensity         (1:op_opalNumHydrogenAbundances))
  allocate (op_energyGroupBoundaries    (1:op_nEnergyGroups + 1))

  allocate (op_speciesOpacities         (ABSORPTION:TRANSPORT , 1:op_opalNumHydrogenAbundances))
!
!
!    ...Inititalize those arrays with zeros, which will contain the same data throughout the run.
!
!
  op_absorptionKind = 0
  op_emissionKind   = 0
  op_transportKind  = 0
!
!
!    ...Set the energy group boundaries.
!
!
  call op_setEnergyGroupBoundaries ()
!
!
!    ...Inititalize the other subunits. The low temperature subunit can only be inititalized
!       after the low temperature boundaries have been set, which in turn requires
!       initialization of the tabulated subunit, so the order matters.
!
!
  call op_initConstant  ()
  call op_initConstcm2g ()
  call op_initTabulated ()
  call op_initIntegrate ()

  ! Check for errors:
  

  if (op_opalNumHydrogenAbundances < 1) then
     ! If you are here, then something
     ! species is not specified properly. Make sure the OPAL tables
     ! are specified in the flash.par file.
     if (op_globalMe == MASTER_PE) then
        call Logfile_stamp(op_opalNumHydrogenAbundances, "[Opacity_init] op_opalNumHydrogenAbundances")
     end if
     call Driver_abortFlash("[opacity_init] No list of hydrogen abundances?")
  end if

!
!
!    ...Write the opacity unit data to file (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeOpacity ()
  end if
!
!
!    ...Ready!
!
!
  call Timers_stop("Opacity_init")
  return
end subroutine Opacity_init
