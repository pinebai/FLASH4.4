!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_prepareNewGravityAccel
!!
!! NAME
!!
!!  hy_uhd_prepareNewGravityAccel
!!
!! SYNOPSIS
!!
!!  call hy_uhd_prepareNewGravityAccel(integer, INTENT(IN)  :: blockCount,
!!                                     integer, INTENT(IN), dimension(blockCount)  :: blocklist,
!!                                     logical(in) :: gcMaskLogged)
!!
!! DESCRIPTION
!!
!!  Call Gravity_potentialListOfBlocks to update the gravitational potential
!!  stored in UNK, for configurations with self-gravity.
!!
!!  Also perform ancillary actions needed in some configurations,
!!  in particular related to sink particles if those are being used.
!!
!!  Functionally, this prepares UNK and other simulation state so that
!!  subsequent matching calls to Gravity_accelOneRow return the proper
!!  gravitational accelerations acting on the fluid (or, more general,
!!  the mass distribution described by the "dens" variable of UNK).
!!
!!  The word "matching" in the above description means, in particular,
!!  agreement on which variable to use for storing the potential (whether "gpot",
!!  "gpoh", or "gpol"), as well as other auxiliary variables if needed.
!!  The reason for having this wrapper in the Hydro implementation
!!  is to achieve this agreement (in cooperation with the hy_uhd_putGravityUnsplit
!!  wrapper).
!!
!! ARGUMENTS
!!
!!   blockCount : block count, passed on to Gravity_potentialListOfBlocks
!!
!!   blocklist : list of blocks, passed on to Gravity_potentialListOfBlocks
!!
!!   gcMaskLogged : a flag to control logging of the guard cell mask for the
!!                  Grid_fillGuardCells call
!!
!! SEE ALSO
!!
!!  Gravity_potentialListOfBlocks
!!  Gravity_accelOneRow
!!***

#define DEBUG_GRID_GCMASK

#include "Flash.h"
#include "constants.h"

subroutine hy_uhd_prepareNewGravityAccel(blockCount,blockList,gcMaskLogged)

  use Grid_interface, ONLY : Grid_fillGuardCells,    &
                             Grid_addToVar
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Gravity_interface, ONLY : Gravity_potentialListOfBlocks

#if defined(GRAVITY) && defined(GPOT_VAR)
  use Particles_interface, ONLY : Particles_sinkAccelGasOnSinksAndSinksOnGas
#endif

  use Hydro_data, ONLY : hy_useGravity,hy_useCosmology, &
                         hy_gpotAlreadyUpToDate, &
                         hy_gpotVar,hy_extraAccelVars

  implicit none

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) ::  blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  logical,intent(in) :: gcMaskLogged
  !! -----------------------------------------------------

  logical :: gcMask(NUNK_VARS)

#ifdef GPOT_VAR
  if (hy_useGravity) then

     gcMask = .FALSE.

     !! Gravity calculation at n+1 by calling Poisson solver
     if (hy_gpotVar .LE. 0) then
        call Gravity_potentialListOfBlocks(blockCount, blockList)
        hy_gpotAlreadyUpToDate = .TRUE.
        !! Fill guardcells for only the gpot variable
        gcMask(GPOT_VAR) = .true.
     else
        call Gravity_potentialListOfBlocks(blockCount, blockList, potentialIndex=hy_gpotVar)
        !! Fill guardcells for only the gpoh or gpol variable
        hy_gpotAlreadyUpToDate = (hy_gpotVar==GPOT_VAR)
        gcMask(hy_gpotVar) = .true.
        if (hy_useCosmology) then
           call Particles_sinkAccelGasOnSinksAndSinksOnGas((/0,0,0/),accelVars=hy_extraAccelVars)
        else
#ifdef SGAX_VAR
           if (hy_extraAccelVars(1).NE.SGAX_VAR .AND. hy_extraAccelVars(1)>0) &
              call Grid_addToVar(SGAX_VAR, hy_extraAccelVars(1), 1.0, .TRUE.)
#endif
#ifdef SGAY_VAR
           if (hy_extraAccelVars(2).NE.SGAY_VAR .AND. hy_extraAccelVars(2)>0) &
              call Grid_addToVar(SGAY_VAR, hy_extraAccelVars(2), 1.0, .TRUE.)
#endif
#ifdef SGAZ_VAR
           if (hy_extraAccelVars(3).NE.SGAZ_VAR .AND. hy_extraAccelVars(3)>0) &
              call Grid_addToVar(SGAZ_VAR, hy_extraAccelVars(3), 1.0, .TRUE.)
#endif
        end if
     end if


#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant[Pot]')
     end if
#endif

     !! Guardcell filling for gpot
     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,&
          maskSize=NUNK_VARS,mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)
  endif
#endif

end subroutine hy_uhd_prepareNewGravityAccel
