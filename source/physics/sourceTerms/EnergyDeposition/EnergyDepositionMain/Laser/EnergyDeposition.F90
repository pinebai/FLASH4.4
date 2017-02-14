!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/EnergyDeposition
!!
!! NAME
!!
!!  EnergyDeposition
!!
!! SYNOPSIS
!!
!!  call EnergyDeposition (integer, intent (in)           :: blockCount, 
!!                         integer, intent (in)           :: blockList (:), 
!!                         real,    intent (in)           :: timeStep,
!!                         real,    intent (in)           :: timeSimulation,
!!                         integer, intent (in), optional :: passSplitDriver)
!!
!! DESCRIPTION
!!
!!  Compute the energy deposited due to irradiation by a laser.
!!  This is the driver routine for computing the energy deposition
!!  during one timestep. It is assumed that the domain is of
!!  dimensions no larger than it takes for light to travel during
!!  this timestep. Hence we follow all rays through the complete
!!  domain during the timestep.
!!
!!  The code operates in a loop where a local routine computes the 
!!  energy deposition and traverses the rays through individual blocks. 
!!  The loop terminates when all the rays have been processed.
!!
!! ARGUMENTS
!!
!!  blockCount      : Number of blocks on current processor
!!  blockList       : All block ID numbers
!!  timeStep        : current timestep value
!!  timeSimulation  : current simulation time
!!  passSplitDriver : indicates first/second half of time step for split driver
!!
!! NOTES
!!          
!!  The current implementation assumes presence of one or more laser 
!!  beams and their path is computed using the geometric optics assumptions.
!!
!!***

#define DEBUG_GRID_GCMASK

subroutine EnergyDeposition (blockCount, blockList, timeStep, timeSimulation, passSplitDriver)

  use EnergyDeposition_data, ONLY : ed_beams,                      &
                                    ed_currentStepNumber,          &
                                    ed_depoDt,                     &
                                    ed_depoReuseMaxSteps,          &
                                    ed_depoVar,                    &
                                    ed_depoVarValid,               &
                                    ed_energyInTimestep,           &
                                    ed_energyOutTimestep,          &
                                    ed_gcMaskSize,                 &
                                    ed_gcMask,                     &   
                                    ed_globalComm,                 &
                                    ed_globalMe,                   &
                                    ed_gradOrder,                  &
                                    ed_irradVar,                   &
                                    ed_laserIONumberOfPositions,   &
                                    ed_laserIONumberOfRaysWritten, &
                                    ed_laserIORayFrequency,        &
                                    ed_laserIORayPositions,        &
                                    ed_laserIORayPower,            &
                                    ed_laserIORayTags,             &
                                    ed_laserIOWrite,               &
                                    ed_prevDt,                     &
                                    ed_printRays,                  &
                                    ed_rayCount,                   &
                                    ed_rays,                       &
                                    ed_savedDepoStepNumber,        &
                                    ed_useEnergyDeposition

  use Driver_interface,      ONLY : Driver_getNStep


  use Logfile_interface,     ONLY : Logfile_stampVarMask

  use Grid_interface,        ONLY : Grid_fillGuardCells,       &
                                    Grid_getBlkPtr,            &
                                    Grid_notifySolnDataUpdate, &
                                    Grid_releaseBlkPtr

  use Timers_interface,      ONLY : Timers_start, &
                                    Timers_stop

  use ed_interface,          ONLY : ed_createRays,        &
                                    ed_createRayTags,     &
                                    ed_checkReuseDepo,    &
                                    ed_laserIOCheckWrite, &
                                    ed_printEnergyStamp,  &
                                    ed_printRaysData,     &
                                    ed_saveRays,          &
                                    ed_traceRays,         &
                                    ed_updatePlasma,      &
                                    ed_updateRays

  use ed_commInterface,      ONLY : ed_commInitComm,          &
                                    ed_commProgressTransport

  use IO_interface,          ONLY : IO_endRayWrite,   &
                                    IO_startRayWrite, &
                                    IO_writeRays

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  include "Flash_mpi.h"

  integer, intent (in)           :: blockCount
  integer, intent (in)           :: blockList (1:blockCount)
  real,    intent (in)           :: timeStep
  real,    intent (in)           :: timeSimulation
  integer, intent (in), optional :: passSplitDriver

  logical :: activeRays
  logical :: createRays
  logical :: moveRays
  logical :: isTransportDone
  logical :: laserIsOn
  logical :: reuseDepo

  integer :: block, blockID

  real    :: timeStepScaleFactor, scaleFact

  real, pointer :: solnData (:,:,:,:)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

!
!
!     ...Check, if energy deposition is needed at all. If not, return at once.
!
!
  if (.not. ed_useEnergyDeposition) then
       return
  end if

  call Timers_start ("EnergyDeposition")
!
!
!     ...Get the current step number of the simulation.
!
!
  call Driver_getNStep (ed_currentStepNumber)
!
!
!     ...Set ed_laserIOWrite logical variable if it is time to do laser IO.
!        Begin laser IO if needed.
!
!
  call ed_laserIOCheckWrite (passSplitDriver)

  laserIsOn  = any (      (timeSimulation >= ed_beams (:) % pulseStartingTime) &
                    .and. (timeSimulation <= ed_beams (:) % pulseEndingTime)   )

  call ed_checkReuseDepo (reuseDepo,laserIsOn)

  if (ed_laserIOWrite) then
      call IO_startRayWrite ()
  end if

  createRays = (laserIsOn .AND. .NOT. reuseDepo)

  if (createRays) then
!
!
!     ...Fills all guard cells in the 'unk' array in all directions.
!
!
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(ed_gcMask, .FALSE., '[EnergyDeposition]', 'gcNeed')
     end if
#endif
     call Grid_fillGuardCells (CENTER,ALLDIR,minLayers=max(1,ed_gradOrder),&
          eosMode=MODE_EOS_NOP,                  &
          maskSize=ed_gcMaskSize, mask=ed_gcMask,&
          makeMaskConsistent=.TRUE.,             &
          doLogMask=.NOT.gcMaskLogged)
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        gcMaskLogged = .TRUE.
     end if
#endif
!
!
!     ...Initialize the deposition variable in 'unk' array.
!
!
     do block = 1,blockCount
        blockID = blockList (block)
        call Grid_getBlkPtr (blockID, solnData, CENTER)
        solnData (ed_depoVar,:,:,:) = 0.0
        if (ed_irradVar > 0) solnData (ed_irradVar,:,:,:) = 0.0
        call Grid_releaseBlkPtr (blockID, solnData, CENTER)
     end do
  end if
  
  if (.NOT. laserIsOn) then
!
!   Just zero out depo var for time intervals when the laser is off
!
     do block = 1,blockCount
        blockID = blockList (block)
        call Grid_getBlkPtr (blockID, solnData, CENTER)
        solnData (ed_depoVar,:,:,:) = 0.0
        if (ed_irradVar > 0) solnData (ed_irradVar,:,:,:) = 0.0
        call Grid_releaseBlkPtr (blockID, solnData, CENTER)
     end do  
  end if 
  
  
!
!     ...Create rays to be inserted into the domain if needed. If not,
!        exit the routine.
!
!
  if (createRays) then

!
!     ...Initialize laser energy accumulation variables.
!
      ed_energyInTimestep  = 0.0
      ed_energyOutTimestep = 0.0
      call Timers_start  ("Create Rays")
      call ed_createRays (blockCount, blockList, timeStep, timeSimulation)
      call Timers_stop   ("Create Rays")

      moveRays = .false.              ! because ray block ID's are valid at creation

      call Timers_start  ("ed_updateRays")
      call ed_updateRays (moveRays) 
      call Timers_stop   ("ed_updateRays")

  else if (reuseDepo) then

      if (ed_globalMe==MASTER_PE) then
         timeStepScaleFactor = timeStep / ed_prevDt
         ed_energyInTimestep  = timeStepScaleFactor * ed_energyInTimestep
         ed_energyOutTimestep = timeStepScaleFactor * ed_energyOutTimestep
      else
         ed_energyInTimestep  = 0.0 !non-master contribution are already included in master's accumulation
         ed_energyOutTimestep = 0.0
      end if

  else ! if (.NOT. laserIsOn) then
      call Timers_stop   ("EnergyDeposition")
      return
  endif


  IF (.NOT. reuseDepo) then
!
!
!    ...Create a unique tag for each ray.
!
!
     call ed_createRayTags (passSplitDriver)
!
!
!     ...Print the rays data created on all processors (if requested).
!
!
     if (ed_printRays) then
        call ed_printRaysData (ed_globalMe)
     end if
!
!
!     ....Now trace the path of the ray through the domain on a cell by cell
!         basis until all rays are processed.
!
!            ed_traceRays       : local routine to trace the path of the
!                                 current rays through cells within a block
!                                 and deposit their energy into the cells.
!
!            IO_writeRays       : perform an IO on the rays (if needed).
!
!            ed_updateRays      : local routine to update the rays info
!                                 after their previous tracing.
!
!            MPI_Allreduce      : inform all processors, if at least one
!                                 ray is still existing somewhere, in which
!                                 case we must still continue processing
!                                 the ray(s).
!
!
     call Timers_start ("Transport Rays")

     call ed_commInitComm ()

     moveRays   = .true.
     activeRays = .true.

     do while (activeRays)       
        !CD: The following timers are commented out because each MPI rank
        !can execute this do while loop a different number of times.
        !call Timers_start  ("ed_traceRays")
        call ed_traceRays  (timeStep)
        !call Timers_stop   ("ed_traceRays")
     
        !call Timers_start  ("ed_saveRays")
        call ed_saveRays   ()
        !call Timers_stop   ("ed_saveRays")

        if (ed_laserIOWrite) then
           call IO_writeRays (ed_laserIONumberOfRaysWritten, &
                ed_laserIORayTags,             &
                ed_laserIORayPositions,        &
                ed_laserIORayPower,            &
                ed_laserIONumberOfPositions    )
        end if

        !call Timers_start  ("ed_updateRays")
        call ed_updateRays (moveRays) 
        !call Timers_stop   ("ed_updateRays")

        call ed_commIsTransportDone(isTransportDone)
        activeRays = .not.isTransportDone
     end do

     call Timers_stop ("Transport Rays")

     ed_savedDepoStepNumber = ed_currentStepNumber
     ed_depoVarValid = .TRUE.
     ed_depoDt = timeStep

  END IF

  ed_prevDt = timeStep
!
!
!     ...Print out the current energy profile.
!
!
  call Timers_start        ("Print Energy Stamp")
  call ed_printEnergyStamp (timeStep, timeSimulation)
  call Timers_stop         ("Print Energy Stamp")
!
!
!     ...Since some cell variables have changed, update the cell data that depends
!        on them.
!
!
  call Timers_start    ("Update Plasma")
  if (reuseDepo) then
     scaleFact = timeStep / ed_depoDt
     call ed_updatePlasma (blockCount, blockList, scaleFact)
     call Grid_notifySolnDataUpdate()
  else
     call ed_updatePlasma (blockCount, blockList)
  end if
  call Timers_stop     ("Update Plasma")
!
!
!     ...If an IO has been performed previosly on the rays,  don't do this the next time
!        step.
!
!
  if (ed_laserIOWrite) then
      ed_laserIOWrite = .false.
      call IO_endRayWrite ()
  end if
  
  call Timers_stop ("EnergyDeposition")

  return
end subroutine EnergyDeposition
