!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commInitComm
!!
!! NAME
!!
!!  ed_commInitComm
!!
!! SYNOPSIS
!!
!!  call ed_commInitComm ()
!!
!! DESCRIPTION
!!
!!  Initializes laser communication.
!!
!! ARGUMENTS
!!
!!  No arguments
!!
!! SIDE EFFECTS
!!
!!  Updates ed_commGlobalRays, ed_commCounterState and ed_commPipelineState.
!!
!! NOTES
!!
!!  This implementation will only work if the sum of the rays over all MPI ranks is
!!  less than the maximum value of a Fortran integer.
!!
!!  Future: The code can be optimized in the following ways:
!!
!!     1). The MPI_AllReduce is only needed once if the number of rays
!!         is constant over a complete run.
!!
!!     2). The pipeline initialization only needs to happen after a
!!         mesh regrid event.
!!
!!***

subroutine ed_commInitComm ()

  use EnergyDeposition_data,   ONLY : ed_maxRayCount, &
                                      ed_meshComm,    &
                                      ed_rayCount

  use ed_commData,             ONLY : ed_commChannelSize,        &
                                      ed_commCounterState,       &
                                      ed_commGlobalRays,         &
                                      ed_commLogUnit,            &
                                      ed_commPipelineState,      &
                                      ed_commRaysBetweenMsgTest, &
                                      ed_commRaysUntilMsgTest,   &
                                      COUNTER_ACTIVE,            &
                                      PIPELINE_ACTIVE

  use Driver_interface,        ONLY : Driver_abortFlash
  use Grid_interface,          ONLY : Grid_getNeighProcList
  use UTCounter_sharedCounter, ONLY : UTCounter_startCounter

  use UTPipeline,              ONLY : UTPipeline_init,     &
                                      UTPipeline_initComm
  implicit none

#include "EnergyDeposition.h"

  include "Flash_mpi.h"

  integer :: error
  integer :: numNeigh

  logical, parameter   :: includeMyProc = .false.
  integer, pointer     :: tmpNeighProcList (:)
  integer, allocatable :: neighProcList (:)
!
!
!   ...Calculate the total number of rays over all MPI ranks.
!
!
  call MPI_AllReduce (ed_rayCount,       &
                      ed_commGlobalRays, &
                      1,                 &
                      FLASH_INTEGER,     &
                      MPI_SUM,           &
                      ed_meshComm,       &
                      error              )

  if (error /= MPI_SUCCESS) then
      call Driver_abortFlash ('[ed_commInitComm] MPI_AllReduce fail')
  end if

  if (ed_commGlobalRays <= 0) then
      call Driver_abortFlash ('[ed_commInitComm] number of globalRays <= 0... overflow?')
  end if
!
!
!   ...Obtain the proc IDs of all surrounding neighbor blocks in order to initialize the
!      communication pipeline.  Note that the Grid subroutine will return a null pointer
!      if all neighboring blocks exist on my MPI rank and includeMyProc = .false.  If all
!      neighboring blocks exist on my MPI rank then we pass a placeholder 1-element array
!      to the pipeline initialization subroutine.
!
!
  call Grid_getNeighProcList (includeMyProc, tmpNeighProcList, numNeigh)

  if (numNeigh > 0) then
      if (.not. associated (tmpNeighProcList)) then
           call Driver_abortFlash ("[ed_commInitComm] Neigh proc list not assoc.?")
     end if

     allocate (neighProcList (1:numNeigh))
     neighProcList (:) = tmpNeighProcList (:)
  else
     allocate (neighProcList(1))
     neighProcList (1) = -1
  end if

  if (associated (tmpNeighProcList)) deallocate (tmpNeighProcList)
  nullify (tmpNeighProcList)

  call UTPipeline_init (RAY_ATTR_COUNT,     &
                        ed_maxRayCount,     &
                        ed_commChannelSize, &
                        ed_meshComm,        &
                        numNeigh,           &
                        neighProcList,      &
                        ed_commLogUnit      )

  deallocate (neighProcList)

  call UTPipeline_initComm    ()
  call UTCounter_startCounter (ed_commGlobalRays)

  ed_commCounterState     = COUNTER_ACTIVE
  ed_commPipelineState    = PIPELINE_ACTIVE
  ed_commRaysUntilMsgTest = ed_commRaysBetweenMsgTest
!
!
!    ...Ready!
!
!  
  return
end subroutine ed_commInitComm
