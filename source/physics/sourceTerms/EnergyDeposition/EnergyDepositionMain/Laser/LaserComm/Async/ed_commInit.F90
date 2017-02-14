!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commInit
!!
!! NAME
!!
!!  ed_commInit
!!
!! SYNOPSIS
!!
!!  call ed_commInit()
!!
!! DESCRIPTION
!!
!!  Initializates the laser runtime parameters.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!***

#include "EnergyDeposition.h"

subroutine ed_commInit()
  use EnergyDeposition_data, ONLY : ed_meshComm, ed_meshMe, &
       ed_laser3Din2D, ed_maxRayCount
  use ed_commData, ONLY : ed_commLog, ed_commDebug, ed_commChannelSize, &
       ed_commRaysBetweenMsgTest, ed_commLogUnit, ed_commRayPosIndex
  !$ use ed_commData, ONLY : ed_commLock
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use UTCounter_sharedCounter, ONLY : UTCounter_init
  implicit none
  integer :: masterRank, test_int

  call RuntimeParameters_get ("ed_commLog", ed_commLog)
  call RuntimeParameters_get ("ed_commDebug", ed_commDebug)
  call RuntimeParameters_get ("ed_commChannelSize", ed_commChannelSize)

  if (ed_commChannelSize > ed_maxRayCount) then
     !Things will work provided ed_maxRayCount is greater than or
     !equal to ed_commChannelSize.  I recommend that ed_maxRayCount
     != 10*channelSize so that there is a high degree of message
     !parallelism
     call Driver_abortFlash('Make ed_commChannelSize <= ed_maxRayCount')
  end if

  call RuntimeParameters_get ("ed_commRaysBetweenMsgTest", &
       ed_commRaysBetweenMsgTest)
  if (ed_commRaysBetweenMsgTest < 1) then
     call Driver_abortFlash("Invalid ed_commRaysBetweenMsgTest")
  end if

  if (ed_commLog) then
     call Logfile_open(ed_commLogUnit,.true.)
  else
     ed_commLogUnit = -1
  end if

  if (ed_meshMe == 0) then
     print *, "WARNING : The global maximum number of rays is", &
          HUGE(test_int)
  end if

  masterRank = 0
  call UTCounter_init(ed_meshComm, masterRank, ed_commLogUnit)
  !This is called before the mesh is allocated

  if (ed_laser3Din2D) then
     ed_commRayPosIndex = (/ RAY_POSX, RAY_POSZ, RAY_POSY /)
  else
     ed_commRayPosIndex = (/ RAY_POSX, RAY_POSY, RAY_POSZ /)
  end if

  !$ call omp_init_lock(ed_commLock)
end subroutine ed_commInit
