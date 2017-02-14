!!****if* source/physics/RadTrans/RadTransMain/RadTrans_init
!!
!!  NAME 
!!
!!  RadTrans_init
!!
!!  SYNOPSIS
!!
!!  call RadTrans_init()
!!
!!  DESCRIPTION 
!!    Initialize radiative transfer unit
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"

subroutine RadTrans_init()
  use RadTrans_data
  use rt_interface, ONLY: rt_init
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, Driver_getComm
  implicit none


  call Driver_getMype(MESH_COMM,rt_meshMe)
  call Driver_getMype(GLOBAL_COMM,rt_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, rt_globalNumProcs)

  call RuntimeParameters_get ("useRadTrans", rt_useRadTrans)

  ! Store physical constants:
  call PhysicalConstants_get("speed of light",rt_speedlt)
  rt_speedlt3 = rt_speedlt / 3.0
  call PhysicalConstants_get("Stefan-Boltzmann",rt_radconst)
  rt_radconst = 4.0 * rt_radconst / rt_speedlt
  call PhysicalConstants_get("Boltzmann", rt_boltz)

  call RuntimeParameters_get("rt_dtFactor", rt_dtFactor)
  call RuntimeParameters_get("meshCopyCount", rt_meshCopyCount)  
  call Driver_getMype(MESH_ACROSS_COMM, rt_acrossMe)
  call Driver_getComm(MESH_ACROSS_COMM, rt_acrossComm)
  call Driver_getComm(GLOBAL_COMM,rt_globalComm)

  rt_dbgContext%step = -1
  rt_dbgContext%group = -1

  call rt_init

  if (rt_meshMe == MASTER_PE) print *, 'RadTrans initialized'
  
  return
end subroutine RadTrans_init
