!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_init
!!
!! NAME
!!  
!!  Pipeline_init
!!
!! SYNOPSIS
!! 
!!  call Pipeline_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the Pipeline unit. Anything that is common to all pipelines
!!  during an application should go in here.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine Pipeline_init ()

  use Pipeline_data

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface,            ONLY : Driver_abortFlash

#include "Flash.h"
#include "constants.h"

  implicit none
!
!
!     ...Get the needed external data.
!
!
  call Driver_getComm        (GLOBAL_COMM,  pl_commGlobal)
  call Driver_getComm        (  MESH_COMM,  pl_commMesh  )

  call RuntimeParameters_get ("basenm",     pl_baseName  )
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_init
