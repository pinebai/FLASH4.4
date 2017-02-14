!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_init
!!
!! NAME
!!  RuntimeParameters_init
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_init(logical(out) :: restart)
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the data need in the Runtime Parameters
!!  Unit.  Reads the parameter file, usually flash.par
!!  and broadcasts parameters to the other processors.
!!  
!!  
!!
!! ARGUMENTS
!!
!! restart - true if run is restarted from checkpoint, false if starting
!!           from scratch
!!           
!!
!!
!!***


#define ENABLE_CHK_FILE_OPTION



subroutine RuntimeParameters_init( restart)

  use RuntimeParameters_data
  use IO_interface, ONLY : IO_initRPsFromCheckpoint
  use RuntimeParameters_interface, ONLY : RuntimeParameters_read, &
    RuntimeParameters_get, RuntimeParameters_set, RuntimeParameters_bcast
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm

  implicit none  

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"


  logical, intent(out) :: restart

  integer :: checkpointFlagExists, parfileFlagExists
  logical :: useParfile
  character(len=MAX_STRING_LENGTH) :: checkpointFile
  character(len=MAX_STRING_LENGTH) :: parmFile
  character(len=20) :: pCheck
  character(len=1000) :: errorMessage = " "
  integer :: flags(2), ierr
  
  rp_numIgnoredParams = 0
  call Driver_getComm(GLOBAL_COMM, rp_globalComm)
  call Driver_getMype(GLOBAL_COMM, rp_globalMe)

  if (rp_globalMe == MASTER_PE) then
     
     ! Find out if there's a -c checkpointFile on the command
     ! line, if so, we're going to try to get the runtime parameters 
     ! from the checkoint file.  Otherwise, get them 
     ! from the flash.par file.
 
     !! just until we get checkpoint implemented
     checkpointFlagExists = 0

#ifdef ENABLE_CHK_FILE_OPTION
     !! DEV: Implemention in process: -chk_file checkpointFile on the command line
     !! set checkpoint_flag to true and return cpname
     pCheck='-chk_file'
     call rp_getOpt(pCheck, 1, checkpointFlagExists, checkpointFile)
#endif

     
     ! Get the parameter file name from the command line, or use
     !'flash.par' if no name is specified.  
     
     pCheck='-par_file'
     call rp_getOpt(pCheck, 1, parfileFlagExists, parmFile)

     if(parfileFlagExists /= 1) parmFile = 'flash.par'

     flags(1) = parfileFlagExists
     flags(2) = checkpointFlagExists

     useParfile = ((checkpointFlagExists==0) .OR. (parfileFlagExists/=0))
  end if

  call MPI_Bcast(flags, 2, FLASH_INTEGER, MASTER_PE, rp_globalComm, ierr)
  if (rp_globalMe .NE. MASTER_PE) then
     checkpointFlagExists = flags(2)
     useParfile = ((flags(2)==0) .OR. (flags(1)/=0))
  end if

  call nameValueLL_initContext(parameter)

  if (checkpointFlagExists == 1) then
     restart = .true.
     call IO_initRPsFromCheckpoint( checkpointFile, ierr)
  end if


  if (rp_globalMe == MASTER_PE) then

     ! Set default values from Config files for all
     ! runtime parameters parsed by the setup script
     ! in a Config file included in the simulation.
     !
     ! If checkpointFlagExists==1, then this will only add
     ! RPs that for some reasons are defined in Config files
     ! but were not in the checkpoint file; it will not modify
     ! values of RPs that were already added above; so RPs from
     ! the checkpoint get precedence. 
     call rp_initParameters(parmFile)

     if (checkpointFlagExists == 1) then
        restart = .true.
        call RuntimeParameters_set('restart', restart)
     else
        restart = .true.
     end if

     if (useParfile) then

        ! Get the parameter file name from the command line, or use
        !'flash.par' if no name is specified.  

        if(parfileFlagExists /= 1) parmFile = 'flash.par'

        ! Read the input parameter file and set the values of all
        ! runtime parameters.
        call RuntimeParameters_read(parmFile)

     endif
     
     ! Broadcast parameters to the other processors.
     call RuntimeParameters_get('restart', restart)
     call RuntimeParameters_bcast ()
     
  else

     ! Non-master processors receive parameter information from the
     ! master.
     call RuntimeParameters_bcast ()
     call RuntimeParameters_get('restart', restart)

  endif

end subroutine RuntimeParameters_init
