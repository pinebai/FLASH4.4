!!****if* source/monitors/Logfile/LogfileMain/Logfile_init
!!
!! NAME
!!  Logfile_init
!!
!! SYNOPSIS
!!  Logfile_init()
!!               
!!                
!!
!! DESCRIPTION
!!  Performs Logfile initializations.  Initializes the runtime
!!  runtime parameters needed by the Logfile Unit and then calls
!!  Logfile_create
!!
!!
!! ARGUMENTS
!!
!! PARAMETERS
!! 
!!   These are the runtime parameters used by Logfile unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!    log_file [STRING]
!!        Name of log file to create
!!    run_comment [STRING]
!!        Comment for run
!!    run_number [STRING]
!!        Identification number for run
!!    basenm
!!
!! NOTES
!!  variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!!
!!***

subroutine Logfile_init()

  use Logfile_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
       RuntimeParameters_set
  use Logfile_interface, ONLY : Logfile_create
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm, Driver_getNumProcs

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: posBlank
  character(len=MAX_STRING_LENGTH) :: strLog, test_baseName = "NotSet_",&
       procID
  integer,dimension(70) :: tempProcID

  integer i,j,k
  
  call RuntimeParameters_get("run_comment", log_runComment)
  call RuntimeParameters_get("run_number", log_runNum)

  call RuntimeParameters_get("log_file", log_fileName)
  call Driver_getMype(GLOBAL_COMM,log_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM,log_globalNumProcs)
  call Driver_getComm(GLOBAL_COMM,log_globalComm )

  i=1
  k=0
  j=log_globalMe
  do while (i<=log_globalNumProcs)
     i=i*10
     k=k+1
     tempProcID(k)=mod(j,10)
     j=j/10
  end do
  procID=char(48+tempProcID(k))

  do i = 1,k-1
     procID=procID(:i)//char(48+tempProcID(k-i))
  end do

#ifdef FLASH_IO
  !  Set the default log_file to something similar to basenm if defaults are not used
  call RuntimeParameters_get("basenm",test_baseName)
  if (test_baseName .NE. "NotSet_") then   ! Must account for IO not being included...
     if ((test_baseName .NE. 'flash_') .AND.                    &
          &      (log_fileName .EQ. "flash.log")) then
        posBlank = index(test_baseName,' ')
        if (posBlank > 2) then  !  Rule out " " and "_" and "A" ...
           if (test_baseName(posBlank-1:posBlank-1)=='_') then ! only when name ends in _ ...
              strLog = test_baseName(:posBlank-2) // '.log'  ! Remove trailing _
              log_fileName = strLog
              call RuntimeParameters_set('log_file',log_fileName)
           end if
        end if
     endif
  endif
#endif

  j=index(log_fileName,'.')-1
  log_fileNameLocal=log_fileName(:j)//procID(:k)//'.log'

#ifdef FLASH_GRID_PARAMESH3OR4
  call RuntimeParameters_get("enableMaskedGCFill", log_enableGcMaskLogging)
#endif

  call Logfile_create ()
  

end subroutine Logfile_init
