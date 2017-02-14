!!****if* source/IO/IOMain/io_prepareSimInfo
!!
!! NAME
!!
!!  io_prepareSimInfo
!!
!!
!! SYNOPSIS
!!
!!  io_prepareSimInfo()
!!          
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  This function prepares the simulation info for io output.
!!  This includes things like, the setup call, file creation time,
!!  build directory and other statistics about the run
!!  
!!  
!!
!! ARGUMENTS
!!   none   
!!
!!***


subroutine io_prepareSimInfo() 
  
  use IO_data, ONLY : io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp
  

  implicit none
  
#include "Flash_mpi.h"
#include "constants.h"

  character(len=80)   :: setup_flashRelease


  io_flashRelease = setup_flashRelease()


  call setup_buildstats(io_buildDate, io_buildDir, io_buildMachine, io_setupCall, &
             io_cflags, io_fflags)

  call setup_buildstamp (io_setupTimeStamp, io_buildTimeStamp, MAX_STRING_LENGTH)
  
  call current_date_time(io_fileCreationTime)

  call removeNullChar(io_buildDate)
  call removeNullChar(io_buildDir)
  call removeNullChar(io_buildMachine)
  call removeNullChar(io_setupCall)
  call removeNullChar(io_cflags)
  call removeNullChar(io_fflags)
  call removeNullChar(io_setupTimeStamp)
  call removeNullChar(io_buildTimeStamp)
  call removeNullChar(io_fileCreationTime)
        


end subroutine io_prepareSimInfo

