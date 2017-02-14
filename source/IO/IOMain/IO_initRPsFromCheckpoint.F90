!!****if* source/IO/IOMain/IO_initRPsFromCheckpoint
!!
!! NAME
!!
!!  IO_initRPsFromCheckpoint
!!
!! SYNOPSIS
!!
!!  call IO_initRPsFromCheckpoint(character(len=*)(IN) :: filename,
!!                                integer(OUT)         :: ierr)
!!
!! DESCRIPTION
!!
!!   Reads runtime parameter information and some integer scalars
!!   from a checkpoint file, then sets the runtime parameters in
!!   the linked list for the current run accordingly.
!!
!!   
!! ARGUMENTS
!!
!!
!!   filename : name of checkpoint file
!!
!!   ierr : if not 0, indicates that something went wrong
!!
!! NOTES
!!
!!   This subroutine is called early during initialization.
!!   The implementation, and subroutines called from it, must
!!   operate correctly before IO_init has been called.
!!
!!   This is not meant to be called by users.
!!   
!!***

subroutine IO_initRPsFromCheckpoint( filename, ierr)
  use RuntimeParameters_interface, ONLY : RuntimeParameters_set, RuntimeParameters_add
  use Driver_interface, ONLY : Driver_abortFlash
  use IO_data, ONLY : io_baseName, io_checkpointFileNumber, &
       io_maxParms, &
       io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_outputSplitNum, io_comm, io_chkptFileID, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels,&
       io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, &
       io_logToIntParmValuesPrev, io_globalMe

  implicit none
#include "constants.h"
  character(len=*),intent(IN) :: filename
  integer,intent(OUT) :: ierr
  integer :: i, nRunNum, ioStatus, tempInt
  integer :: checkpointFileNumber
  logical :: success
  character(len=MAX_STRING_LENGTH)     :: runNumStr


  io_numRealParms = io_maxParms
  io_numIntParms = io_maxParms
  io_numStrParms = io_maxParms
  io_numLogParms = io_maxParms

  !! allocate the space for the Runtime Parameters

  ! as in io_prepareListsWrite
  !! array to store if they changed or not (we dont need this info here)
  allocate (io_realParmNames (io_numrealParms))
  allocate (io_realParmValues (io_numRealParms))
  
  allocate (io_intParmNames (io_numIntParms))
  allocate (io_intParmValues (io_numIntParms))
    
  allocate (io_strParmNames (io_numStrParms))
  allocate (io_strParmValues (io_numStrParms))
  
  allocate (io_logParmNames (io_numLogParms))
  allocate (io_logParmValues (io_numLogParms))
  allocate (io_logToIntParmValues (io_numLogParms))


!! allocate the space for the scalars
  io_numRealScalars = io_maxParms
  io_numIntScalars = io_maxParms
  io_numStrScalars = io_maxParms
  io_numLogScalars = io_maxParms
  allocate (io_realScalarNames (io_numRealScalars))
  allocate (io_realScalarValues (io_numRealScalars))
  
  allocate (io_intScalarNames (io_numIntScalars))
  allocate (io_intScalarValues (io_numIntScalars))
  
  allocate (io_strScalarNames (io_numStrScalars))
  allocate (io_strScalarValues (io_numStrScalars))
  
  allocate (io_logScalarNames (io_numLogScalars))
  allocate (io_logScalarValues (io_numLogScalars))
  allocate (io_logToIntScalarValues (io_numLogScalars))

  call io_readRPsFromCheckpoint(io_globalMe, filename, success, checkpointFileNumber)
  
  if (success .and. (io_globalMe==MASTER_PE)) then

     ! as in io_finalizeListsRead
     do i=1, io_numRealParms
        call removeNullChar(io_realParmNames(i))
        call RuntimeParameters_add(io_realParmNames(i), io_realParmValues(i))
     end do

     do i=1, io_numIntParms
        call removeNullChar(io_intParmNames(i))
        call RuntimeParameters_add(io_intParmNames(i), io_intParmValues(i))
     end do

     do i=1, io_numStrParms
        call removeNullChar(io_strParmNames(i))
        call removeNullChar(io_strParmValues(i))
        call RuntimeParameters_add(io_strParmNames(i), io_strParmValues(i))
        if (io_strParmNames(i)(1:len_trim(io_strParmNames(i))) == 'run_number') then 
           runNumStr = io_strParmValues(i)
#ifdef DEBUG_ALL
           print*,'runNumStr=="',runNumStr,'"...'
#endif
           read (runNumStr, *, iostat=ioStatus) tempInt
           if (ioStatus == 0) then
              tempInt = tempInt + 1
              write (runNumStr, *, iostat=ioStatus) tempInt
              if (ioStatus == 0) then
                 runNumStr = trim(adjustl(runNumStr))
                 call RuntimeParameters_set(io_strParmNames(i), runNumStr)
              end if
           end if
        end if
     end do


     do i=1, io_numLogParms
        call removeNullChar(io_logParmNames(i))
        if (io_logToIntParmValues(i) == 1) then
           io_logParmValues(i) = .true.
           call RuntimeParameters_add(io_logParmNames(i), io_logParmValues(i))
        else if (io_logToIntParmValues(i) == 0) then
           io_logParmValues(i) = .false.
           call RuntimeParameters_add(io_logParmNames(i), io_logParmValues(i))
        else if (io_logToIntParmValues(i) == -1) then   
           !!This was a dummy value do nothing
        else
           call Driver_abortFlash("Error reading LogParmValues")
        end if

     end do

     call RuntimeParameters_set('checkpointFileNumber', checkpointFileNumber)
!!     io_checkpointFileNumber = checkpointFileNumber

     do i=1, io_numIntScalars
        call removeNullChar(io_intScalarNames(i))
        select case (io_intScalarNames(i))
        case('plotfilenumber','forcedplotfilenumber','particlefilenumber')
           call RuntimeParameters_set(io_intScalarNames(i), io_intScalarValues(i))
        case('checkpointfilenumber')
           call RuntimeParameters_set(io_intScalarNames(i), io_intScalarValues(i))
!!           io_checkpointFileNumber = checkpointFileNumber
        end select
     end do

  end if

  ! as in io_finalizeListsWrite
  !! deallocate space for RuntimeParms

  deallocate(io_intParmNames)
  deallocate(io_intParmValues)
  
  deallocate(io_realParmNames)
  deallocate(io_realParmValues)
  
  deallocate(io_strParmNames)
  deallocate(io_strParmValues)

  deallocate(io_logParmNames)
  deallocate(io_logParmValues)
  deallocate(io_logToIntParmValues)


  deallocate(io_intScalarNames)
  deallocate(io_intScalarValues)
  
  deallocate(io_realScalarNames)
  deallocate(io_realScalarValues)
  
  deallocate(io_strScalarNames)
  deallocate(io_strScalarValues)
  
  deallocate(io_logScalarNames)
  deallocate(io_logScalarValues)
  deallocate(io_logToIntScalarValues)


  if (success) then
     ierr = 0
  else
     ierr = 1
     if (io_globalMe==MASTER_PE) print*,'WARNING: reading runtime parameters from file "'//trim(filename)//'" failed.'
  end if

end subroutine IO_initRPsFromCheckpoint
