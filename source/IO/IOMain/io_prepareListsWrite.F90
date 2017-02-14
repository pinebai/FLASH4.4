!!****if* source/IO/IOMain/io_prepareListsWrite
!!
!! NAME
!!
!!  io_prepareListsWrite
!!
!!
!! SYNOPSIS
!!
!!  io_prepareListsWrite()
!!          
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  This function prepares the runtime parameters lists and scalar
!!  lists for output.  This function can be used by all IO routines
!!  which is why it sits at the highest level
!!  
!!
!!  
!!  
!!  
!!  
!!
!! ARGUMENTS
!!  none
!!
!!***


subroutine io_prepareListsWrite() 
  
  use IO_data, ONLY : io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues

  use RuntimeParameters_interface, ONLY : &
    RuntimeParameters_getAllInt, RuntimeParameters_getAllReal, &
    RuntimeParameters_getAllStr, RuntimeParameters_getAllLog
  

  implicit none
  
#include "Flash_mpi.h"
#include "constants.h"

  integer             :: i, maxNumParms
  logical, allocatable, dimension(:) :: changedArray !! Array to store if values changed
  

#ifdef DEBUG_IO
!!     if (io_globalMe==MASTER_PE) then
        print*,'[io_prepareListsWrite]From RP unit, #s params to write are...'
        print 999,io_numRealParms,io_numIntParms,io_numStrParms,io_numLogParms
999     format('                      ... Real:',I3,',Int:',I3,',Str:',I3,',Log:',I3)
!!     end if
#endif

  maxNumParms = max(io_numrealParms,io_numIntParms,io_numStrParms,io_numLogParms)

  !! allocate the space for the Runtime Parameters
  
  !! array to store if they changed or not (we dont need this info here)
  allocate (changedArray(maxNumParms))


  allocate (io_realParmNames (io_numrealParms))
  allocate (io_realParmValues (io_numRealParms))
  
  allocate (io_intParmNames (io_numIntParms))
  allocate (io_intParmValues (io_numIntParms))
    
  allocate (io_strParmNames (io_numStrParms))
  allocate (io_strParmValues (io_numStrParms))
  
  allocate (io_logParmNames (io_numLogParms))
  allocate (io_logParmValues (io_numLogParms))
  allocate (io_logToIntParmValues (io_numLogParms))



  ! get runtime parameters

  call RuntimeParameters_getAllInt(io_numIntParms, io_intParmNames, &
       io_intParmValues, changedArray)
  
  
  call RuntimeParameters_getAllReal(io_numRealParms, &
       io_realParmNames, io_realParmValues, changedArray)
  

  call RuntimeParameters_getAllStr(io_numStrParms, &
       io_strParmNames, io_strParmValues, changedArray)


  call RuntimeParameters_getAllLog(io_numLogParms, &
      io_logParmNames, io_logParmValues, changedArray)
  
  ! translate the logical values to integers for calling c routine
  do i=1, io_numLogParms
     if(io_logParmValues(i)) then
        io_logToIntParmValues(i) = 1
     else
        io_logToIntParmValues(i) = 0
     end if
  end do
  


  call RuntimeParameters_getAllInt(io_numIntParms, &
       io_intParmNames, io_intParmValues, changedArray)
  

  call RuntimeParameters_getAllReal(io_numRealParms, &
       io_realParmNames, io_realParmValues, changedArray)
  

  call RuntimeParameters_getAllStr(io_numStrParms, &
       io_strParmNames, io_strParmValues, changedArray)


  call RuntimeParameters_getAllLog(io_numLogParms, &
      io_logParmNames, io_logParmValues, changedArray)
  
!! deallocate changedArray
   deallocate (changedArray)


  ! translate the logical values to integers for calling c routine
  do i=1, io_numLogParms
     if(io_logParmValues(i)) then
        io_logToIntParmValues(i) = 1
     else
        io_logToIntParmValues(i) = 0
     end if
  end do
  





  !! Now do the scalars!

  
  call io_getNumScalarsReal(io_numRealScalars)
  call io_getNumScalarsInt(io_numIntScalars)
  call io_getNumScalarsStr(io_numStrScalars)
  call io_getNumScalarsLog(io_numLogScalars)
  maxNumParms = max(io_numRealScalars, io_numIntScalars, &
                & io_numStrScalars, io_numLogScalars)
#ifdef DEBUG_IO
!!     if (io_globalMe==MASTER_PE) then
        print*,'[io_prepareListsWrite]From RP unit, #s scalars to write are...'
        print 999,io_numRealScalars,io_numIntScalars,io_numStrScalars,io_numLogScalars

!!     end if
#endif
  allocate(changedArray(maxNumParms))

  if (io_numRealScalars > 0) then
     allocate (io_realScalarNames (io_numRealScalars))
     allocate (io_realScalarValues (io_numRealScalars))
     call io_getAllScalarsReal(io_numRealScalars, &
          io_realScalarNames, io_realScalarValues, changedArray)
  else
     allocate (io_realScalarNames (1))
     allocate (io_realScalarValues (1))
     io_realScalarNames(1) = "dummy"
     io_realScalarValues(1)= -1.
     io_numRealScalars = 1
  end if


  if (io_numIntScalars > 0) then
     allocate (io_intScalarNames (io_numIntScalars))
     allocate (io_intScalarValues (io_numIntScalars))
     call io_getAllScalarsInt(io_numIntScalars, &
          io_intScalarNames, io_intScalarValues, changedArray)
  else
     allocate (io_intScalarNames (1))
     allocate (io_intScalarValues (1))
     io_intScalarNames(1) = "dummy"
     io_intScalarValues(1)= -1
     io_numIntScalars = 1
  end if
  

  if (io_numStrScalars > 0) then
     allocate (io_strScalarNames (io_numStrScalars))
     allocate (io_strScalarValues (io_numStrScalars))
     call io_getAllScalarsStr(io_numStrScalars, &
          io_strScalarNames, io_strScalarValues, changedArray)
  else
     allocate (io_strScalarNames (1))
     allocate (io_strScalarValues (1))
     io_strScalarNames(1) = "dummy"
     io_strScalarValues(1)= "none"
     io_numStrScalars = 1
  end if


  if (io_numLogScalars > 0) then
     allocate (io_logScalarNames (io_numLogScalars))
     allocate (io_logScalarValues (io_numLogScalars))
     allocate (io_logToIntScalarValues (io_numLogScalars))
     call io_getAllScalarsLog(io_numLogScalars, &
          io_logScalarNames, io_logScalarValues, changedArray)
     ! translate the logical values to integers for calling c routine
     do i=1, io_numLogScalars
        if(io_logScalarValues(i)) then
           io_logToIntScalarValues(i) = 1
        else
           io_logToIntScalarValues(i) = 0
        end if
     end do
  else
     allocate (io_logScalarNames (1))
     allocate (io_logScalarValues (1))
     allocate (io_logToIntScalarValues (1))
     io_logScalarNames(1) = "dummy"
     io_logToIntScalarValues(1) = -1
     io_numLogScalars = 1
  end if
    

 deallocate(changedArray)


end subroutine io_prepareListsWrite

