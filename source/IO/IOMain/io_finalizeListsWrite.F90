!!****if* source/IO/IOMain/io_finalizeListsWrite
!!
!! NAME
!!
!!  io_finalizeListsWrite
!!
!!
!! SYNOPSIS
!!
!!  io_finalizeListsWrite() 
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
!! 
!!
!!***


subroutine io_finalizeListsWrite()
  
  use IO_data, ONLY : io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues
  
  
  implicit none
  


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
  

end subroutine io_finalizeListsWrite
