!!****if* source/IO/IOMain/io_prepareListsRead
!!
!! NAME
!!
!!  io_prepareListsRead
!!
!!
!! SYNOPSIS
!!
!!  io_prepareListsRead()
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


subroutine io_prepareListsRead() 
  
  use IO_data, ONLY : io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValuesPrev
  
  
  implicit none
  

  !! allocate the space for the Runtime Parameters
  allocate (io_realParmNamesPrev (io_numRealParmsPrev))
  allocate (io_realParmValuesPrev  (io_numRealParmsPrev))
  
  allocate (io_intParmNamesPrev  (io_numIntParmsPrev))

  allocate (io_intParmValuesPrev  (io_numIntParmsPrev))
    
  allocate (io_strParmNamesPrev  (io_numStrParmsPrev))
  allocate (io_strParmValuesPrev  (io_numStrParmsPrev))
  
  allocate (io_logParmNamesPrev  (io_numLogParmsPrev))
  allocate (io_logParmValuesPrev  (io_numLogParmsPrev))
  allocate (io_logToIntParmValuesPrev  (io_numLogParmsPrev))


!! allocate the space for the scalars
  allocate (io_realScalarNames (io_numRealScalars))
  allocate (io_realScalarValues (io_numRealScalars))
  
  allocate (io_intScalarNames (io_numIntScalars))
  allocate (io_intScalarValues (io_numIntScalars))
  
  allocate (io_strScalarNames (io_numStrScalars))
  allocate (io_strScalarValues (io_numStrScalars))
  
  allocate (io_logScalarNames (io_numLogScalars))
  allocate (io_logScalarValues (io_numLogScalars))
  allocate (io_logToIntScalarValues (io_numLogScalars))
  




end subroutine io_prepareListsRead
