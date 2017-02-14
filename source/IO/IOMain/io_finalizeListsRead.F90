!!****if* source/IO/IOMain/io_finalizeListsRead
!!
!! NAME
!!
!!  io_finalizeListsRead
!!
!!
!! SYNOPSIS
!!
!!  io_finalizeListsRead() 
!!          
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!   Given the runtime parameter and scalar lists read into arrays from the  
!!   restart checkpoint file, this routine sets the values in the IO data structures
!!   for storing the parameters and the scalars and deallocates the arrays 
!!   that were used for reading them in from the checkpoint.
!!  
!!
!! ARGUMENTS
!! 
!!
!!***


subroutine io_finalizeListsRead()
  
  use IO_data, ONLY : io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValuesPrev
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_setPrev
  
  
  implicit none
  

  integer :: i


  do i=1, io_numRealParmsPrev
     call removeNullChar(io_realParmNamesPrev(i))
     call RuntimeParameters_setPrev(io_realParmNamesPrev(i), io_realParmValuesPrev(i))
  end do

  do i=1, io_numIntParmsPrev
     call removeNullChar(io_intParmNamesPrev(i))
     call RuntimeParameters_setPrev(io_intParmNamesPrev(i), io_intParmValuesPrev(i))
  end do

  do i=1, io_numStrParmsPrev
     call removeNullChar(io_strParmNamesPrev(i))
     call removeNullChar(io_strParmValuesPrev(i))
     call RuntimeParameters_setPrev(io_strParmNamesPrev(i), io_strParmValuesPrev(i))
  end do


  do i=1, io_numLogParmsPrev
     call removeNullChar(io_logParmNamesPrev(i))
     if (io_logToIntParmValuesPrev(i) == 1) then
        io_logParmValuesPrev(i) = .true.
        call RuntimeParameters_setPrev(io_logParmNamesPrev(i), io_logParmValuesPrev(i))
     else if (io_logToIntParmValuesPrev(i) == 0) then
        io_logParmValuesPrev(i) = .false.
        call RuntimeParameters_setPrev(io_logParmNamesPrev(i), io_logParmValuesPrev(i))
     else if (io_logToIntParmValuesPrev(i) == -1) then   
        !!This was a dummy value do nothing
     else
        call Driver_abortFlash("Error reading LogParmValues")
     end if

  end do

  !! deallocate space for RuntimeParms
  if (io_numIntParmsPrev > 0) then
     deallocate(io_intParmNamesPrev)
     deallocate(io_intParmValuesPrev)
  end if
  
  if (io_numRealParmsPrev > 0) then
     deallocate(io_realParmNamesPrev)
     deallocate(io_realParmValuesPrev)
  end if
  
  if (io_numStrParmsPrev > 0) then
     deallocate(io_strParmNamesPrev)
     deallocate(io_strParmValuesPrev)
  end if
  
  if (io_numLogParmsPrev > 0) then
     deallocate(io_logParmNamesPrev)
     deallocate(io_logParmValuesPrev)
     deallocate(io_logToIntParmValuesPrev)
  end if


  do i=1, io_numRealScalars
     call removeNullChar(io_realScalarNames(i))
     call io_setPrevScalarReal(io_realScalarNames(i), io_realScalarValues(i))
  end do

 

  do i=1, io_numIntScalars
     call removeNullChar(io_intScalarNames(i))
     call io_setPrevScalarInt(io_intScalarNames(i), io_intScalarValues(i))
  end do

  do i=1, io_numStrScalars
     call removeNullChar(io_strScalarNames(i))
     call removeNullChar(io_strScalarValues(i))
     call io_setPrevScalarStr(io_strScalarNames(i), io_strScalarValues(i))
  end do


  do i=1, io_numLogScalars
     call removeNullChar(io_logScalarNames(i))
     if (io_logToIntScalarValues(i) == 1) then
        io_logScalarValues(i) = .true.
        call io_setPrevScalarLog(io_logScalarNames(i), io_logScalarValues(i))
     else if(io_logToIntScalarValues(i) == 0) then
        io_logScalarValues(i) = .false.
        call io_setPrevScalarLog(io_logScalarNames(i), io_logScalarValues(i))
     else if(io_logToIntScalarValues(i) == -1) then
        !! Do nothing this is a dummy value
     else
        call Driver_abortFlash("Error reading LogScalarValues")
     end if

  end do


  deallocate(io_intScalarNames)
  deallocate(io_intScalarValues)
  
  deallocate(io_realScalarNames)
  deallocate(io_realScalarValues)
  
  deallocate(io_strScalarNames)
  deallocate(io_strScalarValues)
  
  deallocate(io_logScalarNames)
  deallocate(io_logScalarValues)
  deallocate(io_logToIntScalarValues)
  

end subroutine io_finalizeListsRead
