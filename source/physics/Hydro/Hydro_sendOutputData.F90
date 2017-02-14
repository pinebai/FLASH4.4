!!****f* source/physics/Hydro/Hydro_sendOutputData
!!
!! NAME
!!  Hydro_sendOutputData
!!
!! SYNOPSIS
!!  
!!  Hydro_sendOutputData()
!!  
!! DESCRIPTION 
!!   
!!  Sends any data from the Hydro unit (stored in Hydro_data) to the
!!  IO unit, through the IO_setScalar interface, for checkpointing.
!!  Can retrieve the data upon restart in Hydro_init through the 
!!  IO_getScalar interface.
!!  
!!  This is the API stub implementation. 
!!
!! ARGUMENTS 
!!
!! USES
!!  
!!  IO_setScalar, which puts a scalar into the IO unit so it can make 
!!  it into a checkpoint. 
!!
!! USED BY
!!
!!  IO_updateScalars, which calls this function to let Hydro know 
!!  it's time to checkpoint.
!!
!! EXAMPLE
!!
!!  Here's what an implementation of this function might do:
!!
!!  USE IO_interface, ONLY : IO_setScalar
!!  call IO_setScalar("time", simTime)
!!
!! 
!!***

subroutine Hydro_sendOutputData()

implicit none
end subroutine Hydro_sendOutputData

