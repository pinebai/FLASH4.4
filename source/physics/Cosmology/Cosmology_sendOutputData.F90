!!****f* source/physics/Cosmology/Cosmology_sendOutputData
!!
!! NAME
!!  Cosmology_sendOutputData
!!
!! SYNOPSIS
!!
!!  Cosmology_sendOutputData()
!!
!! DESCRIPTION
!!
!!  Sends any data from the Cosmology unit (stored in Cosmology_data) to the
!!  IO unit, through the IO_setScalar interface, for checkpointing.
!!  Can retrieve the data upon restart in Cosmology_init through the
!!  IO_getScalar interface.
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
!!  IO_updateScalars, which calls this function to let Cosmology know that
!!  the scalars need to be updated for a checkpoint or plot file
!!
!!
!!
!!***

subroutine Cosmology_sendOutputData()

implicit none

end subroutine Cosmology_sendOutputData
