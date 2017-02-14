!!****f* source/Grid/Grid_sendOutputData
!!
!! NAME
!!  Grid_sendOutputData
!!
!! SYNOPSIS
!!
!!  call Grid_sendOutputData()
!!  
!! DESCRIPTION 
!!
!!  This routine allows the Grid unit to checkpoint any scalar data
!!  stored in the Grid_data Fortran modules and is called by the
!!  routine IO_updateScalars before checkpointing.  To send data to
!!  the IO unit this routine calls IO_setScalar.  In addition this
!!  routine may prepare any other data the IO unit needs for
!!  checkpointing which the Grid unit owns.
!!
!!  For example, the Grid unit owns the variable which defines the
!!  grid geometry gr_geometry. This value needs to be checkpointed so
!!  that visualization tools can determine if the simulation ran with
!!  cartesian, spherical, cylindrical etc. coordinates.  To send
!!  scalar data such as gr_geometry to the IO unit to be checkpointed,
!!  the routine calls IO_setScalar("geometry", gr_geometry)
!!   
!!
!!  ARGUMENTS  
!!
!!  SEE ALSO
!!   
!!   IO_setScalar, IO_updateScalars
!!
!!***

subroutine Grid_sendOutputData()

implicit none

end subroutine Grid_sendOutputData
