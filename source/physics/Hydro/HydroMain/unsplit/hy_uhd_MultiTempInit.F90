!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_MultiTempInit
!!
!! NAME
!!
!!  hy_uhd_MultiTempInit
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_MultiTempInit()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize some unit scope variables which are typically taken
!!  directly from the runtime parameters and which are specific
!!  to the multiTemp variant of the unsplit Hydro implementation.
!!  This must be called once by Hydro_init.F90.  Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This is a stub, to be replaced with an actual iomplementation (which
!!  may live in a subdirectory).
!!***

subroutine hy_uhd_MultiTempInit()
  implicit none
end subroutine hy_uhd_MultiTempInit
