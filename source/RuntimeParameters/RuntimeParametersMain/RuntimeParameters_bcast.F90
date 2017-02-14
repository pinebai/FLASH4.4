!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_bcast
!!
!! NAME
!!  RuntimeParameters_bcast
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_bcast()
!!
!! DESCRIPTION
!!
!! Broadcasts parameters from the current processor to the other processors.
!! Only the master processor MASTER_PE reads in runtime parameters.
!! 
!! ARGUMENTS
!!
!!        
!!
!!
!!
!!***

subroutine RuntimeParameters_bcast()

  use RuntimeParameters_data, ONLY : rp_globalMe,parameter

  implicit none
  

  call nameValueLL_bcast(parameter, rp_globalMe)

  return

end subroutine RuntimeParameters_bcast


  

