!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_getNumIgn
!!
!! NAME
!!  RuntimeParameters_getNumIgn
!!
!! SYNOPSIS
!!
!!  call RuntimeParameters_getNumIgn(integer(out) :: numIgnoredParams)
!!
!! DESCRIPTION
!!
!!  Returns the number of ignored runtime parameters.
!!
!!  The 'ignored' runtime parameters are those that are found 
!!  in the flash.par but are no recognized because they were
!!  not declared in any Config file included in the setup.  
!!  (Runtime parameters MUST be declared in a Config file with the 
!!  keyword PARAMETER.  They must also be given a default
!!  value in a Config file.  The flash.par file then allows the
!!  default values to be overwritten.)  
!!
!! ARGUMENTS
!!
!!  numIgnoredParams: number of times that invalid PARAMETER lines
!!                    were found while processing flash.par.
!!
!!***

subroutine RuntimeParameters_getNumIgn (numIgnoredParams)

  use RuntimeParameters_data, ONLY : rp_numIgnoredParams

  implicit none
  integer, intent(out)                      :: numIgnoredParams

  numIgnoredParams = rp_numIgnoredParams

  return
  
end subroutine RuntimeParameters_getNumIgn
