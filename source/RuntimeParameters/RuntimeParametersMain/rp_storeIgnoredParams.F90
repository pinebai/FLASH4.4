!!****if* source/RuntimeParameters/RuntimeParametersMain/rp_storeIgnoredParams
!!
!! NAME
!!  rp_storeIgnoredParams
!!
!! SYNOPSIS
!!
!!  rp_storeIgnoredParams(integer :: numIgnoredParams)
!!
!! DESCRIPTION
!!
!!   This routine stores the ignored runtime parameters.  They will later be 
!!   stamped to the logfile.  (They can not be stamped directly as the logfile
!!   doesn't get initialized unitl after the runtime parameters have been read in.)
!!   
!!   The 'ignored' runtime parameters are those parameters that are found 
!!   in the flash.par but never declared in any Config file.  
!!   (Runtime parameters MUST be declared in a Config file with the 
!!   keyword PARAMETER.  They must also be given a default
!!   value in a Config file.  The flash.par file then allows the default values
!!   to be overwritten)  
!!
!!   Sometimes these ignored parameters are just that ignored.  They don't matter
!!   to the simulation or are left over from some other flash.par file and the
!!   user doesn't care.  However, we do want to warn the user that something in
!!   their flash.par file did not make it into the simulation.
!!
!!
!! ARGUMENTS
!!
!! numIgnoredParams:       the number of parameters found in the flash.par file but
!!                         not found in any Config file for the simulation
!!
!!
!!
!!***

subroutine rp_storeIgnoredParams (name)

  use RuntimeParameters_data, only : rp_ignoredParams, rp_numIgnoredParams

  implicit none

#include "constants.h"

  character(len=MAX_STRING_LENGTH), intent(in) :: name

  if (rp_numIgnoredParams < 50) then
     rp_numIgnoredParams = rp_numIgnoredParams + 1 

     rp_ignoredParams(rp_numIgnoredParams) = trim(name)
  else
     print*,'rp_storeIgnoredParams: no space for "',trim(name),'"!'
  end if

  return

end subroutine rp_storeIgnoredParams


  


