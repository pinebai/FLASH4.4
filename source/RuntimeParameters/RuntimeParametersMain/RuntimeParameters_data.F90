!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_data
!!
!!  NAME
!!    RuntimeParameters_data
!!
!!  SYNOPSIS
!!
!!    use RuntimeParameters_data
!!    
!!  DESCRIPTION
!!
!!    RuntimeParameters_data is a fortran module 
!!    In addition to holding runtime data this module holds a 
!!    few overloaded functions as well.  In order for overloading
!!    to work correctly and portably in fortran the functions 
!!    need to be inside of a fortran module     
!!
!!
!!
!!  ARGUMENTS
!!    
!!***



module RuntimeParameters_data

! No ONLY: below because the nameValueLL_data also contains a bunch of interfaces.
!  They are included when you use RuntimeParameters_data -- lets the user avoid
!  using both modules.
  use nameValueLL_data 
  
#include "constants.h" 

  
  type (context_type), save :: parameter

  !store the unknown parameters
  !these are the parameters that are found in the flash.par but not declared in any Config
  !store them and then write them to the logfile.  Just a warning to users. Hopefully there
  !are not more than 50 of them...
  character (len=MAX_STRING_LENGTH), save :: rp_ignoredParams(50)
  integer, save :: rp_numIgnoredParams
  integer, save :: rp_globalMe, rp_globalComm
  
end module RuntimeParameters_data
