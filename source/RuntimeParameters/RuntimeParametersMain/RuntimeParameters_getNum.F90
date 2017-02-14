!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_getNum
!!
!! NAME
!!  RuntimeParameters_getNum
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_getNum(integer(out) :: nparms)
!!                      
!!
!! DESCRIPTION
!!
!!  Returns the number of parameters of a given type.  The
!!  routines below are implemented; they are not overloaded,
!!  so call them directly.
!!
!!  RuntimeParameter_getNumReal
!!  RuntimeParameter_getNumInt
!!  RuntimeParameter_getNumStr
!!  RuntimeParameter_getNumLog
!!
!! ARGUMENTS
!!
!!  nparms:     number of parameters
!!
!!
!!
!!
!!***


   
subroutine RuntimeParameters_getNumReal (nparms)

  use RuntimeParameters_data, ONLY : parameter

  implicit none
  integer, intent(out)                      :: nparms

  call nameValueLL_getNumReal(parameter, nparms)

  return
  
end subroutine RuntimeParameters_getNumReal


   
subroutine RuntimeParameters_getNumInt (nparms)

  use RuntimeParameters_data, ONLY : parameter

  implicit none
  integer, intent(out)                :: nparms

  call nameValueLL_getNumInt(parameter, nparms)

  return
  
end subroutine RuntimeParameters_getNumInt


   
subroutine RuntimeParameters_getNumStr (nparms)

  use RuntimeParameters_data, ONLY : parameter

  implicit none
  integer, intent(out)          :: nparms

  call nameValueLL_getNumStr(parameter, nparms)

  return
  
end subroutine RuntimeParameters_getNumStr


   
subroutine RuntimeParameters_getNumLog (nparms)

  use RuntimeParameters_data, ONLY : parameter

  implicit none
  integer, intent(out)                      :: nparms

  call nameValueLL_getNumLog(parameter, nparms)

  return
  
end subroutine RuntimeParameters_getNumLog



