!!****f* source/RuntimeParameters/RuntimeParameters_getNum
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
!!  RuntimeParameter_getNumReal(integer(out) :: nparms)
!!  RuntimeParameter_getNumInt(integer(out) :: nparms)
!!  RuntimeParameter_getNumStr(integer(out) :: nparms)
!!  RuntimeParameter_getNumLog(integer(out) :: nparms)
!!
!! ARGUMENTS
!!
!!  nparms:     number of parameters
!!
!!
!!
!!***


   
subroutine RuntimeParameters_getNumReal (nparms)
  implicit none
  integer, intent(out)                :: nparms
end subroutine RuntimeParameters_getNumReal


   
subroutine RuntimeParameters_getNumInt (nparms)
  implicit none
  integer, intent(out)               :: nparms
end subroutine RuntimeParameters_getNumInt


   
subroutine RuntimeParameters_getNumStr (nparms)
  implicit none
  integer, intent(out)               :: nparms
end subroutine RuntimeParameters_getNumStr


   
subroutine RuntimeParameters_getNumLog (nparms)
  implicit none
  integer, intent(out)               :: nparms
end subroutine RuntimeParameters_getNumLog



