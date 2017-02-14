!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setPEarrayJmax
!!
!! NAME
!!
!!  op_setPEarrayJmax
!!
!! SYNOPSIS
!!
!!  call op_setPEarrayJmax ()
!!
!! DESCRIPTION
!!
!!  This routine sets the Jmax array for determining the photoelectric
!!  cross sections according to the F.Biggs and R.Lighthill report:
!!
!!       Analytical Approximations for X-Ray Cross Sections II
!!       Frank Biggs and Ruth Lighthill
!!       Weapons Effects Research Department
!!       Sandia Laboratories, December 1971
!!
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEarrayJmax ()

  use op_lowTempData,   ONLY : op_Jmax
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  logical :: arrayExists
  logical :: goodShape
!
!
!   ...Check allocation and dimensions of the j-index delimiter array.
!
!
  arrayExists = allocated (op_Jmax)

  if (.not.arrayExists) then
       call Driver_abortFlash ('[op_setPEarrayJmax] ERROR: array op_Jmax not allocated')
  end if

  goodShape = size (op_Jmax) == 100

  if (.not.goodShape) then
       call Driver_abortFlash ('[op_setPEarrayJmax] ERROR: array op_Jmax has bad shape')
  end if
!
!
!   ...Set the j-index delimiter array
!
!
  op_Jmax (  1) = 8
  op_Jmax (  2) = 8
  op_Jmax (  3) = 7
  op_Jmax (  4) = 7
  op_Jmax (  5) = 7
  op_Jmax (  6) = 7
  op_Jmax (  7) = 7
  op_Jmax (  8) = 6
  op_Jmax (  9) = 6
  op_Jmax ( 10) = 7

  op_Jmax ( 11) = 6
  op_Jmax ( 12) = 6
  op_Jmax ( 13) = 6
  op_Jmax ( 14) = 6
  op_Jmax ( 15) = 6
  op_Jmax ( 16) = 6
  op_Jmax ( 17) = 6
  op_Jmax ( 18) = 6
  op_Jmax ( 19) = 7
  op_Jmax ( 20) = 7

  op_Jmax ( 21) = 6
  op_Jmax ( 22) = 5
  op_Jmax ( 23) = 5
  op_Jmax ( 24) = 5
  op_Jmax ( 25) = 6
  op_Jmax ( 26) = 6
  op_Jmax ( 27) = 6
  op_Jmax ( 28) = 8
  op_Jmax ( 29) = 7
  op_Jmax ( 30) = 8

  op_Jmax ( 31) = 9
  op_Jmax ( 32) = 9
  op_Jmax ( 33) = 9
  op_Jmax ( 34) = 9
  op_Jmax ( 35) = 9
  op_Jmax ( 36) = 9
  op_Jmax ( 37) = 9
  op_Jmax ( 38) = 10
  op_Jmax ( 39) = 10
  op_Jmax ( 40) = 9

  op_Jmax ( 41) = 9
  op_Jmax ( 42) = 9
  op_Jmax ( 43) = 9
  op_Jmax ( 44) = 9
  op_Jmax ( 45) = 9
  op_Jmax ( 46) = 9
  op_Jmax ( 47) = 9
  op_Jmax ( 48) = 8
  op_Jmax ( 49) = 9
  op_Jmax ( 50) = 10

  op_Jmax ( 51) = 10
  op_Jmax ( 52) = 10
  op_Jmax ( 53) = 10
  op_Jmax ( 54) = 11
  op_Jmax ( 55) = 11
  op_Jmax ( 56) = 12
  op_Jmax ( 57) = 12
  op_Jmax ( 58) = 12
  op_Jmax ( 59) = 12
  op_Jmax ( 60) = 11

  op_Jmax ( 61) = 13
  op_Jmax ( 62) = 13
  op_Jmax ( 63) = 13
  op_Jmax ( 64) = 12
  op_Jmax ( 65) = 12
  op_Jmax ( 66) = 12
  op_Jmax ( 67) = 12
  op_Jmax ( 68) = 12
  op_Jmax ( 69) = 12
  op_Jmax ( 70) = 12

  op_Jmax ( 71) = 12
  op_Jmax ( 72) = 12
  op_Jmax ( 73) = 12
  op_Jmax ( 74) = 12
  op_Jmax ( 75) = 12
  op_Jmax ( 76) = 12
  op_Jmax ( 77) = 12
  op_Jmax ( 78) = 12
  op_Jmax ( 79) = 12
  op_Jmax ( 80) = 12

  op_Jmax ( 81) = 12
  op_Jmax ( 82) = 12
  op_Jmax ( 83) = 12
  op_Jmax ( 84) = 12
  op_Jmax ( 85) = 12
  op_Jmax ( 86) = 12
  op_Jmax ( 87) = 12
  op_Jmax ( 88) = 12
  op_Jmax ( 89) = 12
  op_Jmax ( 90) = 12

  op_Jmax ( 91) = 12
  op_Jmax ( 92) = 12
  op_Jmax ( 93) = 12
  op_Jmax ( 94) = 12
  op_Jmax ( 95) = 12
  op_Jmax ( 96) = 12
  op_Jmax ( 97) = 12
  op_Jmax ( 98) = 12
  op_Jmax ( 99) = 12
  op_Jmax (100) = 12
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setPEarrayJmax
