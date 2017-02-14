!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_packBCs
!!
!! NAME
!!  gr_packBCs
!!
!! SYNOPSIS
!!
!!  gr_packBCs(integer, intent(IN) :: bcILeft, &
!!             integer, intent(IN) :: bcIRight, &
!!             integer, intent(IN) :: bcJLeft, &
!!             integer, intent(IN) :: bcJRight, &
!!             integer, intent(IN) :: bcKLeft, &
!!             integer, intent(IN) :: bcKRight)
!!
!!  
!! DESCRIPTION 
!!
!!  This subroutine is used by gr_createDomain to intialise the Paramesh 
!!  array named boundary_index.  This subroutine is called once 
!!  for every single element of boundary_index.
!!
!!
!! ARGUMENTS   
!!
!!   The input arguments describe the boundary for each block
!!   face by using particular constants.
!!   Constants for boundary conditions are defined in
!!   "constants.h". Any boundary conditions in the allowed
!!   range of numerical values -50..-20 can be used if they
!!   are recognized by the Grid_applyBCEge implementation.
!!   In particular, PERIODIC is not allowed here.
!!
!!   bcILeft : LOW IAXIS block face boundary.
!!   bcIRight : HIGH IAXIS block face boundary.
!!   bcJLeft : LOW JAXIS block face boundary.
!!   bcJRight : HIGH JAXIS block face boundary.
!!   bcKLeft : LOW KAXIS block face boundary.
!!   bcKRight : HIGH KAXIS block face boundary.
!!
!!***

#include "Flash.h"
#include "constants.h"

integer function gr_packBCs(bcILeft, bcIRight, bcJLeft, bcJRight, bcKLeft, bcKRight)
  implicit none
  integer,intent(IN) :: bcILeft, bcIRight, bcJLeft, bcJRight, bcKLeft, bcKRight
  logical :: areTheSame
  integer :: xl,xr,yl,yr,zl,zr,out

  areTheSame = .TRUE.
  if (bcIRight .NE. bcILeft) areTheSame = .FALSE.
#if NDIM > 1
  if (areTheSame .AND. bcJLeft .NE. bcILeft) areTheSame = .FALSE.
  if (areTheSame .AND. bcJRight .NE. bcILeft) areTheSame = .FALSE.
#endif
#if NDIM > 2
  if (areTheSame .AND. bcKLeft .NE. bcILeft) areTheSame = .FALSE.
  if (areTheSame .AND. bcKRight .NE. bcILeft) areTheSame = .FALSE.
#endif
  if (areTheSame) then
     gr_packBCs = bcILeft
  else

     xl = abs(bcILeft) - 20
     xr = abs(bcIRight) - 20
     yl = abs(bcJLeft) - 20
     yr = abs(bcJRight) - 20
     zl = abs(bcKLeft) - 20
     zr = abs(bcKRight) - 20
     out = xr*32 + xl
#if NDIM > 1
     out = (yr*32 + yl) * 1024 + out
#endif
#if NDIM > 2
     out = (zr*32 + zl) * 1024*1024 + out
#endif
     gr_packBCs = -1000 - out
  end if
end function gr_packBCs

integer function gr_extractBCForDirection(packedBCs,axis,leftOrRight)
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer,intent(IN) :: packedBCs,axis,leftOrRight
  integer :: p,b
  if (packedBCs > -1000) then
     gr_extractBCForDirection = packedBCs
  else
     p = abs(packedBCs + 1000)
     if (axis == IAXIS) then
        b = mod(p,1024)
     else if (axis == JAXIS) then
        b = mod(p / 1024, 1024)
     else if (axis == KAXIS) then
        b = p / 1024 / 1024
     else
        call Driver_abortFlash('invalid axis in gr_extractBCForDirection!')
     end if

     if (leftOrRight == LOW) then
        gr_extractBCForDirection = -(mod(b,32) + 20) 
     else
        gr_extractBCForDirection = -(b / 32 + 20)
     end if
  end if
end function gr_extractBCForDirection
