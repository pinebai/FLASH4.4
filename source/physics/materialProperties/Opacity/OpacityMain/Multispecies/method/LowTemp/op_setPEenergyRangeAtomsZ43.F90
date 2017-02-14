!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setPEenergyRangeAtomsZ43
!!
!! NAME
!!
!!  op_setPEenergyRangeAtomsZ43
!!
!! SYNOPSIS
!!
!!  call op_setPEenergyRangeAtomsZ43 ()
!!
!! DESCRIPTION
!!
!!  This routine sets the energy range (boundaries) for determining the photoelectric
!!  cross sections according to the F.Biggs and R.Lighthill report:
!!
!!       Analytical Approximations for X-Ray Cross Sections II
!!       Frank Biggs and Ruth Lighthill
!!       Weapons Effects Research Department
!!       Sandia Laboratories, December 1971
!!
!!  The atomic elements set in this routine are in the range: 1 =< Z =< 43.
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEenergyRangeAtomsZ43 ()

  use op_lowTempData,   ONLY : op_PEenergyRange
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

# include "Opacity.h"

  logical :: arrayExists
  logical :: goodShape
!
!
!   ...Check allocation and dimensions of the energy range array.
!
!
  arrayExists = allocated (op_PEenergyRange)

  if (.not.arrayExists) then
       call Driver_abortFlash ('[op_setPEenergyRangeAtomsZ43] ERROR: array op_PEenergyRange not allocated')
  end if

  goodShape =       (size (op_PEenergyRange,1) == HIGH-LOW+1)   &
              .and. (size (op_PEenergyRange,2) == 13)           &
              .and. (size (op_PEenergyRange,3) == 100)

  if (.not.goodShape) then
       call Driver_abortFlash ('[op_setPEenergyRangeAtomsZ43] ERROR: array op_PEenergyRange has bad shape')
  end if
!
!
!   ...Set the energy range array. Last index denotes the
!      atomic order number Z.
!
!
  op_PEenergyRange (LOW:HIGH, 1,  1) = (/    0.0100,    0.0140/)
  op_PEenergyRange (LOW:HIGH, 2,  1) = (/    0.0140,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 3,  1) = (/    0.1000,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 4,  1) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  1) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  1) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  1) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 8,  1) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  2) = (/    0.0100,    0.0250/)
  op_PEenergyRange (LOW:HIGH, 2,  2) = (/    0.0250,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 3,  2) = (/    0.1000,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 4,  2) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  2) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  2) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  2) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 8,  2) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  3) = (/    0.0100,    0.0550/)
  op_PEenergyRange (LOW:HIGH, 2,  3) = (/    0.0550,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3,  3) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  3) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  3) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  3) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  3) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  4) = (/    0.0100,    0.1110/)
  op_PEenergyRange (LOW:HIGH, 2,  4) = (/    0.1110,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3,  4) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  4) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  4) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  4) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  4) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  5) = (/    0.0100,    0.1880/)
  op_PEenergyRange (LOW:HIGH, 2,  5) = (/    0.1880,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3,  5) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  5) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  5) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  5) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  5) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  6) = (/    0.0100,    0.2840/)
  op_PEenergyRange (LOW:HIGH, 2,  6) = (/    0.2840,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3,  6) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  6) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  6) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  6) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  6) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  7) = (/    0.0100,    0.4000/)
  op_PEenergyRange (LOW:HIGH, 2,  7) = (/    0.4000,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3,  7) = (/    0.8000,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  7) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  7) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  7) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7,  7) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  8) = (/    0.0100,    0.5330/)
  op_PEenergyRange (LOW:HIGH, 2,  8) = (/    0.5330,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 3,  8) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  8) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  8) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  8) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,  9) = (/    0.0100,    0.6870/)
  op_PEenergyRange (LOW:HIGH, 2,  9) = (/    0.6870,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 3,  9) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4,  9) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5,  9) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6,  9) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 10) = (/    0.0100,    0.0180/)
  op_PEenergyRange (LOW:HIGH, 2, 10) = (/    0.0180,    0.8670/)
  op_PEenergyRange (LOW:HIGH, 3, 10) = (/    0.8670,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 10) = (/    4.0000,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 10) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 10) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 10) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 11) = (/    0.0100,    0.0320/)
  op_PEenergyRange (LOW:HIGH, 2, 11) = (/    0.0320,    1.0730/)
  op_PEenergyRange (LOW:HIGH, 3, 11) = (/    1.0730,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 11) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 11) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 11) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 12) = (/    0.0100,    0.0500/)
  op_PEenergyRange (LOW:HIGH, 2, 12) = (/    0.0500,    1.3050/)
  op_PEenergyRange (LOW:HIGH, 3, 12) = (/    1.3050,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 12) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 12) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 12) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 13) = (/    0.0100,    0.0730/)
  op_PEenergyRange (LOW:HIGH, 2, 13) = (/    0.0730,    1.5600/)
  op_PEenergyRange (LOW:HIGH, 3, 13) = (/    1.5600,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 13) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 13) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 13) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 14) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 14) = (/    0.1000,    1.8390/)
  op_PEenergyRange (LOW:HIGH, 3, 14) = (/    1.8390,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 14) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 14) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 14) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 15) = (/    0.0100,    0.1300/)
  op_PEenergyRange (LOW:HIGH, 2, 15) = (/    0.1300,    2.1440/)
  op_PEenergyRange (LOW:HIGH, 3, 15) = (/    2.1440,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 15) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 15) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 15) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 16) = (/    0.0100,    0.1650/)
  op_PEenergyRange (LOW:HIGH, 2, 16) = (/    0.1650,    2.4720/)
  op_PEenergyRange (LOW:HIGH, 3, 16) = (/    2.4720,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 16) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 16) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 16) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 17) = (/    0.0100,    0.2020/)
  op_PEenergyRange (LOW:HIGH, 2, 17) = (/    0.2020,    2.8240/)
  op_PEenergyRange (LOW:HIGH, 3, 17) = (/    2.8240,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 17) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 17) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 17) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 18) = (/    0.0100,    0.2450/)
  op_PEenergyRange (LOW:HIGH, 2, 18) = (/    0.2450,    3.2030/)
  op_PEenergyRange (LOW:HIGH, 3, 18) = (/    3.2030,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 18) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 18) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 18) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 19) = (/    0.0100,    0.0180/)
  op_PEenergyRange (LOW:HIGH, 2, 19) = (/    0.0180,    0.2940/)
  op_PEenergyRange (LOW:HIGH, 3, 19) = (/    0.2940,    3.6070/)
  op_PEenergyRange (LOW:HIGH, 4, 19) = (/    3.6070,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 19) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 19) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 19) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 20) = (/    0.0100,    0.0250/)
  op_PEenergyRange (LOW:HIGH, 2, 20) = (/    0.0250,    0.3460/)
  op_PEenergyRange (LOW:HIGH, 3, 20) = (/    0.3460,    4.0370/)
  op_PEenergyRange (LOW:HIGH, 4, 20) = (/    4.0370,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 20) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 20) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 20) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 21) = (/    0.0100,    0.0320/)
  op_PEenergyRange (LOW:HIGH, 2, 21) = (/    0.0320,    0.4010/)
  op_PEenergyRange (LOW:HIGH, 3, 21) = (/    0.4010,    4.4910/)
  op_PEenergyRange (LOW:HIGH, 4, 21) = (/    4.4910,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 21) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 21) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 22) = (/    0.0100,    0.4560/)
  op_PEenergyRange (LOW:HIGH, 2, 22) = (/    0.4560,    4.9660/)
  op_PEenergyRange (LOW:HIGH, 3, 22) = (/    4.9660,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 22) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 22) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 23) = (/    0.0100,    0.5130/)
  op_PEenergyRange (LOW:HIGH, 2, 23) = (/    0.5130,    5.4650/)
  op_PEenergyRange (LOW:HIGH, 3, 23) = (/    5.4650,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 23) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 23) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 24) = (/    0.0100,    0.5750/)
  op_PEenergyRange (LOW:HIGH, 2, 24) = (/    0.5750,    5.9890/)
  op_PEenergyRange (LOW:HIGH, 3, 24) = (/    5.9890,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 4, 24) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 24) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 25) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 25) = (/    0.1000,    0.6400/)
  op_PEenergyRange (LOW:HIGH, 3, 25) = (/    0.6400,    6.5390/)
  op_PEenergyRange (LOW:HIGH, 4, 25) = (/    6.5390,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 25) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 25) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 26) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 26) = (/    0.1000,    0.7080/)
  op_PEenergyRange (LOW:HIGH, 3, 26) = (/    0.7080,    7.1120/)
  op_PEenergyRange (LOW:HIGH, 4, 26) = (/    7.1120,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 26) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 26) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 27) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 27) = (/    0.1000,    0.7790/)
  op_PEenergyRange (LOW:HIGH, 3, 27) = (/    0.7790,    7.7090/)
  op_PEenergyRange (LOW:HIGH, 4, 27) = (/    7.7090,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 5, 27) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 27) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 28) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 28) = (/    0.1000,    0.8540/)
  op_PEenergyRange (LOW:HIGH, 3, 28) = (/    0.8540,    0.8710/)
  op_PEenergyRange (LOW:HIGH, 4, 28) = (/    0.8710,    1.0080/)
  op_PEenergyRange (LOW:HIGH, 5, 28) = (/    1.0080,    8.3320/)
  op_PEenergyRange (LOW:HIGH, 6, 28) = (/    8.3320,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 28) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 28) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 29) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 29) = (/    0.1000,    0.9330/)
  op_PEenergyRange (LOW:HIGH, 3, 29) = (/    0.9330,    1.0960/)
  op_PEenergyRange (LOW:HIGH, 4, 29) = (/    1.0960,    8.9810/)
  op_PEenergyRange (LOW:HIGH, 5, 29) = (/    8.9810,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 6, 29) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 29) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 30) = (/    0.0100,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 2, 30) = (/    0.1000,    1.0200/)
  op_PEenergyRange (LOW:HIGH, 3, 30) = (/    1.0200,    1.0430/)
  op_PEenergyRange (LOW:HIGH, 4, 30) = (/    1.0430,    1.1930/)
  op_PEenergyRange (LOW:HIGH, 5, 30) = (/    1.1930,    9.6590/)
  op_PEenergyRange (LOW:HIGH, 6, 30) = (/    9.6590,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 30) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 30) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 31) = (/    0.0100,    0.0170/)
  op_PEenergyRange (LOW:HIGH, 2, 31) = (/    0.0170,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 3, 31) = (/    0.1000,    1.1150/)
  op_PEenergyRange (LOW:HIGH, 4, 31) = (/    1.1150,    1.1420/)
  op_PEenergyRange (LOW:HIGH, 5, 31) = (/    1.1420,    1.3000/)
  op_PEenergyRange (LOW:HIGH, 6, 31) = (/    1.3000,   10.3670/)
  op_PEenergyRange (LOW:HIGH, 7, 31) = (/   10.3670,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 31) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 31) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 32) = (/    0.0100,    0.0280/)
  op_PEenergyRange (LOW:HIGH, 2, 32) = (/    0.0280,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 3, 32) = (/    0.1000,    1.2170/)
  op_PEenergyRange (LOW:HIGH, 4, 32) = (/    1.2170,    1.2480/)
  op_PEenergyRange (LOW:HIGH, 5, 32) = (/    1.2480,    1.4130/)
  op_PEenergyRange (LOW:HIGH, 6, 32) = (/    1.4130,   11.1040/)
  op_PEenergyRange (LOW:HIGH, 7, 32) = (/   11.1040,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 32) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 32) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 33) = (/    0.0100,    0.0410/)
  op_PEenergyRange (LOW:HIGH, 2, 33) = (/    0.0410,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 3, 33) = (/    0.1000,    1.3230/)
  op_PEenergyRange (LOW:HIGH, 4, 33) = (/    1.3230,    1.3590/)
  op_PEenergyRange (LOW:HIGH, 5, 33) = (/    1.3590,    1.5300/)
  op_PEenergyRange (LOW:HIGH, 6, 33) = (/    1.5300,   11.8670/)
  op_PEenergyRange (LOW:HIGH, 7, 33) = (/   11.8670,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 33) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 33) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 34) = (/    0.0100,    0.0550/)
  op_PEenergyRange (LOW:HIGH, 2, 34) = (/    0.0550,    0.1000/)
  op_PEenergyRange (LOW:HIGH, 3, 34) = (/    0.1000,    1.4340/)
  op_PEenergyRange (LOW:HIGH, 4, 34) = (/    1.4340,    1.4750/)
  op_PEenergyRange (LOW:HIGH, 5, 34) = (/    1.4750,    1.6520/)
  op_PEenergyRange (LOW:HIGH, 6, 34) = (/    1.6520,   12.6580/)
  op_PEenergyRange (LOW:HIGH, 7, 34) = (/   12.6580,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 34) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 34) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 35) = (/    0.0100,    0.0710/)
  op_PEenergyRange (LOW:HIGH, 2, 35) = (/    0.0710,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 35) = (/    0.8000,    1.5510/)
  op_PEenergyRange (LOW:HIGH, 4, 35) = (/    1.5510,    1.5970/)
  op_PEenergyRange (LOW:HIGH, 5, 35) = (/    1.5970,    1.7820/)
  op_PEenergyRange (LOW:HIGH, 6, 35) = (/    1.7820,   13.4740/)
  op_PEenergyRange (LOW:HIGH, 7, 35) = (/   13.4740,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 35) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 35) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 36) = (/    0.0100,    0.0900/)
  op_PEenergyRange (LOW:HIGH, 2, 36) = (/    0.0900,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 36) = (/    0.8000,    1.6750/)
  op_PEenergyRange (LOW:HIGH, 4, 36) = (/    1.6750,    1.7270/)
  op_PEenergyRange (LOW:HIGH, 5, 36) = (/    1.7270,    1.9210/)
  op_PEenergyRange (LOW:HIGH, 6, 36) = (/    1.9210,   14.3230/)
  op_PEenergyRange (LOW:HIGH, 7, 36) = (/   14.3230,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 36) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 36) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 37) = (/    0.0100,    0.0310/)
  op_PEenergyRange (LOW:HIGH, 2, 37) = (/    0.0310,    0.1100/)
  op_PEenergyRange (LOW:HIGH, 3, 37) = (/    0.1100,    1.8050/)
  op_PEenergyRange (LOW:HIGH, 4, 37) = (/    1.8050,    1.8630/)
  op_PEenergyRange (LOW:HIGH, 5, 37) = (/    1.8630,    2.0650/)
  op_PEenergyRange (LOW:HIGH, 6, 37) = (/    2.0650,   15.2000/)
  op_PEenergyRange (LOW:HIGH, 7, 37) = (/   15.2000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 37) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 37) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 38) = (/    0.0100,    0.0200/)
  op_PEenergyRange (LOW:HIGH, 2, 38) = (/    0.0200,    0.1330/)
  op_PEenergyRange (LOW:HIGH, 3, 38) = (/    0.1330,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 4, 38) = (/    0.8000,    1.9400/)
  op_PEenergyRange (LOW:HIGH, 5, 38) = (/    1.9400,    2.0070/)
  op_PEenergyRange (LOW:HIGH, 6, 38) = (/    2.0070,    2.2160/)
  op_PEenergyRange (LOW:HIGH, 7, 38) = (/    2.2160,   16.1050/)
  op_PEenergyRange (LOW:HIGH, 8, 38) = (/   16.1050,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 38) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,10, 38) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 39) = (/    0.0100,    0.0250/)
  op_PEenergyRange (LOW:HIGH, 2, 39) = (/    0.0250,    0.1560/)
  op_PEenergyRange (LOW:HIGH, 3, 39) = (/    0.1560,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 4, 39) = (/    0.8000,    2.0790/)
  op_PEenergyRange (LOW:HIGH, 5, 39) = (/    2.0790,    2.1560/)
  op_PEenergyRange (LOW:HIGH, 6, 39) = (/    2.1560,    2.3730/)
  op_PEenergyRange (LOW:HIGH, 7, 39) = (/    2.3730,   17.0380/)
  op_PEenergyRange (LOW:HIGH, 8, 39) = (/   17.0380,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 39) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,10, 39) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 40) = (/    0.0100,    0.1800/)
  op_PEenergyRange (LOW:HIGH, 2, 40) = (/    0.1800,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 40) = (/    0.8000,    2.2230/)
  op_PEenergyRange (LOW:HIGH, 4, 40) = (/    2.2230,    2.3070/)
  op_PEenergyRange (LOW:HIGH, 5, 40) = (/    2.3070,    2.5330/)
  op_PEenergyRange (LOW:HIGH, 6, 40) = (/    2.5330,   17.9980/)
  op_PEenergyRange (LOW:HIGH, 7, 40) = (/   17.9980,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 40) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 40) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 41) = (/    0.0100,    0.2040/)
  op_PEenergyRange (LOW:HIGH, 2, 41) = (/    0.2040,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 41) = (/    0.8000,    2.3700/)
  op_PEenergyRange (LOW:HIGH, 4, 41) = (/    2.3700,    2.4640/)
  op_PEenergyRange (LOW:HIGH, 5, 41) = (/    2.4640,    2.6980/)
  op_PEenergyRange (LOW:HIGH, 6, 41) = (/    2.6980,   18.9860/)
  op_PEenergyRange (LOW:HIGH, 7, 41) = (/   18.9860,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 41) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 41) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 42) = (/    0.0100,    0.2280/)
  op_PEenergyRange (LOW:HIGH, 2, 42) = (/    0.2280,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 42) = (/    0.8000,    2.5210/)
  op_PEenergyRange (LOW:HIGH, 4, 42) = (/    2.5210,    2.6250/)
  op_PEenergyRange (LOW:HIGH, 5, 42) = (/    2.6250,    2.8670/)
  op_PEenergyRange (LOW:HIGH, 6, 42) = (/    2.8670,   20.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 42) = (/   20.0000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 42) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 42) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 43) = (/    0.0100,    0.2530/)
  op_PEenergyRange (LOW:HIGH, 2, 43) = (/    0.2530,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 43) = (/    0.8000,    2.6770/)
  op_PEenergyRange (LOW:HIGH, 4, 43) = (/    2.6770,    2.7930/)
  op_PEenergyRange (LOW:HIGH, 5, 43) = (/    2.7930,    3.0430/)
  op_PEenergyRange (LOW:HIGH, 6, 43) = (/    3.0430,   21.0440/)
  op_PEenergyRange (LOW:HIGH, 7, 43) = (/   21.0440,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 43) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 43) = (/  500.0000,10000.0000/)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setPEenergyRangeAtomsZ43
