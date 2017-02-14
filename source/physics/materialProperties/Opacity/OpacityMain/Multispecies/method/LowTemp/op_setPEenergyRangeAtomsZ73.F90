!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setPEenergyRangeAtomsZ73
!!
!! NAME
!!
!!  op_setPEenergyRangeAtomsZ73
!!
!! SYNOPSIS
!!
!!  call op_setPEenergyRangeAtomsZ73 ()
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
!!  The atomic elements set in this routine are in the range: 44 =< Z =< 73.
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEenergyRangeAtomsZ73 ()

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
       call Driver_abortFlash ('[op_setPEenergyRangeAtomsZ73] ERROR: array op_PEenergyRange not allocated')
  end if

  goodShape =       (size (op_PEenergyRange,1) == HIGH-LOW+1)   &
              .and. (size (op_PEenergyRange,2) == 13)           &
              .and. (size (op_PEenergyRange,3) == 100)

  if (.not.goodShape) then
       call Driver_abortFlash ('[op_setPEenergyRangeAtomsZ73] ERROR: array op_PEenergyRange has bad shape')
  end if
!
!
!   ...Set the energy range array. Last index denotes the
!      atomic order number Z.
!
!
  op_PEenergyRange (LOW:HIGH, 1, 44) = (/    0.0100,    0.2800/)
  op_PEenergyRange (LOW:HIGH, 2, 44) = (/    0.2800,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 44) = (/    0.8000,    2.8380/)
  op_PEenergyRange (LOW:HIGH, 4, 44) = (/    2.8380,    2.9670/)
  op_PEenergyRange (LOW:HIGH, 5, 44) = (/    2.9670,    3.2240/)
  op_PEenergyRange (LOW:HIGH, 6, 44) = (/    3.2240,   22.1170/)
  op_PEenergyRange (LOW:HIGH, 7, 44) = (/   22.1170,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 44) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 44) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 45) = (/    0.0100,    0.3070/)
  op_PEenergyRange (LOW:HIGH, 2, 45) = (/    0.3070,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 45) = (/    0.8000,    3.0040/)
  op_PEenergyRange (LOW:HIGH, 4, 45) = (/    3.0040,    3.1460/)
  op_PEenergyRange (LOW:HIGH, 5, 45) = (/    3.1460,    3.4120/)
  op_PEenergyRange (LOW:HIGH, 6, 45) = (/    3.4120,   23.2200/)
  op_PEenergyRange (LOW:HIGH, 7, 45) = (/   23.2200,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 45) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 45) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 46) = (/    0.0100,    0.3350/)
  op_PEenergyRange (LOW:HIGH, 2, 46) = (/    0.3350,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 46) = (/    0.8000,    3.1740/)
  op_PEenergyRange (LOW:HIGH, 4, 46) = (/    3.1740,    3.3300/)
  op_PEenergyRange (LOW:HIGH, 5, 46) = (/    3.3300,    3.6050/)
  op_PEenergyRange (LOW:HIGH, 6, 46) = (/    3.6050,   24.3500/)
  op_PEenergyRange (LOW:HIGH, 7, 46) = (/   24.3500,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 46) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 46) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 47) = (/    0.0100,    0.3670/)
  op_PEenergyRange (LOW:HIGH, 2, 47) = (/    0.3670,    0.8000/)
  op_PEenergyRange (LOW:HIGH, 3, 47) = (/    0.8000,    3.3510/)
  op_PEenergyRange (LOW:HIGH, 4, 47) = (/    3.3510,    3.5240/)
  op_PEenergyRange (LOW:HIGH, 5, 47) = (/    3.5240,    3.8060/)
  op_PEenergyRange (LOW:HIGH, 6, 47) = (/    3.8060,   25.5140/)
  op_PEenergyRange (LOW:HIGH, 7, 47) = (/   25.5140,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 47) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 47) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 48) = (/    0.0100,    0.4040/)
  op_PEenergyRange (LOW:HIGH, 2, 48) = (/    0.4040,    3.5370/)
  op_PEenergyRange (LOW:HIGH, 3, 48) = (/    3.5370,    3.7270/)
  op_PEenergyRange (LOW:HIGH, 4, 48) = (/    3.7270,    4.0180/)
  op_PEenergyRange (LOW:HIGH, 5, 48) = (/    4.0180,   26.7110/)
  op_PEenergyRange (LOW:HIGH, 6, 48) = (/   26.7110,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 48) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 48) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 49) = (/    0.0100,    0.0160/)
  op_PEenergyRange (LOW:HIGH, 2, 49) = (/    0.0160,    0.4430/)
  op_PEenergyRange (LOW:HIGH, 3, 49) = (/    0.4430,    3.7300/)
  op_PEenergyRange (LOW:HIGH, 4, 49) = (/    3.7300,    3.9380/)
  op_PEenergyRange (LOW:HIGH, 5, 49) = (/    3.9380,    4.2380/)
  op_PEenergyRange (LOW:HIGH, 6, 49) = (/    4.2380,   27.9400/)
  op_PEenergyRange (LOW:HIGH, 7, 49) = (/   27.9400,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 8, 49) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 49) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 50) = (/    0.0100,    0.0240/)
  op_PEenergyRange (LOW:HIGH, 2, 50) = (/    0.0240,    0.4850/)
  op_PEenergyRange (LOW:HIGH, 3, 50) = (/    0.4850,    0.7140/)
  op_PEenergyRange (LOW:HIGH, 4, 50) = (/    0.7140,    3.9290/)
  op_PEenergyRange (LOW:HIGH, 5, 50) = (/    3.9290,    4.1560/)
  op_PEenergyRange (LOW:HIGH, 6, 50) = (/    4.1560,    4.4650/)
  op_PEenergyRange (LOW:HIGH, 7, 50) = (/    4.4650,   29.2000/)
  op_PEenergyRange (LOW:HIGH, 8, 50) = (/   29.2000,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 50) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,10, 50) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 51) = (/    0.0100,    0.0310/)
  op_PEenergyRange (LOW:HIGH, 2, 51) = (/    0.0310,    0.5280/)
  op_PEenergyRange (LOW:HIGH, 3, 51) = (/    0.5280,    0.7660/)
  op_PEenergyRange (LOW:HIGH, 4, 51) = (/    0.7660,    4.1320/)
  op_PEenergyRange (LOW:HIGH, 5, 51) = (/    4.1320,    4.3810/)
  op_PEenergyRange (LOW:HIGH, 6, 51) = (/    4.3810,    4.6980/)
  op_PEenergyRange (LOW:HIGH, 7, 51) = (/    4.6980,   30.4910/)
  op_PEenergyRange (LOW:HIGH, 8, 51) = (/   30.4910,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 51) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,10, 51) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 52) = (/    0.0100,    0.0400/)
  op_PEenergyRange (LOW:HIGH, 2, 52) = (/    0.0400,    0.5720/)
  op_PEenergyRange (LOW:HIGH, 3, 52) = (/    0.5720,    1.0060/)
  op_PEenergyRange (LOW:HIGH, 4, 52) = (/    1.0060,    4.3410/)
  op_PEenergyRange (LOW:HIGH, 5, 52) = (/    4.3410,    4.6120/)
  op_PEenergyRange (LOW:HIGH, 6, 52) = (/    4.6120,    4.9390/)
  op_PEenergyRange (LOW:HIGH, 7, 52) = (/    4.9390,   31.8140/)
  op_PEenergyRange (LOW:HIGH, 8, 52) = (/   31.8140,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 52) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,10, 52) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 53) = (/    0.0100,    0.0500/)
  op_PEenergyRange (LOW:HIGH, 2, 53) = (/    0.0500,    0.6190/)
  op_PEenergyRange (LOW:HIGH, 3, 53) = (/    0.6190,    1.0720/)
  op_PEenergyRange (LOW:HIGH, 4, 53) = (/    1.0720,    4.5570/)
  op_PEenergyRange (LOW:HIGH, 5, 53) = (/    4.5570,    4.8520/)
  op_PEenergyRange (LOW:HIGH, 6, 53) = (/    4.8520,    5.1880/)
  op_PEenergyRange (LOW:HIGH, 7, 53) = (/    5.1880,   33.1700/)
  op_PEenergyRange (LOW:HIGH, 8, 53) = (/   33.1700,  100.0000/)
  op_PEenergyRange (LOW:HIGH, 9, 53) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,10, 53) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 54) = (/    0.0100,    0.0630/)
  op_PEenergyRange (LOW:HIGH, 2, 54) = (/    0.0630,    0.6720/)
  op_PEenergyRange (LOW:HIGH, 3, 54) = (/    0.6720,    0.9360/)
  op_PEenergyRange (LOW:HIGH, 4, 54) = (/    0.9360,    1.1430/)
  op_PEenergyRange (LOW:HIGH, 5, 54) = (/    1.1430,    4.7820/)
  op_PEenergyRange (LOW:HIGH, 6, 54) = (/    4.7820,    5.1020/)
  op_PEenergyRange (LOW:HIGH, 7, 54) = (/    5.1020,    5.4450/)
  op_PEenergyRange (LOW:HIGH, 8, 54) = (/    5.4450,   34.5610/)
  op_PEenergyRange (LOW:HIGH, 9, 54) = (/   34.5610,  100.0000/)
  op_PEenergyRange (LOW:HIGH,10, 54) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,11, 54) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 55) = (/    0.0100,    0.0790/)
  op_PEenergyRange (LOW:HIGH, 2, 55) = (/    0.0790,    0.7260/)
  op_PEenergyRange (LOW:HIGH, 3, 55) = (/    0.7260,    1.0650/)
  op_PEenergyRange (LOW:HIGH, 4, 55) = (/    1.0650,    1.2170/)
  op_PEenergyRange (LOW:HIGH, 5, 55) = (/    1.2170,    5.0120/)
  op_PEenergyRange (LOW:HIGH, 6, 55) = (/    5.0120,    5.3600/)
  op_PEenergyRange (LOW:HIGH, 7, 55) = (/    5.3600,    5.7130/)
  op_PEenergyRange (LOW:HIGH, 8, 55) = (/    5.7130,   35.9850/)
  op_PEenergyRange (LOW:HIGH, 9, 55) = (/   35.9850,  100.0000/)
  op_PEenergyRange (LOW:HIGH,10, 55) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,11, 55) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 56) = (/    0.0100,    0.0920/)
  op_PEenergyRange (LOW:HIGH, 2, 56) = (/    0.0920,    0.7800/)
  op_PEenergyRange (LOW:HIGH, 3, 56) = (/    0.7800,    1.0610/)
  op_PEenergyRange (LOW:HIGH, 4, 56) = (/    1.0610,    1.1350/)
  op_PEenergyRange (LOW:HIGH, 5, 56) = (/    1.1350,    1.2910/)
  op_PEenergyRange (LOW:HIGH, 6, 56) = (/    1.2910,    5.2470/)
  op_PEenergyRange (LOW:HIGH, 7, 56) = (/    5.2470,    5.6230/)
  op_PEenergyRange (LOW:HIGH, 8, 56) = (/    5.6230,    5.9870/)
  op_PEenergyRange (LOW:HIGH, 9, 56) = (/    5.9870,   37.4410/)
  op_PEenergyRange (LOW:HIGH,10, 56) = (/   37.4410,  100.0000/)
  op_PEenergyRange (LOW:HIGH,11, 56) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 56) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 57) = (/    0.0100,    0.0990/)
  op_PEenergyRange (LOW:HIGH, 2, 57) = (/    0.0990,    0.8320/)
  op_PEenergyRange (LOW:HIGH, 3, 57) = (/    0.8320,    1.1240/)
  op_PEenergyRange (LOW:HIGH, 4, 57) = (/    1.1240,    1.2040/)
  op_PEenergyRange (LOW:HIGH, 5, 57) = (/    1.2040,    1.3630/)
  op_PEenergyRange (LOW:HIGH, 6, 57) = (/    1.3630,    5.4840/)
  op_PEenergyRange (LOW:HIGH, 7, 57) = (/    5.4840,    5.8910/)
  op_PEenergyRange (LOW:HIGH, 8, 57) = (/    5.8910,    6.2660/)
  op_PEenergyRange (LOW:HIGH, 9, 57) = (/    6.2660,   38.9250/)
  op_PEenergyRange (LOW:HIGH,10, 57) = (/   38.9250,  100.0000/)
  op_PEenergyRange (LOW:HIGH,11, 57) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 57) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 58) = (/    0.0100,    0.1050/)
  op_PEenergyRange (LOW:HIGH, 2, 58) = (/    0.1050,    0.8830/)
  op_PEenergyRange (LOW:HIGH, 3, 58) = (/    0.8830,    1.1850/)
  op_PEenergyRange (LOW:HIGH, 4, 58) = (/    1.1850,    1.2730/)
  op_PEenergyRange (LOW:HIGH, 5, 58) = (/    1.2730,    1.4350/)
  op_PEenergyRange (LOW:HIGH, 6, 58) = (/    1.4350,    5.7230/)
  op_PEenergyRange (LOW:HIGH, 7, 58) = (/    5.7230,    6.1640/)
  op_PEenergyRange (LOW:HIGH, 8, 58) = (/    6.1640,    6.5490/)
  op_PEenergyRange (LOW:HIGH, 9, 58) = (/    6.5490,   40.4430/)
  op_PEenergyRange (LOW:HIGH,10, 58) = (/   40.4430,  100.0000/)
  op_PEenergyRange (LOW:HIGH,11, 58) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 58) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 59) = (/    0.0100,    0.1120/)
  op_PEenergyRange (LOW:HIGH, 2, 59) = (/    0.1120,    0.9310/)
  op_PEenergyRange (LOW:HIGH, 3, 59) = (/    0.9310,    1.2420/)
  op_PEenergyRange (LOW:HIGH, 4, 59) = (/    1.2420,    1.3370/)
  op_PEenergyRange (LOW:HIGH, 5, 59) = (/    1.3370,    1.5050/)
  op_PEenergyRange (LOW:HIGH, 6, 59) = (/    1.5050,    5.9640/)
  op_PEenergyRange (LOW:HIGH, 7, 59) = (/    5.9640,    6.4400/)
  op_PEenergyRange (LOW:HIGH, 8, 59) = (/    6.4400,    6.8350/)
  op_PEenergyRange (LOW:HIGH, 9, 59) = (/    6.8350,   41.9910/)
  op_PEenergyRange (LOW:HIGH,10, 59) = (/   41.9910,  100.0000/)
  op_PEenergyRange (LOW:HIGH,11, 59) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 59) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 60) = (/    0.0100,    0.1170/)
  op_PEenergyRange (LOW:HIGH, 2, 60) = (/    0.1170,    0.9780/)
  op_PEenergyRange (LOW:HIGH, 3, 60) = (/    0.9780,    1.2980/)
  op_PEenergyRange (LOW:HIGH, 4, 60) = (/    1.2980,    1.5750/)
  op_PEenergyRange (LOW:HIGH, 5, 60) = (/    1.5750,    6.2080/)
  op_PEenergyRange (LOW:HIGH, 6, 60) = (/    6.2080,    6.7220/)
  op_PEenergyRange (LOW:HIGH, 7, 60) = (/    6.7220,    7.1280/)
  op_PEenergyRange (LOW:HIGH, 8, 60) = (/    7.1280,   43.5690/)
  op_PEenergyRange (LOW:HIGH, 9, 60) = (/   43.5690,  100.0000/)
  op_PEenergyRange (LOW:HIGH,10, 60) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,11, 60) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 61) = (/    0.0100,    0.1220/)
  op_PEenergyRange (LOW:HIGH, 2, 61) = (/    0.1220,    1.0270/)
  op_PEenergyRange (LOW:HIGH, 3, 61) = (/    1.0270,    1.0520/)
  op_PEenergyRange (LOW:HIGH, 4, 61) = (/    1.0520,    1.3570/)
  op_PEenergyRange (LOW:HIGH, 5, 61) = (/    1.3570,    1.4710/)
  op_PEenergyRange (LOW:HIGH, 6, 61) = (/    1.4710,    1.6480/)
  op_PEenergyRange (LOW:HIGH, 7, 61) = (/    1.6480,    6.4590/)
  op_PEenergyRange (LOW:HIGH, 8, 61) = (/    6.4590,    7.0130/)
  op_PEenergyRange (LOW:HIGH, 9, 61) = (/    7.0130,    7.4280/)
  op_PEenergyRange (LOW:HIGH,10, 61) = (/    7.4280,   45.1840/)
  op_PEenergyRange (LOW:HIGH,11, 61) = (/   45.1840,  100.0000/)
  op_PEenergyRange (LOW:HIGH,12, 61) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,13, 61) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 62) = (/    0.0100,    0.1270/)
  op_PEenergyRange (LOW:HIGH, 2, 62) = (/    0.1270,    1.0780/)
  op_PEenergyRange (LOW:HIGH, 3, 62) = (/    1.0780,    1.1060/)
  op_PEenergyRange (LOW:HIGH, 4, 62) = (/    1.1060,    1.4190/)
  op_PEenergyRange (LOW:HIGH, 5, 62) = (/    1.4190,    1.5410/)
  op_PEenergyRange (LOW:HIGH, 6, 62) = (/    1.5410,    1.7230/)
  op_PEenergyRange (LOW:HIGH, 7, 62) = (/    1.7230,    6.7160/)
  op_PEenergyRange (LOW:HIGH, 8, 62) = (/    6.7160,    7.3120/)
  op_PEenergyRange (LOW:HIGH, 9, 62) = (/    7.3120,    7.7360/)
  op_PEenergyRange (LOW:HIGH,10, 62) = (/    7.7360,   46.8340/)
  op_PEenergyRange (LOW:HIGH,11, 62) = (/   46.8340,  100.0000/)
  op_PEenergyRange (LOW:HIGH,12, 62) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,13, 62) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 63) = (/    0.0100,    0.1340/)
  op_PEenergyRange (LOW:HIGH, 2, 63) = (/    0.1340,    1.1310/)
  op_PEenergyRange (LOW:HIGH, 3, 63) = (/    1.1310,    1.1610/)
  op_PEenergyRange (LOW:HIGH, 4, 63) = (/    1.1610,    1.4810/)
  op_PEenergyRange (LOW:HIGH, 5, 63) = (/    1.4810,    1.6140/)
  op_PEenergyRange (LOW:HIGH, 6, 63) = (/    1.6140,    1.8000/)
  op_PEenergyRange (LOW:HIGH, 7, 63) = (/    1.8000,    6.9770/)
  op_PEenergyRange (LOW:HIGH, 8, 63) = (/    6.9770,    7.6180/)
  op_PEenergyRange (LOW:HIGH, 9, 63) = (/    7.6180,    8.0520/)
  op_PEenergyRange (LOW:HIGH,10, 63) = (/    8.0520,   48.5190/)
  op_PEenergyRange (LOW:HIGH,11, 63) = (/   48.5190,  100.0000/)
  op_PEenergyRange (LOW:HIGH,12, 63) = (/  100.0000,  500.0000/)
  op_PEenergyRange (LOW:HIGH,13, 63) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 64) = (/    0.0100,    0.1400/)
  op_PEenergyRange (LOW:HIGH, 2, 64) = (/    0.1400,    1.1850/)
  op_PEenergyRange (LOW:HIGH, 3, 64) = (/    1.1850,    1.2170/)
  op_PEenergyRange (LOW:HIGH, 4, 64) = (/    1.2170,    1.5440/)
  op_PEenergyRange (LOW:HIGH, 5, 64) = (/    1.5440,    1.6880/)
  op_PEenergyRange (LOW:HIGH, 6, 64) = (/    1.6880,    1.8810/)
  op_PEenergyRange (LOW:HIGH, 7, 64) = (/    1.8810,    7.2430/)
  op_PEenergyRange (LOW:HIGH, 8, 64) = (/    7.2430,    7.9300/)
  op_PEenergyRange (LOW:HIGH, 9, 64) = (/    7.9300,    8.3750/)
  op_PEenergyRange (LOW:HIGH,10, 64) = (/    8.3750,   50.2390/)
  op_PEenergyRange (LOW:HIGH,11, 64) = (/   50.2390,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 64) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 65) = (/    0.0100,    0.1470/)
  op_PEenergyRange (LOW:HIGH, 2, 65) = (/    0.1470,    1.2400/)
  op_PEenergyRange (LOW:HIGH, 3, 65) = (/    1.2400,    1.2740/)
  op_PEenergyRange (LOW:HIGH, 4, 65) = (/    1.2740,    1.6100/)
  op_PEenergyRange (LOW:HIGH, 5, 65) = (/    1.6100,    1.7650/)
  op_PEenergyRange (LOW:HIGH, 6, 65) = (/    1.7650,    1.9630/)
  op_PEenergyRange (LOW:HIGH, 7, 65) = (/    1.9630,    7.5140/)
  op_PEenergyRange (LOW:HIGH, 8, 65) = (/    7.5140,    8.2520/)
  op_PEenergyRange (LOW:HIGH, 9, 65) = (/    8.2520,    8.7080/)
  op_PEenergyRange (LOW:HIGH,10, 65) = (/    8.7080,   51.9960/)
  op_PEenergyRange (LOW:HIGH,11, 65) = (/   51.9960,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 65) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 66) = (/    0.0100,    0.1540/)
  op_PEenergyRange (LOW:HIGH, 2, 66) = (/    0.1540,    1.2950/)
  op_PEenergyRange (LOW:HIGH, 3, 66) = (/    1.2950,    1.3320/)
  op_PEenergyRange (LOW:HIGH, 4, 66) = (/    1.3320,    1.6760/)
  op_PEenergyRange (LOW:HIGH, 5, 66) = (/    1.6760,    1.8420/)
  op_PEenergyRange (LOW:HIGH, 6, 66) = (/    1.8420,    2.0460/)
  op_PEenergyRange (LOW:HIGH, 7, 66) = (/    2.0460,    7.7900/)
  op_PEenergyRange (LOW:HIGH, 8, 66) = (/    7.7900,    8.5800/)
  op_PEenergyRange (LOW:HIGH, 9, 66) = (/    8.5800,    9.0460/)
  op_PEenergyRange (LOW:HIGH,10, 66) = (/    9.0460,   53.7880/)
  op_PEenergyRange (LOW:HIGH,11, 66) = (/   53.7880,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 66) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 67) = (/    0.0100,    0.1610/)
  op_PEenergyRange (LOW:HIGH, 2, 67) = (/    0.1610,    1.3510/)
  op_PEenergyRange (LOW:HIGH, 3, 67) = (/    1.3510,    1.3920/)
  op_PEenergyRange (LOW:HIGH, 4, 67) = (/    1.3920,    1.7430/)
  op_PEenergyRange (LOW:HIGH, 5, 67) = (/    1.7430,    1.9230/)
  op_PEenergyRange (LOW:HIGH, 6, 67) = (/    1.9230,    2.1300/)
  op_PEenergyRange (LOW:HIGH, 7, 67) = (/    2.1300,    8.0720/)
  op_PEenergyRange (LOW:HIGH, 8, 67) = (/    8.0720,    8.9180/)
  op_PEenergyRange (LOW:HIGH, 9, 67) = (/    8.9180,    9.3940/)
  op_PEenergyRange (LOW:HIGH,10, 67) = (/    9.3940,   55.6180/)
  op_PEenergyRange (LOW:HIGH,11, 67) = (/   55.6180,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 67) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 68) = (/    0.0100,    0.1690/)
  op_PEenergyRange (LOW:HIGH, 2, 68) = (/    0.1690,    1.4090/)
  op_PEenergyRange (LOW:HIGH, 3, 68) = (/    1.4090,    1.4530/)
  op_PEenergyRange (LOW:HIGH, 4, 68) = (/    1.4530,    1.8120/)
  op_PEenergyRange (LOW:HIGH, 5, 68) = (/    1.8120,    2.0060/)
  op_PEenergyRange (LOW:HIGH, 6, 68) = (/    2.0060,    2.2170/)
  op_PEenergyRange (LOW:HIGH, 7, 68) = (/    2.2170,    8.3580/)
  op_PEenergyRange (LOW:HIGH, 8, 68) = (/    8.3580,    9.2640/)
  op_PEenergyRange (LOW:HIGH, 9, 68) = (/    9.2640,    9.7520/)
  op_PEenergyRange (LOW:HIGH,10, 68) = (/    9.7520,   57.4860/)
  op_PEenergyRange (LOW:HIGH,11, 68) = (/   57.4860,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 68) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 69) = (/    0.0100,    0.1770/)
  op_PEenergyRange (LOW:HIGH, 2, 69) = (/    0.1770,    1.4680/)
  op_PEenergyRange (LOW:HIGH, 3, 69) = (/    1.4680,    1.5150/)
  op_PEenergyRange (LOW:HIGH, 4, 69) = (/    1.5150,    1.8810/)
  op_PEenergyRange (LOW:HIGH, 5, 69) = (/    1.8810,    2.0900/)
  op_PEenergyRange (LOW:HIGH, 6, 69) = (/    2.0900,    2.3060/)
  op_PEenergyRange (LOW:HIGH, 7, 69) = (/    2.3060,    8.6480/)
  op_PEenergyRange (LOW:HIGH, 8, 69) = (/    8.6480,    9.6170/)
  op_PEenergyRange (LOW:HIGH, 9, 69) = (/    9.6170,   10.1160/)
  op_PEenergyRange (LOW:HIGH,10, 69) = (/   10.1160,   59.3900/)
  op_PEenergyRange (LOW:HIGH,11, 69) = (/   59.3900,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 69) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 70) = (/    0.0100,    0.1880/)
  op_PEenergyRange (LOW:HIGH, 2, 70) = (/    0.1880,    1.5280/)
  op_PEenergyRange (LOW:HIGH, 3, 70) = (/    1.5280,    1.5770/)
  op_PEenergyRange (LOW:HIGH, 4, 70) = (/    1.5770,    1.9500/)
  op_PEenergyRange (LOW:HIGH, 5, 70) = (/    1.9500,    2.1750/)
  op_PEenergyRange (LOW:HIGH, 6, 70) = (/    2.1750,    2.3980/)
  op_PEenergyRange (LOW:HIGH, 7, 70) = (/    2.3980,    8.9430/)
  op_PEenergyRange (LOW:HIGH, 8, 70) = (/    8.9430,    9.9780/)
  op_PEenergyRange (LOW:HIGH, 9, 70) = (/    9.9780,   10.4890/)
  op_PEenergyRange (LOW:HIGH,10, 70) = (/   10.4890,   61.3320/)
  op_PEenergyRange (LOW:HIGH,11, 70) = (/   61.3320,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 70) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 71) = (/    0.0100,    0.1990/)
  op_PEenergyRange (LOW:HIGH, 2, 71) = (/    0.1990,    1.5910/)
  op_PEenergyRange (LOW:HIGH, 3, 71) = (/    1.5910,    1.6410/)
  op_PEenergyRange (LOW:HIGH, 4, 71) = (/    1.6410,    2.0240/)
  op_PEenergyRange (LOW:HIGH, 5, 71) = (/    2.0240,    2.2640/)
  op_PEenergyRange (LOW:HIGH, 6, 71) = (/    2.2640,    2.4940/)
  op_PEenergyRange (LOW:HIGH, 7, 71) = (/    2.4940,    9.2450/)
  op_PEenergyRange (LOW:HIGH, 8, 71) = (/    9.2450,   10.3490/)
  op_PEenergyRange (LOW:HIGH, 9, 71) = (/   10.3490,   10.8740/)
  op_PEenergyRange (LOW:HIGH,10, 71) = (/   10.8740,   63.3160/)
  op_PEenergyRange (LOW:HIGH,11, 71) = (/   63.3160,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 71) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 72) = (/    0.0100,    0.2130/)
  op_PEenergyRange (LOW:HIGH, 2, 72) = (/    0.2130,    1.6620/)
  op_PEenergyRange (LOW:HIGH, 3, 72) = (/    1.6620,    1.7160/)
  op_PEenergyRange (LOW:HIGH, 4, 72) = (/    1.7160,    2.1080/)
  op_PEenergyRange (LOW:HIGH, 5, 72) = (/    2.1080,    2.3640/)
  op_PEenergyRange (LOW:HIGH, 6, 72) = (/    2.3640,    2.6000/)
  op_PEenergyRange (LOW:HIGH, 7, 72) = (/    2.6000,    9.5600/)
  op_PEenergyRange (LOW:HIGH, 8, 72) = (/    9.5600,   10.7390/)
  op_PEenergyRange (LOW:HIGH, 9, 72) = (/   10.7390,   11.2720/)
  op_PEenergyRange (LOW:HIGH,10, 72) = (/   11.2720,   65.3450/)
  op_PEenergyRange (LOW:HIGH,11, 72) = (/   65.3450,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 72) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 73) = (/    0.0100,    0.2280/)
  op_PEenergyRange (LOW:HIGH, 2, 73) = (/    0.2280,    1.7350/)
  op_PEenergyRange (LOW:HIGH, 3, 73) = (/    1.7350,    1.7930/)
  op_PEenergyRange (LOW:HIGH, 4, 73) = (/    1.7930,    2.1940/)
  op_PEenergyRange (LOW:HIGH, 5, 73) = (/    2.1940,    2.4690/)
  op_PEenergyRange (LOW:HIGH, 6, 73) = (/    2.4690,    2.7090/)
  op_PEenergyRange (LOW:HIGH, 7, 73) = (/    2.7090,    9.8800/)
  op_PEenergyRange (LOW:HIGH, 8, 73) = (/    9.8800,   11.1360/)
  op_PEenergyRange (LOW:HIGH, 9, 73) = (/   11.1360,   11.6800/)
  op_PEenergyRange (LOW:HIGH,10, 73) = (/   11.6800,   67.4160/)
  op_PEenergyRange (LOW:HIGH,11, 73) = (/   67.4160,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 73) = (/  500.0000,10000.0000/)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setPEenergyRangeAtomsZ73
