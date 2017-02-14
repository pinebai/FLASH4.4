!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setPEenergyRangeAtomsZ100
!!
!! NAME
!!
!!  op_setPEenergyRangeAtomsZ100
!!
!! SYNOPSIS
!!
!!  call op_setPEenergyRangeAtomsZ100 ()
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
!!  The atomic elements set in this routine are in the range: 74 =< Z =< 100.
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEenergyRangeAtomsZ100 ()

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
       call Driver_abortFlash ('[op_setPEenergyRangeAtomsZ100] ERROR: array op_PEenergyRange not allocated')
  end if

  goodShape =       (size (op_PEenergyRange,1) == HIGH-LOW+1)   &
              .and. (size (op_PEenergyRange,2) == 13)           &
              .and. (size (op_PEenergyRange,3) == 100)

  if (.not.goodShape) then
       call Driver_abortFlash ('[op_setPEenergyRangeAtomsZ100] ERROR: array op_PEenergyRange has bad shape')
  end if
!
!
!   ...Set the energy range array. Last index denotes the
!      atomic order number Z.
!
!
  op_PEenergyRange (LOW:HIGH, 1, 74) = (/    0.0100,    0.2420/)
  op_PEenergyRange (LOW:HIGH, 2, 74) = (/    0.2420,    1.8090/)
  op_PEenergyRange (LOW:HIGH, 3, 74) = (/    1.8090,    1.8710/)
  op_PEenergyRange (LOW:HIGH, 4, 74) = (/    1.8710,    2.2810/)
  op_PEenergyRange (LOW:HIGH, 5, 74) = (/    2.2810,    2.5750/)
  op_PEenergyRange (LOW:HIGH, 6, 74) = (/    2.5750,    2.8200/)
  op_PEenergyRange (LOW:HIGH, 7, 74) = (/    2.8200,   10.2040/)
  op_PEenergyRange (LOW:HIGH, 8, 74) = (/   10.2040,   11.5410/)
  op_PEenergyRange (LOW:HIGH, 9, 74) = (/   11.5410,   12.0980/)
  op_PEenergyRange (LOW:HIGH,10, 74) = (/   12.0980,   69.5250/)
  op_PEenergyRange (LOW:HIGH,11, 74) = (/   69.5250,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 74) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 75) = (/    0.0100,    0.2580/)
  op_PEenergyRange (LOW:HIGH, 2, 75) = (/    0.2580,    1.8830/)
  op_PEenergyRange (LOW:HIGH, 3, 75) = (/    1.8830,    1.9500/)
  op_PEenergyRange (LOW:HIGH, 4, 75) = (/    1.9500,    2.3680/)
  op_PEenergyRange (LOW:HIGH, 5, 75) = (/    2.3680,    2.6820/)
  op_PEenergyRange (LOW:HIGH, 6, 75) = (/    2.6820,    2.9340/)
  op_PEenergyRange (LOW:HIGH, 7, 75) = (/    2.9340,   10.5340/)
  op_PEenergyRange (LOW:HIGH, 8, 75) = (/   10.5340,   11.9570/)
  op_PEenergyRange (LOW:HIGH, 9, 75) = (/   11.9570,   12.5280/)
  op_PEenergyRange (LOW:HIGH,10, 75) = (/   12.5280,   71.6760/)
  op_PEenergyRange (LOW:HIGH,11, 75) = (/   71.6760,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 75) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 76) = (/    0.0100,    0.2740/)
  op_PEenergyRange (LOW:HIGH, 2, 76) = (/    0.2740,    1.9600/)
  op_PEenergyRange (LOW:HIGH, 3, 76) = (/    1.9600,    2.0310/)
  op_PEenergyRange (LOW:HIGH, 4, 76) = (/    2.0310,    2.4570/)
  op_PEenergyRange (LOW:HIGH, 5, 76) = (/    2.4570,    2.7920/)
  op_PEenergyRange (LOW:HIGH, 6, 76) = (/    2.7920,    3.0520/)
  op_PEenergyRange (LOW:HIGH, 7, 76) = (/    3.0520,   10.8710/)
  op_PEenergyRange (LOW:HIGH, 8, 76) = (/   10.8710,   12.3850/)
  op_PEenergyRange (LOW:HIGH, 9, 76) = (/   12.3850,   12.9690/)
  op_PEenergyRange (LOW:HIGH,10, 76) = (/   12.9690,   73.8710/)
  op_PEenergyRange (LOW:HIGH,11, 76) = (/   73.8710,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 76) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 77) = (/    0.0100,    0.2930/)
  op_PEenergyRange (LOW:HIGH, 2, 77) = (/    0.2930,    2.0400/)
  op_PEenergyRange (LOW:HIGH, 3, 77) = (/    2.0400,    2.1160/)
  op_PEenergyRange (LOW:HIGH, 4, 77) = (/    2.1160,    2.5510/)
  op_PEenergyRange (LOW:HIGH, 5, 77) = (/    2.5510,    2.9080/)
  op_PEenergyRange (LOW:HIGH, 6, 77) = (/    2.9080,    3.1730/)
  op_PEenergyRange (LOW:HIGH, 7, 77) = (/    3.1730,   11.2150/)
  op_PEenergyRange (LOW:HIGH, 8, 77) = (/   11.2150,   12.8240/)
  op_PEenergyRange (LOW:HIGH, 9, 77) = (/   12.8240,   13.4190/)
  op_PEenergyRange (LOW:HIGH,10, 77) = (/   13.4190,   76.1110/)
  op_PEenergyRange (LOW:HIGH,11, 77) = (/   76.1110,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 77) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 78) = (/    0.0100,    0.3130/)
  op_PEenergyRange (LOW:HIGH, 2, 78) = (/    0.3130,    2.1220/)
  op_PEenergyRange (LOW:HIGH, 3, 78) = (/    2.1220,    2.2020/)
  op_PEenergyRange (LOW:HIGH, 4, 78) = (/    2.2020,    2.6450/)
  op_PEenergyRange (LOW:HIGH, 5, 78) = (/    2.6450,    3.0270/)
  op_PEenergyRange (LOW:HIGH, 6, 78) = (/    3.0270,    3.2970/)
  op_PEenergyRange (LOW:HIGH, 7, 78) = (/    3.2970,   11.5640/)
  op_PEenergyRange (LOW:HIGH, 8, 78) = (/   11.5640,   13.2730/)
  op_PEenergyRange (LOW:HIGH, 9, 78) = (/   13.2730,   13.8800/)
  op_PEenergyRange (LOW:HIGH,10, 78) = (/   13.8800,   78.3950/)
  op_PEenergyRange (LOW:HIGH,11, 78) = (/   78.3950,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 78) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 79) = (/    0.0100,    0.3340/)
  op_PEenergyRange (LOW:HIGH, 2, 79) = (/    0.3340,    2.2060/)
  op_PEenergyRange (LOW:HIGH, 3, 79) = (/    2.2060,    2.2910/)
  op_PEenergyRange (LOW:HIGH, 4, 79) = (/    2.2910,    2.7430/)
  op_PEenergyRange (LOW:HIGH, 5, 79) = (/    2.7430,    3.1500/)
  op_PEenergyRange (LOW:HIGH, 6, 79) = (/    3.1500,    3.4250/)
  op_PEenergyRange (LOW:HIGH, 7, 79) = (/    3.4250,   11.9190/)
  op_PEenergyRange (LOW:HIGH, 8, 79) = (/   11.9190,   13.7340/)
  op_PEenergyRange (LOW:HIGH, 9, 79) = (/   13.7340,   14.3530/)
  op_PEenergyRange (LOW:HIGH,10, 79) = (/   14.3530,   80.7250/)
  op_PEenergyRange (LOW:HIGH,11, 79) = (/   80.7250,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 79) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 80) = (/    0.0100,    0.3600/)
  op_PEenergyRange (LOW:HIGH, 2, 80) = (/    0.3600,    2.2950/)
  op_PEenergyRange (LOW:HIGH, 3, 80) = (/    2.2950,    2.3850/)
  op_PEenergyRange (LOW:HIGH, 4, 80) = (/    2.3850,    2.8470/)
  op_PEenergyRange (LOW:HIGH, 5, 80) = (/    2.8470,    3.2800/)
  op_PEenergyRange (LOW:HIGH, 6, 80) = (/    3.2800,    3.5620/)
  op_PEenergyRange (LOW:HIGH, 7, 80) = (/    3.5620,   12.2830/)
  op_PEenergyRange (LOW:HIGH, 8, 80) = (/   12.2830,   14.2090/)
  op_PEenergyRange (LOW:HIGH, 9, 80) = (/   14.2090,   14.8420/)
  op_PEenergyRange (LOW:HIGH,10, 80) = (/   14.8420,   83.1020/)
  op_PEenergyRange (LOW:HIGH,11, 80) = (/   83.1020,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 80) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 81) = (/    0.0100,    0.3860/)
  op_PEenergyRange (LOW:HIGH, 2, 81) = (/    0.3860,    2.3890/)
  op_PEenergyRange (LOW:HIGH, 3, 81) = (/    2.3890,    2.4850/)
  op_PEenergyRange (LOW:HIGH, 4, 81) = (/    2.4850,    2.9560/)
  op_PEenergyRange (LOW:HIGH, 5, 81) = (/    2.9560,    3.4160/)
  op_PEenergyRange (LOW:HIGH, 6, 81) = (/    3.4160,    3.7040/)
  op_PEenergyRange (LOW:HIGH, 7, 81) = (/    3.7040,   12.6560/)
  op_PEenergyRange (LOW:HIGH, 8, 81) = (/   12.6560,   14.6970/)
  op_PEenergyRange (LOW:HIGH, 9, 81) = (/   14.6970,   15.3460/)
  op_PEenergyRange (LOW:HIGH,10, 81) = (/   15.3460,   85.5300/)
  op_PEenergyRange (LOW:HIGH,11, 81) = (/   85.5300,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 81) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 82) = (/    0.0100,    0.4130/)
  op_PEenergyRange (LOW:HIGH, 2, 82) = (/    0.4130,    2.4840/)
  op_PEenergyRange (LOW:HIGH, 3, 82) = (/    2.4840,    2.5860/)
  op_PEenergyRange (LOW:HIGH, 4, 82) = (/    2.5860,    3.0660/)
  op_PEenergyRange (LOW:HIGH, 5, 82) = (/    3.0660,    3.5540/)
  op_PEenergyRange (LOW:HIGH, 6, 82) = (/    3.5540,    3.8510/)
  op_PEenergyRange (LOW:HIGH, 7, 82) = (/    3.8510,   13.0350/)
  op_PEenergyRange (LOW:HIGH, 8, 82) = (/   13.0350,   15.2000/)
  op_PEenergyRange (LOW:HIGH, 9, 82) = (/   15.2000,   15.8610/)
  op_PEenergyRange (LOW:HIGH,10, 82) = (/   15.8610,   88.0040/)
  op_PEenergyRange (LOW:HIGH,11, 82) = (/   88.0040,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 82) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 83) = (/    0.0100,    0.4400/)
  op_PEenergyRange (LOW:HIGH, 2, 83) = (/    0.4400,    2.5810/)
  op_PEenergyRange (LOW:HIGH, 3, 83) = (/    2.5810,    2.6890/)
  op_PEenergyRange (LOW:HIGH, 4, 83) = (/    2.6890,    3.1770/)
  op_PEenergyRange (LOW:HIGH, 5, 83) = (/    3.1770,    3.6960/)
  op_PEenergyRange (LOW:HIGH, 6, 83) = (/    3.6960,    4.0000/)
  op_PEenergyRange (LOW:HIGH, 7, 83) = (/    4.0000,   13.4200/)
  op_PEenergyRange (LOW:HIGH, 8, 83) = (/   13.4200,   15.7140/)
  op_PEenergyRange (LOW:HIGH, 9, 83) = (/   15.7140,   16.3910/)
  op_PEenergyRange (LOW:HIGH,10, 83) = (/   16.3910,   90.5260/)
  op_PEenergyRange (LOW:HIGH,11, 83) = (/   90.5260,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 83) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 84) = (/    0.0100,    0.4730/)
  op_PEenergyRange (LOW:HIGH, 2, 84) = (/    0.4730,    2.6830/)
  op_PEenergyRange (LOW:HIGH, 3, 84) = (/    2.6830,    2.7980/)
  op_PEenergyRange (LOW:HIGH, 4, 84) = (/    2.7980,    3.2950/)
  op_PEenergyRange (LOW:HIGH, 5, 84) = (/    3.2950,    3.8490/)
  op_PEenergyRange (LOW:HIGH, 6, 84) = (/    3.8490,    4.1560/)
  op_PEenergyRange (LOW:HIGH, 7, 84) = (/    4.1560,   13.8140/)
  op_PEenergyRange (LOW:HIGH, 8, 84) = (/   13.8140,   16.2440/)
  op_PEenergyRange (LOW:HIGH, 9, 84) = (/   16.2440,   16.9360/)
  op_PEenergyRange (LOW:HIGH,10, 84) = (/   16.9360,   93.1050/)
  op_PEenergyRange (LOW:HIGH,11, 84) = (/   93.1050,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 84) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 85) = (/    0.0100,    0.5060/)
  op_PEenergyRange (LOW:HIGH, 2, 85) = (/    0.5060,    2.7870/)
  op_PEenergyRange (LOW:HIGH, 3, 85) = (/    2.7870,    2.9090/)
  op_PEenergyRange (LOW:HIGH, 4, 85) = (/    2.9090,    3.4160/)
  op_PEenergyRange (LOW:HIGH, 5, 85) = (/    3.4160,    4.0060/)
  op_PEenergyRange (LOW:HIGH, 6, 85) = (/    4.0060,    4.3170/)
  op_PEenergyRange (LOW:HIGH, 7, 85) = (/    4.3170,   14.2140/)
  op_PEenergyRange (LOW:HIGH, 8, 85) = (/   14.2140,   16.7850/)
  op_PEenergyRange (LOW:HIGH, 9, 85) = (/   16.7850,   17.4910/)
  op_PEenergyRange (LOW:HIGH,10, 85) = (/   17.4910,   95.7300/)
  op_PEenergyRange (LOW:HIGH,11, 85) = (/   95.7300,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 85) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 86) = (/    0.0100,    0.5410/)
  op_PEenergyRange (LOW:HIGH, 2, 86) = (/    0.5410,    2.8920/)
  op_PEenergyRange (LOW:HIGH, 3, 86) = (/    2.8920,    3.0220/)
  op_PEenergyRange (LOW:HIGH, 4, 86) = (/    3.0220,    3.5380/)
  op_PEenergyRange (LOW:HIGH, 5, 86) = (/    3.5380,    4.1640/)
  op_PEenergyRange (LOW:HIGH, 6, 86) = (/    4.1640,    4.4820/)
  op_PEenergyRange (LOW:HIGH, 7, 86) = (/    4.4820,   14.6190/)
  op_PEenergyRange (LOW:HIGH, 8, 86) = (/   14.6190,   17.3370/)
  op_PEenergyRange (LOW:HIGH, 9, 86) = (/   17.3370,   18.0550/)
  op_PEenergyRange (LOW:HIGH,10, 86) = (/   18.0550,   98.4040/)
  op_PEenergyRange (LOW:HIGH,11, 86) = (/   98.4040,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 86) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 87) = (/    0.0100,    0.5760/)
  op_PEenergyRange (LOW:HIGH, 2, 87) = (/    0.5760,    3.0000/)
  op_PEenergyRange (LOW:HIGH, 3, 87) = (/    3.0000,    3.1360/)
  op_PEenergyRange (LOW:HIGH, 4, 87) = (/    3.1360,    3.6640/)
  op_PEenergyRange (LOW:HIGH, 5, 87) = (/    3.6640,    4.3250/)
  op_PEenergyRange (LOW:HIGH, 6, 87) = (/    4.3250,    4.6520/)
  op_PEenergyRange (LOW:HIGH, 7, 87) = (/    4.6520,   15.0300/)
  op_PEenergyRange (LOW:HIGH, 8, 87) = (/   15.0300,   17.9040/)
  op_PEenergyRange (LOW:HIGH, 9, 87) = (/   17.9040,   18.6390/)
  op_PEenergyRange (LOW:HIGH,10, 87) = (/   18.6390,  101.1370/)
  op_PEenergyRange (LOW:HIGH,11, 87) = (/  101.1370,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 87) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 88) = (/    0.0100,    0.6100/)
  op_PEenergyRange (LOW:HIGH, 2, 88) = (/    0.6100,    3.1090/)
  op_PEenergyRange (LOW:HIGH, 3, 88) = (/    3.1090,    3.2530/)
  op_PEenergyRange (LOW:HIGH, 4, 88) = (/    3.2530,    3.7910/)
  op_PEenergyRange (LOW:HIGH, 5, 88) = (/    3.7910,    4.4900/)
  op_PEenergyRange (LOW:HIGH, 6, 88) = (/    4.4900,    4.8240/)
  op_PEenergyRange (LOW:HIGH, 7, 88) = (/    4.8240,   15.4460/)
  op_PEenergyRange (LOW:HIGH, 8, 88) = (/   15.4460,   18.4840/)
  op_PEenergyRange (LOW:HIGH, 9, 88) = (/   18.4840,   19.2370/)
  op_PEenergyRange (LOW:HIGH,10, 88) = (/   19.2370,  103.9220/)
  op_PEenergyRange (LOW:HIGH,11, 88) = (/  103.9220,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 88) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 89) = (/    0.0100,    0.6440/)
  op_PEenergyRange (LOW:HIGH, 2, 89) = (/    0.6440,    3.2190/)
  op_PEenergyRange (LOW:HIGH, 3, 89) = (/    3.2190,    3.3710/)
  op_PEenergyRange (LOW:HIGH, 4, 89) = (/    3.3710,    3.9180/)
  op_PEenergyRange (LOW:HIGH, 5, 89) = (/    3.9180,    4.6580/)
  op_PEenergyRange (LOW:HIGH, 6, 89) = (/    4.6580,    5.0020/)
  op_PEenergyRange (LOW:HIGH, 7, 89) = (/    5.0020,   15.8700/)
  op_PEenergyRange (LOW:HIGH, 8, 89) = (/   15.8700,   19.0830/)
  op_PEenergyRange (LOW:HIGH, 9, 89) = (/   19.0830,   19.8450/)
  op_PEenergyRange (LOW:HIGH,10, 89) = (/   19.8450,  106.7590/)
  op_PEenergyRange (LOW:HIGH,11, 89) = (/  106.7590,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 89) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 90) = (/    0.0100,    0.6760/)
  op_PEenergyRange (LOW:HIGH, 2, 90) = (/    0.6760,    3.3320/)
  op_PEenergyRange (LOW:HIGH, 3, 90) = (/    3.3320,    3.4900/)
  op_PEenergyRange (LOW:HIGH, 4, 90) = (/    3.4900,    4.0460/)
  op_PEenergyRange (LOW:HIGH, 5, 90) = (/    4.0460,    4.8300/)
  op_PEenergyRange (LOW:HIGH, 6, 90) = (/    4.8300,    5.1820/)
  op_PEenergyRange (LOW:HIGH, 7, 90) = (/    5.1820,   16.3000/)
  op_PEenergyRange (LOW:HIGH, 8, 90) = (/   16.3000,   19.6930/)
  op_PEenergyRange (LOW:HIGH, 9, 90) = (/   19.6930,   20.4660/)
  op_PEenergyRange (LOW:HIGH,10, 90) = (/   20.4660,  109.6510/)
  op_PEenergyRange (LOW:HIGH,11, 90) = (/  109.6510,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 90) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 91) = (/    0.0100,    0.7060/)
  op_PEenergyRange (LOW:HIGH, 2, 91) = (/    0.7060,    3.4420/)
  op_PEenergyRange (LOW:HIGH, 3, 91) = (/    3.4420,    3.6090/)
  op_PEenergyRange (LOW:HIGH, 4, 91) = (/    3.6090,    4.1740/)
  op_PEenergyRange (LOW:HIGH, 5, 91) = (/    4.1740,    5.0030/)
  op_PEenergyRange (LOW:HIGH, 6, 91) = (/    5.0030,    5.3640/)
  op_PEenergyRange (LOW:HIGH, 7, 91) = (/    5.3640,   16.7330/)
  op_PEenergyRange (LOW:HIGH, 8, 91) = (/   16.7330,   20.3140/)
  op_PEenergyRange (LOW:HIGH, 9, 91) = (/   20.3140,   21.1050/)
  op_PEenergyRange (LOW:HIGH,10, 91) = (/   21.1050,  112.6010/)
  op_PEenergyRange (LOW:HIGH,11, 91) = (/  112.6010,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 91) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 92) = (/    0.0100,    0.7400/)
  op_PEenergyRange (LOW:HIGH, 2, 92) = (/    0.7400,    3.5520/)
  op_PEenergyRange (LOW:HIGH, 3, 92) = (/    3.5520,    3.7280/)
  op_PEenergyRange (LOW:HIGH, 4, 92) = (/    3.7280,    4.3040/)
  op_PEenergyRange (LOW:HIGH, 5, 92) = (/    4.3040,    5.1810/)
  op_PEenergyRange (LOW:HIGH, 6, 92) = (/    5.1810,    5.5480/)
  op_PEenergyRange (LOW:HIGH, 7, 92) = (/    5.5480,   17.1700/)
  op_PEenergyRange (LOW:HIGH, 8, 92) = (/   17.1700,   20.9480/)
  op_PEenergyRange (LOW:HIGH, 9, 92) = (/   20.9480,   21.7590/)
  op_PEenergyRange (LOW:HIGH,10, 92) = (/   21.7590,  115.6060/)
  op_PEenergyRange (LOW:HIGH,11, 92) = (/  115.6060,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 92) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 93) = (/    0.0100,    0.7720/)
  op_PEenergyRange (LOW:HIGH, 2, 93) = (/    0.7720,    3.6640/)
  op_PEenergyRange (LOW:HIGH, 3, 93) = (/    3.6640,    3.8500/)
  op_PEenergyRange (LOW:HIGH, 4, 93) = (/    3.8500,    4.4350/)
  op_PEenergyRange (LOW:HIGH, 5, 93) = (/    4.4350,    5.3660/)
  op_PEenergyRange (LOW:HIGH, 6, 93) = (/    5.3660,    5.7350/)
  op_PEenergyRange (LOW:HIGH, 7, 93) = (/    5.7350,   17.6130/)
  op_PEenergyRange (LOW:HIGH, 8, 93) = (/   17.6130,   21.6000/)
  op_PEenergyRange (LOW:HIGH, 9, 93) = (/   21.6000,   22.4270/)
  op_PEenergyRange (LOW:HIGH,10, 93) = (/   22.4270,  118.6700/)
  op_PEenergyRange (LOW:HIGH,11, 93) = (/  118.6700,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 93) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 94) = (/    0.0100,    0.8010/)
  op_PEenergyRange (LOW:HIGH, 2, 94) = (/    0.8010,    3.7780/)
  op_PEenergyRange (LOW:HIGH, 3, 94) = (/    3.7780,    3.9730/)
  op_PEenergyRange (LOW:HIGH, 4, 94) = (/    3.9730,    4.5680/)
  op_PEenergyRange (LOW:HIGH, 5, 94) = (/    4.5680,    5.5550/)
  op_PEenergyRange (LOW:HIGH, 6, 94) = (/    5.5550,    5.9270/)
  op_PEenergyRange (LOW:HIGH, 7, 94) = (/    5.9270,   18.0630/)
  op_PEenergyRange (LOW:HIGH, 8, 94) = (/   18.0630,   22.2700/)
  op_PEenergyRange (LOW:HIGH, 9, 94) = (/   22.2700,   23.1090/)
  op_PEenergyRange (LOW:HIGH,10, 94) = (/   23.1090,  121.7970/)
  op_PEenergyRange (LOW:HIGH,11, 94) = (/  121.7970,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 94) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 95) = (/    0.0100,    0.8280/)
  op_PEenergyRange (LOW:HIGH, 2, 95) = (/    0.8280,    3.8940/)
  op_PEenergyRange (LOW:HIGH, 3, 95) = (/    3.8940,    4.1000/)
  op_PEenergyRange (LOW:HIGH, 4, 95) = (/    4.1000,    4.7030/)
  op_PEenergyRange (LOW:HIGH, 5, 95) = (/    4.7030,    5.7480/)
  op_PEenergyRange (LOW:HIGH, 6, 95) = (/    5.7480,    6.1220/)
  op_PEenergyRange (LOW:HIGH, 7, 95) = (/    6.1220,   18.5190/)
  op_PEenergyRange (LOW:HIGH, 8, 95) = (/   18.5190,   22.9580/)
  op_PEenergyRange (LOW:HIGH, 9, 95) = (/   22.9580,   23.8120/)
  op_PEenergyRange (LOW:HIGH,10, 95) = (/   23.8120,  124.9900/)
  op_PEenergyRange (LOW:HIGH,11, 95) = (/  124.9900,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 95) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 96) = (/    0.0100,    0.8530/)
  op_PEenergyRange (LOW:HIGH, 2, 96) = (/    0.8530,    4.0120/)
  op_PEenergyRange (LOW:HIGH, 3, 96) = (/    4.0120,    4.2300/)
  op_PEenergyRange (LOW:HIGH, 4, 96) = (/    4.2300,    4.8390/)
  op_PEenergyRange (LOW:HIGH, 5, 96) = (/    4.8390,    5.9450/)
  op_PEenergyRange (LOW:HIGH, 6, 96) = (/    5.9450,    6.3220/)
  op_PEenergyRange (LOW:HIGH, 7, 96) = (/    6.3220,   18.9820/)
  op_PEenergyRange (LOW:HIGH, 8, 96) = (/   18.9820,   23.6630/)
  op_PEenergyRange (LOW:HIGH, 9, 96) = (/   23.6630,   24.5350/)
  op_PEenergyRange (LOW:HIGH,10, 96) = (/   24.5350,  128.2530/)
  op_PEenergyRange (LOW:HIGH,11, 96) = (/  128.2530,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 96) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 97) = (/    0.0100,    0.8770/)
  op_PEenergyRange (LOW:HIGH, 2, 97) = (/    0.8770,    4.1320/)
  op_PEenergyRange (LOW:HIGH, 3, 97) = (/    4.1320,    4.3640/)
  op_PEenergyRange (LOW:HIGH, 4, 97) = (/    4.3640,    4.9770/)
  op_PEenergyRange (LOW:HIGH, 5, 97) = (/    4.9770,    6.1470/)
  op_PEenergyRange (LOW:HIGH, 6, 97) = (/    6.1470,    6.5260/)
  op_PEenergyRange (LOW:HIGH, 7, 97) = (/    6.5260,   19.4520/)
  op_PEenergyRange (LOW:HIGH, 8, 97) = (/   19.4520,   24.3850/)
  op_PEenergyRange (LOW:HIGH, 9, 97) = (/   24.3850,   25.2750/)
  op_PEenergyRange (LOW:HIGH,10, 97) = (/   25.2750,  131.5900/)
  op_PEenergyRange (LOW:HIGH,11, 97) = (/  131.5900,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 97) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 98) = (/    0.0100,    0.9020/)
  op_PEenergyRange (LOW:HIGH, 2, 98) = (/    0.9020,    4.2540/)
  op_PEenergyRange (LOW:HIGH, 3, 98) = (/    4.2540,    4.5020/)
  op_PEenergyRange (LOW:HIGH, 4, 98) = (/    4.5020,    5.1170/)
  op_PEenergyRange (LOW:HIGH, 5, 98) = (/    5.1170,    6.3530/)
  op_PEenergyRange (LOW:HIGH, 6, 98) = (/    6.3530,    6.7350/)
  op_PEenergyRange (LOW:HIGH, 7, 98) = (/    6.7350,   19.9290/)
  op_PEenergyRange (LOW:HIGH, 8, 98) = (/   19.9290,   25.1250/)
  op_PEenergyRange (LOW:HIGH, 9, 98) = (/   25.1250,   26.0300/)
  op_PEenergyRange (LOW:HIGH,10, 98) = (/   26.0300,  135.0050/)
  op_PEenergyRange (LOW:HIGH,11, 98) = (/  135.0050,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 98) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1, 99) = (/    0.0100,    0.9270/)
  op_PEenergyRange (LOW:HIGH, 2, 99) = (/    0.9270,    4.3780/)
  op_PEenergyRange (LOW:HIGH, 3, 99) = (/    4.3780,    4.6440/)
  op_PEenergyRange (LOW:HIGH, 4, 99) = (/    4.6440,    5.2590/)
  op_PEenergyRange (LOW:HIGH, 5, 99) = (/    5.2590,    6.5640/)
  op_PEenergyRange (LOW:HIGH, 6, 99) = (/    6.5640,    6.9490/)
  op_PEenergyRange (LOW:HIGH, 7, 99) = (/    6.9490,   20.4140/)
  op_PEenergyRange (LOW:HIGH, 8, 99) = (/   20.4140,   25.8830/)
  op_PEenergyRange (LOW:HIGH, 9, 99) = (/   25.8830,   26.8030/)
  op_PEenergyRange (LOW:HIGH,10, 99) = (/   26.8030,  138.5020/)
  op_PEenergyRange (LOW:HIGH,11, 99) = (/  138.5020,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12, 99) = (/  500.0000,10000.0000/)
 
  op_PEenergyRange (LOW:HIGH, 1,100) = (/    0.0100,    0.9520/)
  op_PEenergyRange (LOW:HIGH, 2,100) = (/    0.9520,    4.5040/)
  op_PEenergyRange (LOW:HIGH, 3,100) = (/    4.5040,    4.7900/)
  op_PEenergyRange (LOW:HIGH, 4,100) = (/    4.7900,    5.4030/)
  op_PEenergyRange (LOW:HIGH, 5,100) = (/    5.4030,    6.7800/)
  op_PEenergyRange (LOW:HIGH, 6,100) = (/    6.7800,    7.1680/)
  op_PEenergyRange (LOW:HIGH, 7,100) = (/    7.1680,   20.9070/)
  op_PEenergyRange (LOW:HIGH, 8,100) = (/   20.9070,   26.6590/)
  op_PEenergyRange (LOW:HIGH, 9,100) = (/   26.6590,   27.5940/)
  op_PEenergyRange (LOW:HIGH,10,100) = (/   27.5940,  142.0850/)
  op_PEenergyRange (LOW:HIGH,11,100) = (/  142.0850,  500.0000/)
  op_PEenergyRange (LOW:HIGH,12,100) = (/  500.0000,10000.0000/)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setPEenergyRangeAtomsZ100
