!!****if* source/Multispecies/MultispeciesMain/unitTest/Multispecies_unitTest
!!
!! NAME
!!
!!  Multispecies_unitTest
!!
!! SYNOPSIS
!!
!!  Multispecies_unitTest(integer, intent(in)::fileUnit,
!!                        logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the Multispecies
!!  unit.  Normally called from Driver_evolveFlash within a Simulation unitTest
!!  See, for example, source/Simulation/SimulationMain/unitTest/Multispecies
!!
!!  Performs testing by calculating and outputting values using all possible
!!  arguments.  Also shows some erroneous calls (commented out).
!!
!! ARGUMENTS
!!
!!    fileUnit  -- number of file unit for diagnostic output
!!    perfect   -- if .true., unitTest has returned correctly
!! 
!! NOTES
!!  In the Simulation unit, you must set your Config file to 
!!    REQUIRES Multispecies/MultispeciesMain/unitTest
!!
!!***

subroutine Multispecies_unitTest(fileUnit,perfect)

  use Multispecies_data !, ONLY:   multispecies_type, ms_Array
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Multispecies_interface, ONLY : Multispecies_list, &
       Multispecies_getProperty, Multispecies_setProperty, &
       Multispecies_getAvg, Multispecies_getSum, Multispecies_getSumFrac, &
       Multispecies_getSumInv, Multispecies_getSumSqr

  implicit none

  !! has to follow these statements because it includes an interface block
#include "constants.h"
#include "Flash.h"  !! defines species names
#include "Multispecies.h"

  integer, intent(in)         :: fileUnit
  logical, intent(inout)      :: perfect
  !------------------------------------------------------------------------
  type(multispecies_type)     :: singleMS
  integer                     :: i
  character(len=4)            :: speciesName
  real                        :: value, valueAvg,valueSum,                    &
       &          valueSumFrac,valueSumInv,valueSumSqr, valueEB, valueE
  real, dimension(1)          :: weightsSingle
  real, dimension(NSPECIES)   :: weightsAll
  integer, dimension(NSPECIES):: maskAll
  real                        :: weightsOne=1.0
  real, dimension(2)          :: weightsDimTwo
  real, dimension(4)          :: weightsTwo
  integer, dimension(4)       :: maskOne,maskTwo
  integer                     :: maskSingle

  !------------------------------------------------------------------
  ! this diagnostic outputs the entire array -- there are 4 members = NSPECIES

  ! Initialize parameters in Simulation_init()    
  !! Output their values here
  perfect = .true.
  write(fileUnit,900)NSPECIES
  do i=1,NSPECIES
     singleMS = ms_Array(i)
     !! Since species are defined in the Config file, they are actually NUMBERS
     call Simulation_mapIntToStr(singleMS%name,speciesName,MAPBLOCK_UNK)
     write(fileUnit,901)i,speciesName, &
          ms_Array(i)%name, &
          ms_Array(i)%numTotal, &
          ms_Array(i)%numPositive, &
          ms_Array(i)%numNeutral, &
          ms_Array(i)%numNegative, &
          ms_Array(i)%bindingEnergy, &
          ms_Array(i)%adiabaticIndex, &
          ms_Array(i)%eosType, &
          ms_Array(i)%eosZfreeTableFile, &
          ms_Array(i)%eosEnerTableFile
  enddo
  write(fileUnit,*)

  !!  Now see if internal routine gives same results
  write (fileUnit,*)'Now results from Multispecies_list'
  call Multispecies_list(fileUnit)
  !! test ability to get the properties
  write(fileUnit,910)
  call Multispecies_getProperty(BIG_SPEC,A,value)
  write(fileUnit,911)"A",value
  call Multispecies_getProperty(BIG_SPEC,Z,value)
  write(fileUnit,911)"Z",value
  call Multispecies_getProperty(BIG_SPEC,N,value)
  write(fileUnit,911)"N",value
  call Multispecies_getProperty(BIG_SPEC,E,value)
  write(fileUnit,911)"E",value
  call Multispecies_getProperty(BIG_SPEC,EB,value)
  write(fileUnit,911)"EB",value
  call Multispecies_getProperty(BIG_SPEC,GAMMA,value)
  write(fileUnit,911)"GAMMA",value

  !! test ability to set properties
  call Multispecies_setProperty(BIG_SPEC,E,45.0)
  write(fileUnit,912)
  call Multispecies_getProperty(BIG_SPEC,E,valueE)
  write(fileUnit,911)"E",valueE
  !! Note!  The interface block DOES catch this erroneous (integer) call
  !       call Multispecies_setProperty(BIG_SPEC,EB,45)
  !       call Multispecies_getProperty(BIG_SPEC,EB,valueEB)
  !       write(fileUnit,911)"EB",valueEB
  if (abs(valueE-45.0).gt. 1.0e-6) perfect = .false.

  ! Now test all get routines without masks or weights
  do i=1,NSPECIES
     weightsAll(i) = 1.0
  enddo
  write(fileUnit,920)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg)
  call Multispecies_getSum(A,valueSum)
  call Multispecies_getSumFrac(A,valueSumFrac)
  call Multispecies_getSumInv(A,valueSumInv)
  call Multispecies_getSumSqr(A,valueSumSqr)
  write(fileUnit,922)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr

  ! Now test with mask of all species, but straight order
  maskAll(2)=SMAL_SPEC
  maskAll(4)=NEG_SPEC
  maskAll(3)=ZERO_SPEC
  maskAll(1)=BIG_SPEC
  write(fileUnit,823)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,speciesMask=maskAll)
  call Multispecies_getSum(A,valueSum,speciesMask=maskAll)
  call Multispecies_getSumFrac(A,valueSumFrac,speciesMask=maskAll)
  call Multispecies_getSumInv(A,valueSumInv,speciesMask=maskAll)
  call Multispecies_getSumSqr(A,valueSumSqr,speciesMask=maskAll)
  write(fileUnit,932)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr,maskAll


  ! Now test with mask of all species, but shift them around
  maskAll(1)=SMAL_SPEC
  maskAll(2)=NEG_SPEC
  maskAll(3)=ZERO_SPEC
  maskAll(4)=BIG_SPEC
  write(fileUnit,923)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,speciesMask=maskAll)
  call Multispecies_getSum(A,valueSum,speciesMask=maskAll)
  call Multispecies_getSumFrac(A,valueSumFrac,speciesMask=maskAll)
  call Multispecies_getSumInv(A,valueSumInv,speciesMask=maskAll)
  call Multispecies_getSumSqr(A,valueSumSqr,speciesMask=maskAll)
  write(fileUnit,932)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr,maskAll

  ! Now test with mask of two species
  maskTwo(1)=BIG_SPEC
  maskTwo(2)=SMAL_SPEC
  maskTwo(3)=UNDEFINED_INT
  maskTwo(4)=UNDEFINED_INT
  write(fileUnit,924)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,speciesMask=maskTwo)
  call Multispecies_getSum(A,valueSum,speciesMask=maskTwo)
  call Multispecies_getSumFrac(A,valueSumFrac,speciesMask=maskTwo)
  call Multispecies_getSumInv(A,valueSumInv,speciesMask=maskTwo)
  call Multispecies_getSumSqr(A,valueSumSqr,speciesMask=maskTwo)
  write(fileUnit,932)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr,maskTwo

  ! Now test with mask of just one species
  maskOne(1)=UNDEFINED_INT
  maskOne(2)=ZERO_SPEC
  maskOne(3)=UNDEFINED_INT
  maskOne(4)=UNDEFINED_INT
  write(fileUnit,925)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,speciesMask=maskOne)
  call Multispecies_getSum(A,valueSum,speciesMask=maskOne)
  call Multispecies_getSumFrac(A,valueSumFrac,speciesMask=maskOne)
  call Multispecies_getSumInv(A,valueSumInv,speciesMask=maskOne)
  call Multispecies_getSumSqr(A,valueSumSqr,speciesMask=maskOne)
  write(fileUnit,932)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr,maskOne

  ! Now test with varying weights:  zero for everything but one entry
  weightsAll(1) = 1.0 
  weightsAll(2) = 0.0  
  weightsAll(3) = 0.0
  weightsAll(4) = 0.0
  write(fileUnit,926)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,weights=weightsAll)
  call Multispecies_getSum(A,valueSum,weights=weightsAll)
  call Multispecies_getSumFrac(A,valueSumFrac,weights=weightsAll)
  call Multispecies_getSumInv(A,valueSumInv,weights=weightsAll)
  call Multispecies_getSumSqr(A,valueSumSqr,weights=weightsAll)
  write(fileUnit,942)valueAvg,valueSum,valueSumFrac,valueSumInv,          &
       &          valueSumSqr,weightsAll

  ! Now test with varying weights:  single value
  ! This failed in compilation ONLY on a SGI compiler, where the compiler
  !   complains that weights(NSPECIES) doesn't match weightsSingle
  ! "The overall size of the dummy argument array is greater than the size of this actual argument"
  ! That's a pretty valid complaint, although all other compilers/execution were happy with this mismatch.
  ! So code was changed to weights(:) to allow either size 1 or size NSPECIES

  weightsSingle = 2.0
  write(fileUnit,927)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,weights=weightsSingle)
  call Multispecies_getSum(A,valueSum,weights=weightsSingle)
  call Multispecies_getSumFrac(A,valueSumFrac,weights=weightsSingle)
  call Multispecies_getSumInv(A,valueSumInv,weights=weightsSingle)
  call Multispecies_getSumSqr(A,valueSumSqr,weights=weightsSingle)
  write(fileUnit,942)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr, weightsAll

  !  Weights and mask
  weightsAll(1) = 1.0 
  weightsAll(2) = 2.0  
  weightsAll(3) = 3.0
  weightsAll(4) = 4.0
  maskAll(1) = UNDEFINED_INT
  maskAll(2) = BIG_SPEC
  maskAll(3) = ZERO_SPEC
  maskAll(4) = UNDEFINED_INT
  write(fileUnit,928)
  write(fileUnit,921)
  call Multispecies_getAvg(A,valueAvg,speciesMask=maskAll,weights=weightsAll)
  call Multispecies_getSum(A,valueSum,speciesMask=maskAll,weights=weightsAll)
  call Multispecies_getSumFrac(A,valueSumFrac,speciesMask=maskAll,weights=weightsAll)
  call Multispecies_getSumInv(A,valueSumInv,speciesMask=maskAll,weights=weightsAll)
  call Multispecies_getSumSqr(A,valueSumSqr,speciesMask=maskAll,weights=weightsAll)
  write(fileUnit,942)valueAvg,valueSum,valueSumFrac,valueSumInv,              &
       &          valueSumSqr, weightsAll


  !! Should test here for erroneous calls
  write(fileUnit,930)
  write(fileUnit,*)"Call 1:  Incorrect species name (won't compile)"
  !!     call Multispecies_getProperty(NONE_SPEC,A)
  write(fileUnit,*)"Call 2:  Incorrect property name (won't compile)"
  !!     call Multispecies_getProperty(BIG_SPEC,Q)
  write(fileUnit,*)"Call 3:  Incorrect mask dimension (won't compile)"
  maskSingle = 0
  !     call Multispecies_getSum(A,valueSum,weightsAll,maskSingle)
  write(fileUnit,*)"Call 4:  Incorrect weights dimension (returns error msg)"
  weightsDimTwo(1) = 1.0
  weightsDimTwo(2) = 1.0
  !     call Multispecies_getSum(A,valueSum,weightsDimTwo,maskAll)


  print*,'Done with testing Multispecies'

  return 

  !------------------------------------------------------------------------

900 format("Within Multispecies_unitTest, initial values of ",I5," species",/,        &
       &         "I",T5,"Name",T10,"Index",T25,                                           &
       &         "Total",T35,"Positive",T45,"Neutral",T55,"Negative",T65,         &
       &         "bind Ener",T75,"Gamma",4x,"eosType", &
                 2(1x,"table file                      "))
901 format(I3,T5,A4,T15,I5,5X,6(ES9.2,1x), &
         I6, &
         2(1X,A32))
910 format(/,"Testing Multispecies_getProperty for BIG_SPEC")
911 format("Property ",A5," has value ",ES9.2)
912 format(/,"Testing Multispecies_setProperty for BIG_SPEC")
920 format(/,"Testing Multispecies_get* for all species, A, no masks, no weights")
921 format("Average",T10,"Sum",T20,"SumFrac",T30,"SumInv",T40,"SumSqr")
922 format(5(ES9.2),1x)
932 format( (5(ES9.2),1x),/,(4(I5,5x)))
942 format( (5(ES9.2),1x),/,(4(ES9.2,1x)))
823 format(/,"Testing Multispecies_get* for A, with straight mask for all")
923 format(/,"Testing Multispecies_get* for A, with shifted mask for all")
924 format(/,"Testing Multispecies_get* for BIG,SMALL, A, with mask for two")
925 format(/,"Testing Multispecies_get* for ZERO, A, with mask for one")
926 format(/,"Testing Multispecies_get* for A, with no mask, weight for BIG_SPEC")
927 format(/,"Testing Multispecies_get* for A, with no mask, weight of 2.0")
928 format(/,"Testing Multispecies_get* for A, with mask, and weights")
930 format(/,"Testing for erroneous calls")

end subroutine Multispecies_unitTest
