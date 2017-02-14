!!****if* source/PhysicalConstants/PhysicalConstantsMain/unitTest/PhysicalConstants_unitTest
!!
!! NAME
!!
!!  PhysicalConstants_unitTest
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_unitTest(integer, intent(in)::fileUnit,
!!                             logical, intent(out)::perfect  )
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the PhysicalConstants
!!      unit.  Normally called from Driver_evolveFlash within a Simulation unitTest
!!  See, for example, source/Simulation/SimulationMain/unitTest/PhysConst
!!
!! ARGUMENTS
!!      
!!     fileUnit - integer(in)   file for diagnostic output
!!     perfect  - logical(out)  TRUE if all tests are passed
!!
!! NOTES
!!
!!***

subroutine PhysicalConstants_unitTest(fileUnit,perfect)

  use PhysicalConstants_data, ONLY:   pc_nameSISystem, pc_SISystem, pc_globalMe
  use Driver_interface, ONLY : Driver_abortFlash
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get, &
       PhysicalConstants_list, PhysicalConstants_listUnits
  implicit none

  !! has to follow these statements because it includes an interface block
#include "constants.h"

  
  integer, intent(in)         :: fileUnit
  logical, intent(out)        :: perfect
  !------------------------------------------------------------------------
  integer                     ::  i
  real                        ::  value
  character(len=MAX_STRING_LENGTH) :: errorstring
  integer                     ::  isError
  
  
  perfect = .true.

  !------------------------------------------------------------------
  ! Initialize database in Simulation_init()    

  write(fileUnit,901)pc_nameSISystem
  write(fileUnit,902)pc_SISystem
  ! Test database setup
  if (pc_globalMe.EQ.MASTER_PE)  print *,'Now listing Units and Constants databases'


  write(fileUnit, 906)
  call PhysicalConstants_listUnits(fileUnit)

  write(fileUnit, 922)
  call PhysicalConstants_list(fileUnit)


  !  Now test getting constants with different units
  if (pc_globalMe.EQ.MASTER_PE) print *,'Now calling PhysicalConstants_get with various parameters'
  write(fileUnit,910)"Length","Time","Mass","Temp","Charge"
  
  ! test lowercase retrieval of invariant
  call PhysicalConstants_get("pi",value)
  write(fileUnit,920)"pi",value
  if (value .NE. PI) perfect = .FALSE.

  ! test uppercase retrieval of invariant
  call PhysicalConstants_get("PI",value)
  write(fileUnit,920)"PI",value
  if (value .NE. PI) perfect = .FALSE.

  ! test retrieval of unknown constant
  call PhysicalConstants_get("Unknown",value)
  write(fileUnit,920)"Unknown",value

  ! test retrieval of invariant with different units
  call PhysicalConstants_get("Pi", value,unitLength="m",unitCharge="C")
  write(fileUnit,921)"Pi",value,"m","-","-","-","C"
  if (value .NE. PI) perfect = .FALSE.

  ! test retrieval of invariant with wrong  units
  !! WORKS!  as does Driver_abortFlash.
  !        call PhysicalConstants_get("Pi", value,unitLength="m",temp_unit="C")
  !        write(fileUnit,921)"Pi",value,"m","-","-","C","-"
  ! test retrieval of variable constant in default units
  call PhysicalConstants_get("Speed of light",value)
  write(fileUnit,920)"Speed of light",value

  ! test retrieval of variable constant in different (relevant) units
  call PhysicalConstants_get("Speed of light",value,                      &
       &           unitLength="m",unitTime="yr")
  write(fileUnit,921)"Speed of light",value,"m","yr","-","-","-"

  ! test retrieval of variable constant in different (irrelevant) units
  call PhysicalConstants_get("Speed of light",value,unitMass="kg",       &
       &               unitLength="m",unitTime="yr")
  write(fileUnit,921)"Speed of light",value,"m","yr","kg","-","-"

  !  for comparison with MKS below
  call PhysicalConstants_get("Newton",value)
  write(fileUnit,921)"Newton",value,"-","-","-","-","-"
  
  
  ! set default units to MKS instead of CGS
  !LBR removed to eliminate cross-unit calls.  But demonstrates a hack to 
  !LBR  allow on-the-fly changing of runtime parameters.
  !        call RuntimeParameters_read('flash_mks.par')
  !        call RuntimeParameters_get('pc_unitsBase',cgsORmks)
  !        write(fileUnit,'("Using base units ",A3)')cgsORmks

!!!!!!!!
  call pc_checkCGSMKS("mKS",isError)
  if (isError /= 0) then
     write(errorString,980)"mKS"
     call Driver_abortFlash(errorString)
     return
  endif
!!!!!!!!
  write(fileUnit,901)pc_nameSISystem
  write(fileUnit,902)pc_SISystem


  ! test retrieval of variable in default MKS units
  call PhysicalConstants_get("Newton",value)
  write(fileUnit,921)"Newton",value,"-","-","-","-","-"

  ! test retrieval in requested units that agree with MKS
  call PhysicalConstants_get("Newton",value,unitMass="kg",unitLength="m")
  write(fileUnit,921)"Newton",value,"m","-","kg","-","-"

  ! test retrieval in requested units conflicting with MKS
  call PhysicalConstants_get("Newton",value,unitMass="kg",unitLength="cm")
  write(fileUnit,921)"Newton",value,"cm","-","kg","-","-"

  ! test retrieval in irrelevant units
  call PhysicalConstants_get("Newton",value,unitCharge="esu")
  write(fileUnit,921)"Newton",value,"-","-","-","esu","-"

  ! test retrieval with irrelevant units
  call PhysicalConstants_get("Speed of light",value,unitMass="kg",       &
       &               unitTime="yr")
  write(fileUnit,921)"Speed of light",value,"-","yr","kg","-","-"

  ! test actual value in MKS units
  call PhysicalConstants_get("ideal gas constant",value)
  write(fileUnit,921)"ideal gas constant",value,"-","-","-","-","-"
  if (abs(value - 8.3145).GT.0.0005) then
     write(fileUnit,*) "Bad value!"
     perfect = .FALSE.
  end if
  if (pc_globalMe.EQ.MASTER_PE) print*,'Done with calling PhysicalConstants_get'
  
  return 
  
  !------------------------------------------------------------------------
  
900 format("Testing PhysicalConstants")
901 format("Defined base units ",'(',2('"',A3,'"',1x),')')
902 format("Using base unit system #",I2)
904 format("---------------List all Units ------------",/,              &
         &              T15,"Unit",T45,"CGS Value",T70,"Base Unit")
906 format("---------------List all Units---------------")

910 format("Now testing retrieving constants with different units",/,       &
         &       T5,"Constant Name",T25,"Constant Value",T35,5(4X,A6))
920 format(A15,ES15.5,T35,5(9x,"-"))
921 format(A15,ES15.5,T35,5(5x,A5))
922 format("-----------List all CGS Constants----------------")  
980 format("PhysicalConstants_init: pc_unitsBase ", A6, &
       &         " invalid, must be CGS/MKS")

  
end subroutine PhysicalConstants_unitTest
