!!****h* source/PhysicalConstants/PhysicalConstants_interface
!!
!! This is the header file for the PhysicalConstants module that defines its
!! public interfaces.
!!
!!***
Module PhysicalConstants_interface

#include "constants.h"
#include "Flash.h"

  interface PhysicalConstants_get
    subroutine PhysicalConstants_get (name, value,                          &
                    &            unitLength, unitTime, unitMass, unitTemp, unitCharge, unitSubstAmount)
      character(len=*), intent(in)                :: name
      real, intent(out)                           :: value
      character(len=*), optional, intent(in)      :: unitLength, unitTime,    & 
                    &                                unitMass, unitTemp,      & 
                    &                                unitCharge, unitSubstAmount
    end subroutine PhysicalConstants_get
  end interface

  interface PhysicalConstants_init
    subroutine PhysicalConstants_init()       
    end subroutine PhysicalConstants_init
  end interface

  interface PhysicalConstants_list
    subroutine PhysicalConstants_list(fileUnit)
      integer, intent(in)                :: fileUnit
    end subroutine PhysicalConstants_list
  end interface

  interface PhysicalConstants_listUnits
    subroutine PhysicalConstants_listUnits (fileUnit)
      integer, intent(in)                 :: fileUnit
    end subroutine PhysicalConstants_listUnits
  end interface

  interface PhysicalConstants_unitTest
    subroutine PhysicalConstants_unitTest(fileUnit,perfect)
      integer, intent(in)         :: fileUnit
      logical, intent(out)        :: perfect
    end subroutine PhysicalConstants_unitTest
  end interface

end Module PhysicalConstants_interface
