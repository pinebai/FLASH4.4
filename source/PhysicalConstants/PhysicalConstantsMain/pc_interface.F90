!!****if* source/PhysicalConstants/PhysicalConstantsMain/pc_interface
!!
!! NAME
!!  pc_interface
!!
!! SYNOPSIS
!!
!!  use pc_interface, ONLY: pc_makeLowercase,pc_checkCGSMKS,pc_addConstant,pc_addUnit,pc_findConstant,pc_findUnit
!!
!! DESCRIPTION
!!
!! Interface module for internal subroutines to do miscellaneous things in the
!! PhysicalConstants unit.
!!
!! NOTE
!!
!!  This interface should be considered private to the RuntimeParameters unit.
!!
!!***

Module pc_interface
#include "constants.h"

  interface
!*****************************************************************************
!  Routine:    pc_makeLowercase()
 
!  Description: Convert a string to lowercase.
 
     SUBROUTINE pc_makeLowercase (str, strLower)
       implicit none
 !------------------------------------------------------------------------
       character(len=*), intent(in)                :: str
       character(len=*), intent(out)               :: strLower
     END SUBROUTINE pc_makeLowercase
  end interface

!************************************************************************

  interface
!  Routine:  pc_checkCGSMKS()

!  Description:  checks validity of cgsORmks input.
!  Returns isError = -1 if input invalid
     SUBROUTINE pc_checkCGSMKS(cgsORmks,isError)
!------------------------------------------------------------------------
       implicit none
!------------------------------------------------------------------------
       character(len=3), intent(in)     :: cgsORmks
       integer, intent(out)             :: isError

       character(len=3)                 :: cgsORmksLower

!------------------------------------------------------------------------

     END SUBROUTINE pc_checkCGSMKS
  end interface

 !************************************************************************

  interface
!  Routine:     pc_addConstant ()

!  Description: Add a physical constant to the constants database.  You
!               must specify the name of the constant (a string), its
!               value in CGS units, and the exponents for its scaling with
!               length, time, mass, temperature, and charge.


     SUBROUTINE pc_addConstant (name, cgsValue, exponentLength,              & 
     &                                 exponentTime, exponentMass,              & 
     &                                 exponentTemp, exponentCharge,            &
                                       exponentSubstAmount)

       implicit none
!------------------------------------------------------------------------
       character(len=*), intent(in)      :: name
       real, intent(in)                  :: cgsValue, exponentLength,          & 
     &                                       exponentTime, exponentMass,        & 
     &                                       exponentTemp, exponentCharge
       real, intent(in), OPTIONAL         ::               exponentSubstAmount
     
     END SUBROUTINE pc_addConstant
  end interface

!************************************************************************

  interface
!  Routine:     pc_addUnit ()

!  Description: Add a unit of measurement to the units database.  You
!               must specify the name of the unit (a string), its base
!               unit type ("length", "time", "mass", "temperature", or
!               "charge"), and its value in CGS units.


     subroutine pc_addUnit (baseUnitKind, name, cgsValue)
       implicit none
!------------------------------------------------------------------------
       character(len=*), intent(in)        :: name, baseUnitKind
       real, intent(in)                    :: cgsValue
     END SUBROUTINE pc_addUnit
  end interface

 !************************************************************************

  interface
!  Routine:     pc_findConstant ()

!  Description: Find a constant in the array using its name.
!               Return the index of the constant's node if found,
!               otherwise return zero.

     
     SUBROUTINE pc_findConstant (name, index)

       implicit none
!------------------------------------------------------------------------    
       character(len=*), intent(in)        :: name
       integer, intent(out)                :: index
     END SUBROUTINE pc_findConstant
  end interface

!************************************************************************

  interface
!  Routine:     pc_findUnit ()

!  Description: Find a unit in a linked list using its name.
!               Return the index to the unit's node if it's found,
!               otherwise return zero.  In looking for the
!               unit, we also check to make sure that the base unit
!               type matches.
     
     
     SUBROUTINE pc_findUnit (name, baseUnit, index)

       implicit none
!------------------------------------------------------------------------
       
       character(len=*), intent(in)               :: name
       integer, intent(in)                        :: baseUnit
       integer, intent(out)                       :: index
     END SUBROUTINE pc_findUnit
  end interface

end Module pc_interface
