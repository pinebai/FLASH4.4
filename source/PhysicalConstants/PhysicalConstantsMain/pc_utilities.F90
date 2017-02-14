!!****if* source/PhysicalConstants/PhysicalConstantsMain/pc_utilities
!!
!! NAME
!!  pc_utilities
!!
!! SYNOPSIS
!!
!!  pc_makeLowercase()
!!  pc_checkCGSMKS()
!!  pc_addConstant()
!!  pc_addUnit()
!!  pc_findConstant()
!!  pc_findUnit()
!!
!! DESCRIPTION
!!
!! Internal subroutines to do miscellaneous things for the PhysicalConstants
!!      modules
!!
!! ARGUMENTS
!!
!! None
!!
!! NOTE
!!
!!  These routines should not be called by the user
!!
!!***

!*****************************************************************************
!  Routine:    pc_makeLowercase()
 
!  Description: Convert a string to lowercase.
 
        SUBROUTINE pc_makeLowercase (str, strLower)
        implicit none
 !------------------------------------------------------------------------
        character(len=*), intent(in)                :: str
        character(len=*), intent(out)               :: strLower
        integer                                     :: i
 
 !------------------------------------------------------------------------
 
        strLower(1:len(str)) = str 
        do i = 1, len_trim(str)
          if (lge(str(i:i), 'A') .and. lle(str(i:i), 'Z'))                      & 
     &      strLower(i:i) = achar( iachar(str(i:i)) + 32 )
        enddo
 
 !------------------------------------------------------------------------
 
        return 
        END SUBROUTINE pc_makeLowercase

!************************************************************************

!  Routine:  pc_checkCGSMKS()

!  Description:  checks validity of cgsORmks input.
!  Returns isError = -1 if input invalid
        SUBROUTINE pc_checkCGSMKS(cgsORmks,isError)
!------------------------------------------------------------------------
        use PhysicalConstants_data, ONLY : pc_nameSISystem, pc_SISystem
        implicit none
!------------------------------------------------------------------------
        character(len=3), intent(in)     :: cgsORmks
        integer, intent(out)             :: isError

        character(len=3)                 :: cgsORmksLower

!------------------------------------------------------------------------

! initializations        
        isError = 0
        call pc_makeLowercase(cgsORmks,cgsORmksLower)
! Check against allowed values        
        if (cgsORmksLower .eq. pc_nameSISystem(1)) then
            pc_SISystem = 1
        elseif (cgsORmksLower .eq. pc_nameSISystem(2)) then
            pc_SISystem = 2
        else
            pc_SISystem = 0
            isError = -1
        endif

!------------------------------------------------------------------------

        return
        END SUBROUTINE pc_checkCGSMKS

 
!************************************************************************

!  Routine:     pc_addConstant ()

!  Description: Add a physical constant to the constants database.  You
!               must specify the name of the constant (a string), its
!               value in CGS units, and the exponents for its scaling with
!               length, time, mass, temperature, and charge.


        SUBROUTINE pc_addConstant (name, cgsValue, exponentLength,              & 
     &                                 exponentTime, exponentMass,              & 
     &                                 exponentTemp, exponentCharge,            &
                                       exponentSubstAmount)

!------------------------------------------------------------------------
        use PhysicalConstants_data, ONLY  : pc_typeConstant, pc_sizeConstant,   &
     &                                       pc_arrayConstant,                  &
     &          pc_baseUnitLength, pc_baseUnitTime, pc_baseUnitTemp,                  &
     &          pc_baseUnitCharge, pc_baseUnitMass, pc_baseUnitSubstAmount
        implicit none
!------------------------------------------------------------------------
        character(len=*), intent(in)      :: name
        real, intent(in)                  :: cgsValue, exponentLength,          & 
     &                                       exponentTime, exponentMass,        & 
     &                                       exponentTemp, exponentCharge
        real, intent(in),OPTIONAL         ::               exponentSubstAmount
     
        type (pc_typeConstant)            :: node
        integer                           :: istat
        character(len=len(name))          :: nameLower
        real                              :: value

!------------------------------------------------------------------------

! initializations
        value = 0.0

!     Check to make sure the constant isn't already in the  database.
        call pc_findConstant (name, istat)
       
        if (istat /= 0) then       
            write(*,*)' PhysicalConstants_addConstant: Constant ',name,         &
    &               ' already exists.'
            return
!     If not, add to array   
        else
            node%name                     = name
            node%cgsValue                = cgsValue
            node%unitExponent(pc_baseUnitLength)  = exponentLength
            node%unitExponent(pc_baseUnitTime) = exponentTime
            node%unitExponent(pc_baseUnitMass) = exponentMass
            node%unitExponent(pc_baseUnitTemp) = exponentTemp
            node%unitExponent(pc_baseUnitCharge)  = exponentCharge
            if(present(exponentSubstAmount)) then
               node%unitExponent(pc_baseUnitSubstAmount) = exponentSubstAmount
            else
               node%unitExponent(pc_baseUnitSubstAmount) = 0
            end if
            !! increase number of values
            pc_sizeConstant = pc_sizeConstant + 1
            !! add to end of array
            pc_arrayConstant(pc_sizeConstant) = node
        endif
 
!------------------------------------------------------------------------

        return
        END SUBROUTINE pc_addConstant

!************************************************************************

!  Routine:     pc_addUnit ()

!  Description: Add a unit of measurement to the units database.  You
!               must specify the name of the unit (a string), its base
!               unit type ("length", "time", "mass", "temperature", or
!               "charge"), and its value in CGS units.


        subroutine pc_addUnit (baseUnitKind, name, cgsValue)
!------------------------------------------------------------------------
        use PhysicalConstants_data, ONLY    : pc_typeUnit, pc_sizeUnit,         &
     &                                        pc_arrayUnit,                     &
     &          pc_baseUnitLength, pc_baseUnitTime, pc_baseUnitTemp,            &
     &          pc_baseUnitCharge, pc_baseUnitMass, pc_baseUnitSubstAmount
        implicit none
!------------------------------------------------------------------------
        character(len=*), intent(in)        :: name, baseUnitKind
        real, intent(in)                    :: cgsValue

        type (pc_typeUnit)                  :: node 
        integer                             :: istat, index
        integer                             :: baseUnit
        character(len=len(name))            :: nameLower
        character(len=len(baseUnitKind))    :: baseUnitKindLower

!------------------------------------------------------------------------

!      Is the base unit kind reasonable?

        baseUnit = 0
        call pc_makeLowercase(baseUnitKind, baseUnitKindLower)
        if (baseUnitKindLower == "length")      baseUnit = pc_baseUnitLength
        if (baseUnitKindLower == "time")        baseUnit = pc_baseUnitTime
        if (baseUnitKindLower == "mass")        baseUnit = pc_baseUnitMass
        if (baseUnitKindLower == "temperature") baseUnit = pc_baseUnitTemp
        if (baseUnitKindLower == "charge")      baseUnit = pc_baseUnitCharge
        if (baseUnitKindLower == "substance amount") baseUnit = pc_baseUnitSubstAmount

        if (baseUnit == 0) then
          write (*,*) 'PhysicalConstants_addUnit:  base unit ',                 &
     &              baseUnitKind, ' unknown'
          return
        endif

!       Check to make sure the unit isn't already in the database.
        call pc_findUnit (name, baseUnit, istat)
        if (istat /= 0) then
          write (*,*) 'PhysicalConstants_addUnit: Unit ', name,                 &
     &             ' already exists.'
          return
        else
!        If it isn't, create a node and add it to the unit database.
            node%name                    = name
            node%cgsValue                = cgsValue
            node%baseUnit                = baseUnit

            pc_sizeUnit = pc_sizeUnit + 1
            pc_arrayUnit(pc_sizeUnit) = node
        endif

!------------------------------------------------------------------------

        return
        END SUBROUTINE pc_addUnit

 !************************************************************************

!  Routine:     pc_findConstant ()

!  Description: Find a constant in the array using its name.
!               Return the index of the constant's node if found,
!               otherwise return zero.


        SUBROUTINE pc_findConstant (name, index)

!------------------------------------------------------------------------
#include "constants.h"
        use PhysicalConstants_data, ONLY    : pc_typeConstant, pc_sizeConstant, &
     &                                         pc_arrayConstant
        implicit none
!------------------------------------------------------------------------    
        character(len=*), intent(in)        :: name
        integer, intent(out)                :: index
        
        type (pc_typeConstant)              :: node
        character(len=len(name))            :: nameLower
        character(len=MAX_STRING_LENGTH)    :: nameSearch
        integer                             :: i
        
!------------------------------------------------------------------------

! Initializations
        index = 0
        call pc_makeLowercase(name, nameLower)
! Loop through all constants and compare names       
        do i=1, pc_sizeConstant
            node = pc_arrayConstant(i)
            call pc_makeLowercase(node%name, nameSearch)
            if (trim(nameSearch) == nameLower) then
                index = i
                return
            endif
        enddo

!------------------------------------------------------------------------

        return
        END SUBROUTINE pc_findConstant

!************************************************************************

!  Routine:     pc_findUnit ()

!  Description: Find a unit in a linked list using its name.
!               Return the index to the unit's node if it's found,
!               otherwise return zero.  In looking for the
!               unit, we also check to make sure that the base unit
!               type matches.


        SUBROUTINE pc_findUnit (name, baseUnit, index)

!------------------------------------------------------------------------

        use PhysicalConstants_data, ONLY           :  pc_typeUnit, pc_sizeUnit, &
     &                        pc_arrayUnit, PC_NBASEUNITS
        implicit none
!------------------------------------------------------------------------

        character(len=*), intent(in)               :: name
        integer, intent(in)                        :: baseUnit
        integer, intent(out)                       :: index
        
        type (pc_typeUnit)                         :: node
        character(len=len(name))                   :: nameLower
        character(len=MAX_STRING_LENGTH)           :: nameSearch
        character(len=6), dimension(PC_NBASEUNITS) :: pc_nameUnitsType =        &
     &      (/"length","time  ","mass  ","temp  ","charge","amount"/)
        integer                                      :: i
 
!------------------------------------------------------------------------
 
!  Initializations
        index = 0
        call pc_makeLowercase(name,nameLower)
!  Loop through all units and compare names                 
        do i=1, pc_sizeUnit
            node = pc_arrayUnit(i)
            call pc_makeLowercase(node%name, nameSearch)
            nameSearch = trim(nameSearch)
            !! check name concordance
            if (nameSearch == nameLower) then
                !! now check baseUnit types
                if (node%baseUnit == baseUnit) then
                   index = i
                   return
                else
                    write(*,900)name, pc_nameUnitsType(baseUnit),              &
     &                        pc_nameUnitsType(node%baseUnit)
                    return ! with value zero
                endif
            endif
        enddo

!------------------------------------------------------------------------

        return
  800   format('name lower ',A20,' name search ',a20)
  900   format('ERROR Unit ',A20,' is NOT of type ',A10,'; should be of type',A10)
        END SUBROUTINE pc_findUnit
        
