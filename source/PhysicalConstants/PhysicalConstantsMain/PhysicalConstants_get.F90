!!****if* source/PhysicalConstants/PhysicalConstantsMain/PhysicalConstants_get
!!
!! NAME
!!  PhysicalConstants_get
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_get(character(:,IN) :: name,
!!            real(OUT) :: value,
!!            real(IN,OPTIONAL) :: unitLength,
!!            real(IN,OPTIONAL) :: unitTime,
!!            real(IN,OPTIONAL) :: unitMass,
!!            real(IN,OPTIONAL) :: unitTemp,
!!            real(IN,OPTIONAL) :: unitCharge,
!!            real(IN,OPTIONAL) :: unitSubstAmount)
!!
!! DESCRIPTION
!!    Retrieve the value of a physical constant.  Returns in CGS or MKS units, 
!!      depending upon status of runtime parameter pc_unitsBase.
!!    If called with optional unit arguments, the basic constant is converted and 
!!      returned with the requested units.
!!    If the constant name is not recognized, or some other error, 
!!      returns value of 0
!!
!! ARGUMENTS
!!
!!  name:       character(:,IN)   title of the physical constant required
!!  value:      real(OUT)         returned value of the physical constant in appropriate units
!!  unitLength: real(IN,OPTIONAL) value will be returned with this requested unit for length
!!  unitMass:   real(IN,OPTIONAL) value will be returned with this requested unit for mass
!!  unitTime:   real(IN,OPTIONAL) value will be returned with this requested unit for time
!!  unitTemp:   real(IN,OPTIONAL) value will be returned with this requested unit for temperature
!!  unitCharge: real(IN,OPTIONAL) value will be returned with this requested unit for charge
!!  unitSubstAmount: real(IN,OPTIONAL) value will be returned with this requested unit for
!!                                     substance amount
!!
!! EXAMPLE
!!    Example of usage:
!!    use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
!!    ............
!!    real ::   pi, newton
!!    ............
!!    call PhysicalConstants_get("Pi",pi)     ! deprecated, use PI from constants.h instead
!!    call PhysicalConstants_get("newton",newton,unitMass="kg",unitLength="cm")
!!
!! NOTES
!!   
!!   Because PhysicalConstants_get has optional arguments in the interface,
!!   the interface file PhysicalConstants_interface must be used in the calling routine.
!!  
!!   If you see a constant named "Newton", what is meant is probably not the SI unit
!!   of force which has dimensions LENGTH^1 * TIME^(-2) * MASS^1, but 
!!   Newton's gravitational constant G.
!!
!!***

!************************************************************************

        subroutine PhysicalConstants_get (name, value,                          &
     &            unitLength, unitTime, unitMass, unitTemp, unitCharge, unitSubstAmount)

!------------------------------------------------------------------------
#include "constants.h"
        use PhysicalConstants_data, ONLY:   pc_arrayUnit, pc_arrayConstant,     &
     &          pc_typeConstant, pc_initialized, pc_SISystem,                   &
     &          pc_baseUnitLength, pc_baseUnitTime, pc_baseUnitTemp,            &
     &          pc_baseUnitCharge, pc_baseUnitMass, pc_baseUnitSubstAmount
        use PhysicalConstants_interface, ONLY : PhysicalConstants_init

        implicit none
        
 !------------------------------------------------------------------------ 
        character(len=*), intent(in)                :: name
        real, intent(out)                           :: value
        character(len=*), optional, intent(in)      :: unitLength, unitTime,    & 
     &                                                 unitMass, unitTemp,      &
     &                                                 unitCharge, unitSubstAmount

        type (pc_typeConstant)          :: constant
        real                            :: scaleLength, scaleTime, & 
     &                                     scaleMass, scaleTemp, & 
     &                                     scaleCharge, scaleSubstAmount
        integer                            :: indexConstant, indexUnit, MyPE, numProcs

!------------------------------------------------------------------------

! initializations
        if (.not. pc_initialized) then
            call PhysicalConstants_init()
        endif
        value = 0.0
 
!  Find the constant in the database.  If it's not found,
!       exit with a value of zero.

        call pc_findConstant (name, indexConstant)
        if (indexConstant == 0) then
            return
        else
            constant = pc_arrayConstant(indexConstant)
        endif

!  Establish the unit scaling.  If a unit is not found (or is illegally 
!   specified, e.g., unitTime="Mpc"), exit with a value of zero.
!   If MKS requested (pc_SISystem==2) then return in MKS units

        if (present(unitLength)) then
          call pc_findUnit (unitLength, pc_baseUnitLength, indexUnit)
          if (indexUnit == 0) return
          scaleLength = pc_arrayUnit(indexUnit)%cgsValue
        else
          scaleLength = 1.e0
          if (pc_SISystem == 2) scaleLength = 1.E2  ! meters
        endif

        if (present(unitTime)) then
          call pc_findUnit (unitTime, pc_baseUnitTime, indexUnit)
          if (indexUnit == 0) return
          scaleTime = pc_arrayUnit(indexUnit)%cgsValue
        else
          scaleTime = 1.e0
          if (pc_SISystem == 2) scaleTime = 1.0  ! seconds
        endif

        if (present(unitMass)) then
          call pc_findUnit (unitMass, pc_baseUnitMass, indexUnit)
          if (indexUnit == 0) return
          scaleMass = pc_arrayUnit(indexUnit)%cgsValue
        else
          scaleMass = 1.e0
          if (pc_SISystem == 2) scaleMass = 1.E3  ! kilograms
        endif

        if (present(unitTemp)) then
          call pc_findUnit (unitTemp, pc_baseUnitTemp, indexUnit)
          if (indexUnit == 0) return
          scaleTemp = pc_arrayUnit(indexUnit)%cgsValue
        else
          scaleTemp = 1.e0
          if (pc_SISystem == 2) scaleTemp = 1.e0  ! Kelvin
        endif

        if (present(unitCharge)) then
           call pc_findUnit (unitCharge, pc_baseUnitCharge, indexUnit)
           if (indexUnit == 0) return
           scaleCharge = pc_arrayUnit(indexUnit)%cgsValue
        else
          scaleCharge = 1.e0
           if (pc_SISystem == 2) scaleCharge = 2.99792458E9  ! Coulombs
        endif
   
        if (present(unitSubstAmount)) then
          call pc_findUnit (unitSubstAmount, pc_baseUnitSubstAmount, indexUnit)
          if (indexUnit == 0) return
          scaleSubstAmount = pc_arrayUnit(indexUnit)%cgsValue
        else
          scaleSubstAmount = 1.e0
          if (pc_SISystem == 2) scaleSubstAmount = 1.e0  ! mole
        endif


 
!    Convert the constant to the requested units and return it.

        value = constant%cgsValue / & 
     &          ( (scaleLength ** constant%unitExponent(pc_baseUnitLength)) *   & 
     &            (scaleTime ** constant%unitExponent(pc_baseUnitTime)) *       & 
     &            (scaleMass ** constant%unitExponent(pc_baseUnitMass)) *       & 
     &            (scaleTemp ** constant%unitExponent(pc_baseUnitTemp)) *       & 
     &            (scaleCharge ** constant%unitExponent(pc_baseUnitCharge)) *       & 
     &            (scaleSubstAmount ** constant%unitExponent(pc_baseUnitSubstAmount)) )

!------------------------------------------------------------------------

        return
        end subroutine PhysicalConstants_get

