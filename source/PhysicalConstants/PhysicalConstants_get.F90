!!****f* source/PhysicalConstants/PhysicalConstants_get
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
!! PARAMETERS
!!
!!  pc_unitsBase   [default "CGS"] set the default system of units, either "CGS" or "MKS"
!!                 CGS:  centimeters, grams, seconds, charge=esu
!!                 MKS:  meters, kilometers, seconds, charge=Coloumb
!!                     both systems have temperature in Kelvin
!!
!! EXAMPLE
!!    Example of usage:
!!    use PhysicalConstants_interface,ONLY:  PhysicalConstants_get
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
!!   If you see a constant named "Newton", what is meant is not the SI unit
!!   of force which has dimensions LENGTH^1 * TIME^(-2) * MASS^1, but 
!!   Newton's gravitational constant G.
!!
!!***

!************************************************************************

        subroutine PhysicalConstants_get (name, value,                          &
     &            unitLength, unitTime, unitMass, unitTemp, unitCharge, unitSubstAmount)

!------------------------------------------------------------------------
        implicit none
        
!------------------------------------------------------------------------ 

        character(len=*), intent(in)                :: name
        real, intent(out)                           :: value
        character(len=*), optional, intent(in)      :: unitLength, unitTime,    & 
     &                                                 unitMass, unitTemp,      &
     &                                                 unitCharge, unitSubstAmount
!------------------------------------------------------------------------

        value = 0.0
        return
        end subroutine PhysicalConstants_get

