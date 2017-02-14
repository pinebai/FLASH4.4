!!****f* source/PhysicalConstants/PhysicalConstants_listUnits
!!
!! NAME
!!  PhysicalConstants_listUnits
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_listUnits(integer(in) :: fileUnit)            
!!
!! DESCRIPTION
!!
!!  Writes the units of measurement to standard output
!! 
!! ARGUMENTS
!!
!!     fileUnit - integer(in)  file number to dump information
!!
!! PARAMETERS
!!
!!  pc_unitsBase   [default "CGS"] set the default system of units, either "CGS" or "MKS"
!!                 CGS:  centimeters, grams, seconds, charge=esu
!!                 MKS:  meters,  kilograms, seconds, charge=Coloumb
!!                     both systems have temperature in Kelvin
!!
!! NOTES
!!
!!***            

subroutine PhysicalConstants_listUnits (fileUnit)

  implicit none
  
  integer, intent(in)                 :: fileUnit

end subroutine PhysicalConstants_listUnits
