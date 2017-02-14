!!****f* source/PhysicalConstants/PhysicalConstants_list
!!
!! NAME
!!  PhysicalConstants_list
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_list(integer(in) :: fileUnit)            
!!
!! DESCRIPTION
!!
!!  Writes the physical constants to standard output
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
  

subroutine PhysicalConstants_list(fileUnit)
  
  implicit none

  integer, intent(in)                :: fileUnit

end subroutine PhysicalConstants_list

