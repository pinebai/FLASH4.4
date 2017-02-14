!!****f* source/PhysicalConstants/PhysicalConstants_init
!!
!! NAME
!!  PhysicalConstants_init
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_init()
!!
!! DESCRIPTION
!!
!! This subroutine initializes the Physical Constants databases for
!! units and constants.  Must be called in the Simulation_init
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  pc_unitsBase   [default "CGS"] set the default system of units, either "CGS" or "MKS"
!!                 CGS:  centimeters, grams, seconds, charge=esu
!!                 MKS:  meters,  kilogramks, seconds, charge=Coloumb
!!                     both systems have temperature in Kelvin
!!
!! NOTES
!!
!!***            

subroutine PhysicalConstants_init()
  implicit none 
  return
end subroutine PhysicalConstants_init

