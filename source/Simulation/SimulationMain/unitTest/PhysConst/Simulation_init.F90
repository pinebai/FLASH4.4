!!****if* source/Simulation/SimulationMain/unitTest/PhysConst/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for the unitTest of Physical Constants
!!    by calling PhysicalConstants_init()
!!
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!  pc_unitsBase   [default "CGS"] set the default system of units, either "CGS" or "MKS"
!!                 CGS:  centimeters, grams, seconds, charge=esu
!!                 MKS:  meters, kilometers, seconds, charge=Coloumb
!!                     both systems have temperature in Kelvin
!!
!!***

subroutine Simulation_init()
    use PhysicalConstants_interface, ONLY : PhysicalConstants_init
  
    implicit none

#include "constants.h"
  
   

   ! Physical Constants is initialized in Driver_initFlash

end subroutine Simulation_init







