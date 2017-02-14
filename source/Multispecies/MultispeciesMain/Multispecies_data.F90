!!****if* source/Multispecies/MultispeciesMain/Multispecies_data
!!
!! NAME
!!
!!  Multispecies_data
!!
!! SYNOPSIS
!!
!!  use Multispecies_data
!!
!! DESCRIPTION
!!
!!  Holds data variables for the Multispecies unit.
!!  When used in routines, variable names begin with ms_
!!
!! ARGUMENTS
!!
!!
!!***

module Multispecies_data

  real,parameter      :: ms_tinyNumber=1.0e-12

#include "Flash.h"

  ! Every unit should know these.....
  integer, save       :: ms_globalMe
  
#include "Multispecies.h"

  type multispecies_type
     integer :: name !integer because names are defined in .h file
     real :: numTotal
     real :: numPositive
     real :: numNeutral
     real :: numNegative
     real :: bindingEnergy
     real :: adiabaticIndex
     real :: zMin
     real :: opacityLowTemp
     integer :: numElems
     integer :: zElems(MS_MAXELEMS)
     real :: aElems(MS_MAXELEMS)
     real :: fractions(MS_MAXELEMS)

     integer :: eosType
     integer :: eosSubtype

     character(len=MS_STRINGLEN) :: eosZfreeTableFile
     character(len=MS_STRINGLEN) :: eosEnerTableFile
     character(len=MS_STRINGLEN) :: eosPresTableFile
     character(len=MS_STRINGLEN) :: eosGroupName
  end type multispecies_type


  type(multispecies_type), save :: ms_Array(NSPECIES)

CONTAINS
    function ms_isZero(smallReal) RESULT (yesno)
        implicit none
        real, intent(in)        ::  smallReal
        logical                 ::  yesno
        
        yesno = .false.
        if (abs(smallReal) .lt. ms_tinyNumber) yesno = .true.
    end function ms_isZero
    

end module Multispecies_data
