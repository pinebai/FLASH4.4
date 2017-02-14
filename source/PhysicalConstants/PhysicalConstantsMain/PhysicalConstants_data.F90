!!****if* source/PhysicalConstants/PhysicalConstantsMain/PhysicalConstants_data
!!
!! NAME
!!  PhysicalConstants_data
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_data()
!!
!! DESCRIPTION
!!
!! This module holds the data structures for the Physical Constants unit
!!
!! ARGUMENTS
!!
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
!!***!            

MODULE PhysicalConstants_data

#include "constants.h"

!-------------------------------------------------------------------------
  ! processor numbers
  integer, save :: pc_globalMe, pc_globalNumProcs, pr_globalComm

  !  Base units:  length, time, mass, temperature, and charge
  integer, parameter :: PC_NBASEUNITS = 6
  integer, parameter :: pc_baseUnitLength = 1,  pc_baseUnitTime = 2, & 
       &                        pc_baseUnitMass = 3, pc_baseUnitTemp = 4, & 
       &                        pc_baseUnitCharge = 5, pc_baseUnitSubstAmount = 6

  !  Using CGS or MKS units?
  integer, parameter :: PC_NSISYSTEMS = 2
  character(len=3), dimension(PC_NSISYSTEMS), save ::                     &
       &         pc_nameSISystem = (/ 'cgs','mks' /)
  integer, save ::    pc_SISystem
  !  Has PhysicalConstansts been initialized?
  logical, save ::    pc_initialized = .false.        
  !  CGS units:   cm=meter, s=second, g=gram, K=kelvin, esu=electrostatic unit, mol
  !  MKS units:  m=meter; s=second, kg=kilogram, K=kelvin, C=coulomb,           mol
  ! NOTE -- an "esu" is part of the electrostatic CGS system.  It is sometimes
  ! called a "unit of Franklin."  1 esu = 3.3356E-10 Coulombs 
  ! NOTE -- need to use RESHAPE for some compilers for 2d initialization                                              
  character(len=3), dimension(PC_NSISYSTEMS,PC_NBASEUNITS), save ::       &
       &      pc_nameUnitsBase =  RESHAPE(                                        &
       &      (/ 'cm ', 'm  ',               & ! length - centimeter / meter
       &         's  ', 's  ',               & ! time - second / second                             
       &         'g  ', 'kg ',               & ! mass - gram / kilogram                                   
       &         'K  ', 'K  ',               & ! temperature - Kelvin / Kelvin
       &         'esu', 'C  ',               & ! charge:unit of Franklin/Coulomb
       &         'mol', 'mol'               /), & ! amount of substance
       &          (/ 2,6 /))

  !   Derived type declarations for physical constants.

  !  A constant is defined by its name, value in CGS units, and the
  !   exponent of the unit.  e.g. the speed of light is 2.99792458E10 cm / sec
  !   and is defined in this format as name="speed of light", cgsValue="2.99E10"
  !   unitExponent(1=length)=1, unitExponent(2=time)=-1, unitExponent(3=mass)=0, 
  !   unitExponent(4=temperature)=0, unitExponent(5=charge)=0
  type pc_typeConstant
     character(len=MAX_STRING_LENGTH)    :: name
     real                                :: cgsValue
     real                                :: unitExponent(PC_NBASEUNITS)
  end type pc_typeConstant

  !  A unit is defined by its abbreviated name, base unit (must be "length",
  !   "time","mass","temperature" or "charge"), and the value.  e.g. a Coulomb is 
  !   a unit of charge equal to 2.99E9 esu. and is represented here by 
  !   name="C",baseUnit="charge", and cgsValue="2.99E9"
  type pc_typeUnit
     character(len=MAX_STRING_LENGTH)    :: name
     integer                             :: baseUnit
     real                                :: cgsValue
  end type pc_typeUnit

  !! Number of entries in the PhysicalConstants unit and constant array 
  integer, parameter :: PC_NUNITS = 29
  integer, parameter :: PC_NCONSTANTS = 16
  TYPE(pc_typeConstant), dimension(1:PC_NCONSTANTS)    :: pc_arrayConstant
  TYPE(pc_typeUnit), dimension(1:PC_NUNITS)            :: pc_arrayUnit
  ! actual number of array units used
  integer, save                           :: pc_sizeConstant
  integer, save                           :: pc_sizeUnit


END MODULE PhysicalConstants_data

