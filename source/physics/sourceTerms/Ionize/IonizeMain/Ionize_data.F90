!!****if* source/physics/sourceTerms/Ionize/IonizeMain/Ionize_data
!!
!! NAME
!!
!!  Ionize_data
!!
!! SYNOPSIS
!!
!!  use Ionize_data
!!  
!! DESCRIPTION
!!
!!  Data module for Ionize.
!!   
!!***

Module Ionize_data
#include "Flash.h"
#include "Ionize.h"
  implicit none
  integer, parameter :: ion_ELEC = ELEC_SPEC - SPECIES_BEGIN + 1
  integer, parameter :: ion_HYD = H_SPEC - SPECIES_BEGIN + 1

!..ion_nion   = number of ionization species for the element considered
!..ion_den    = electron number density
!..ion_cz     = Ionization coefficients
!..ion_al     = Recombination coefficients
  integer, dimension(ION_NELEM), save :: ion_idx, ion_nion12
  real, dimension(NSPECIES), save :: ion_xfrac
  real, save :: ion_tneimin, ion_tneimax, ion_dneimin, ion_dneimax
  real, save :: ion_emass
  real, save :: ion_smallx

!..recombination coefficients derived from Summers.
!
!
  real,save,dimension(ION_NTEMP) :: ion_tp
  real,save,dimension(ION_NELEM,ION_NTEMP,0:ION_NIMAX)::ion_cfinz
  real,save,dimension(ION_NELEM,ION_NTEMP,ION_NIMAX):: ion_cfric
  real,save ::  ion_btemp,ion_den
  integer,save :: ion_nion
  

  real,dimension(0:ION_NIMAX):: ion_cz
  real,save,dimension(ION_NIMAX) :: ion_al
  integer,save :: ion_ifirst
  logical, save :: ion_useIonize

  !This array defines the elements available in the summers_den_1e8.rates table,
  !and the order in which they appear in the table.
  character (len=2), save, dimension(ION_NELEM) :: ion_symbols
  integer, save, dimension(ION_NELEM) :: ion_nelect

end Module Ionize_data
