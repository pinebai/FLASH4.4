!!****if* source/physics/sourceTerms/Ionize/IonizeMain/ion_readTable
!!
!! NAME
!!
!!  ion_readTable
!!
!! SYNOPSIS
!!
!!  ion_readTable()
!!  
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION 
!!  
!!  Reads the table containing ionization and recombination 
!!  coefficients of Summers.
!!
!! SIDE EFFECTS
!!
!!  ion_symbols - List of elements available in Summers data file. 
!!  ion_nelect - Number of electrons for each available element.
!!  ion_nion12 - Number of ionization states.
!!  ion_cfinz - Ionization coefficients.
!!  ion_cfric - Recombination coefficients.
!!  ion_tp - Temperatures at the evaluated coefficients.
!!
!!  MODIFICATION HISTORY:
!!          WRITTEN BY: S. ORLANDO, AUGUST 2001
!!
!!***

subroutine ion_readTable()
#include "Ionize.h"
  use ion_interface, ONLY : ion_selectElements
  use Ionize_data, ONLY : ion_nion12, ion_tp, &
       ion_cfinz, ion_cfric, ion_nelect, ion_symbols
  implicit none

  !..declare
  integer :: nel,i,j


  !Additions by Chris for FLASH3 - 14 July 2008.
  !-----------------------------------------------------------------
  !In this implementation of ion_readTable, we have data for the 
  !following elements stored in file in this order:
  ion_symbols(1)='he'; ion_symbols(2)='c '; ion_symbols(3)='n '
  ion_symbols(4)='o '; ion_symbols(5)='ne'; ion_symbols(6)='mg'
  ion_symbols(7)='si'; ion_symbols(8)='s '; ion_symbols(9)='ar'
  ion_symbols(10)='ca'; ion_symbols(11)='fe'; ion_symbols(12)='ni'
 
  ion_nelect(1) = 2; ion_nelect(2) = 6; ion_nelect(3) = 7; 
  ion_nelect(4) = 8; ion_nelect(5) = 10; ion_nelect(6) = 12; 
  ion_nelect(7) = 14; ion_nelect(8) = 16; ion_nelect(9) = 18; 
  ion_nelect(10) = 20; ion_nelect(11) = 26; ion_nelect(12) = 28; 

  do nel = 1,ION_NELEM
     ion_nion12(nel) = ion_nelect(nel)+1
  enddo


  !We need to check that the desired simulation can be
  !evolved in the Ionize unit.  We are limited to the elements that are  
  !stored in summers_den_1e8.rates table.  
  ! *** Abort if we have requested a fluid that has no ionize specific data. ***
  !If we have data for each species in the simulation, initialise 
  !several arrays of use to ionize unit: ion_xfrac, ion_idx.
  call ion_selectElements()
  !-----------------------------------------------------------------


  !Read data for the 12 elements.  Later routines will reference 
  !ion_idx to find out whether the elements are requested by the 
  !simulation.

  !DEV:  When we perform an interpolation in intcoeff we touch 
  !uninitialised data.  This also happens in FLASH2.  As I am unsure 
  !what intcoeff is actually doing, I will initialise our arrays to 0.0.
  !I figure interpolating with a zero value is better than with a 
  !random real value!
  ion_cfinz(:,:,:) = 0.0
  ion_cfric(:,:,:) = 0.0

  open(8,file='summers_den_1e8.rates',status= 'unknown')
  do i = 1,ION_NTEMP
     read(8,*) ion_tp(i)
     do nel = 1,ION_NELEM
        read(8,*)(ion_cfinz(nel,i,j),j=1,ion_nion12(nel))
        read(8,*)(ion_cfric(nel,i,j),j=1,ion_nion12(nel))
     end do
  end do
  close (8)


  return
end subroutine ion_readTable

