!!****if* source/physics/sourceTerms/Ionize/IonizeMain/ion_selectElements
!!
!! NAME
!!  
!!  ion_selectElements
!!
!!
!! SYNOPSIS
!! 
!!  call ion_selectElements()
!!
!!  
!! DESCRIPTION
!!
!! Initialises the module level array ion_idx(:) which specifies which
!! elements are being used in this simulation.  Copies the relevant
!! data from sim_specXfrac(:) to ion_xfrac(:).  xfrac is a quantity which
!! is defined as: mass fraction / population fraction.
!!
!! This is here because the flexibility that is available to the 
!! user is limited to the elements for which we have data.
!!
!!***

subroutine ion_selectElements()

  use Simulation_speciesData, ONLY : sim_specElementSymbol, sim_specXfrac, &
       sim_specElement, sim_specNumElect, sim_specSelected
  use Ionize_data, ONLY : ion_idx, ion_xfrac, ion_symbols, &
        ion_nion12, ion_nelect
  use Driver_interface, ONLY : Driver_abortFlash

  ! *****************************************************************
  !           ion_idx = vector with the elements ID
  !
  !           ion_idx(1)    include Helium    if = 1
  !           ion_idx(2)    include Carbon    if = 1
  !           ion_idx(3)    include Nitrogen  if = 1
  !           ion_idx(4)    include Oxygen    if = 1
  !           ion_idx(5)    include Neon      if = 1
  !           ion_idx(6)    include Magnesium if = 1
  !           ion_idx(7)    include Silicon   if = 1
  !           ion_idx(8)    include Sulfur    if = 1
  !           ion_idx(9)    include Argon     if = 1
  !           ion_idx(10)   include Calcium   if = 1
  !           ion_idx(11)   include Iron      if = 1
  !           ion_idx(12)   include Nickel    if = 1
  !
  ! *****************************************************************

  implicit none

#include "Flash.h"
#include "Ionize.h"

  !Leave enough buffer space when we are copying around stored simulation symbols.
  character (len=100) :: simElement, lowerCaseSimElement
  integer :: iSim, iIon, copyFromIndex, copyToIndex, copySize, iEachIon
  logical :: foundElement

  ion_xfrac(:) = 0.0

  !It is assumed that every simulation runs with the elements hydrogen & electron.
  ion_xfrac(1:2) = sim_specXfrac(1:2)

  !ion_nion12, ion_cfinz, ion_cfric, ion_tp are intialised in ion_readTable().  However, 
  !the data is stored in these arrays in a set element order.  We specify which elements 
  !to consider by the values in the ion_idx(:) array.  It is possible that the 
  !simulation elements are in a different order, so the code below ensures we activate 
  !the correct element in ion_idx(:) array.
  !------------------------------------------------------------------

  ion_idx(:) = 0  !Use no elements until we find them being requested by the simulation.

  !We must check that we recognise each of the requested simulation elements:
  do iSim = lbound(sim_specElementSymbol,1), ubound(sim_specElementSymbol,1), 1

     if (sim_specSelected(iSim) .eqv. .true.) then

        !Make simulation element string lower case to be in same format as
        !the ion_symbols(:) data.
        simElement = sim_SpecElementSymbol(iSim)
        call pc_makeLowercase(simElement,lowerCaseSimElement)


        !Iterate over each of the known elements in summers_den_le8.rates table.
        iIon = 1; foundElement = .false.
        do while ( (iIon <= ION_NELEM) .and. (foundElement .eqv. .false.) )

           if ( trim(lowerCaseSimElement) == trim(ion_symbols(iIon)) ) then

#ifdef DEBUG_IONIZE
              print *,  "Matched simluation element: ", &
                   sim_SpecElementSymbol(iSim), " at position", iSim, "with stored element: ", &
                   ion_symbols(iIon), " at position", iIon
#endif

              !Check that the simulation specifies the correct number of electrons.
              if (ion_nelect(iIon) /= sim_specNumElect(iSim)) then
                 print *, "The element:", sim_SpecElementSymbol(iSim), "has:", &
                      sim_specNumElect(iSim), "electrons in the simulation, but:", &
                      ion_nelect(iIon), "hard-coded in the Ionize unit"
                 call Driver_abortFlash("Number of electrons mismatch")
              end if

              ion_idx(iIon) = 1  !This means we want to use this element.


              !We managed to identify the element.  Store it in Ionize data structures 
              !in a manner compatible with the solvers ported from FLASH2.           
              !Copy data into the fixed position ion_xfrac array.
              !----------------------------------------------------------------------
              copyFromIndex = sim_specElement(iSim) - SPECIES_BEGIN + 1

              copyToIndex = 3  !Include Electron + Hydrogen + 1
              do iEachIon = 1, iIon-1, 1
                 copyToIndex = copyToIndex + (ion_nion12(iEachIon) * ion_idx(iEachIon))
              end do

              copySize = ion_nelect(iIon)

#ifdef DEBUG_IONIZE
              print *, "Copy from:", copyFromIndex, " to:", copyToIndex, " for size:", copySize
#endif

              ion_xfrac( copyToIndex : (copyToIndex+copySize) ) = &
                   sim_specXfrac( copyFromIndex : (copyFromIndex+copySize) )
              !----------------------------------------------------------------------


              foundElement = .true.
           end if

           iIon = iIon + 1
        end do


        !Sorry user, we know nothing about the element in your simulation file.
        if (foundElement .eqv. .false.) then
           print *, "Can't find: ", lowerCaseSimElement
           call Driver_abortFlash("Element does not exist in summers_den_1e8.rates table")
        end if

     end if

  end do


  return

end subroutine ion_selectElements
