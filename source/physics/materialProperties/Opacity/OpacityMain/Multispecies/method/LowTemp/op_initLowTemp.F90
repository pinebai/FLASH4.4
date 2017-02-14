!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_initLowTemp
!!
!! NAME
!!
!!  op_initLowTemp
!!
!! SYNOPSIS
!!
!!  call op_initLowTemp ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the section of the low temperature opacity unit. Several arrays are
!!  allocated here and initialized with all the data. The number of atomic elements
!!  and their corresponding atomic numbers for the current opacity run must be known
!!  at this stage, otherwise the program will stop with a message.
!!
!!  The goal of this routine is to get and hold only the data corresponding to the
!!  atomic elements needed in memeory. The big arrays containing the data for all the
!!  elements will only be alive during the time spent in this routine.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initLowTemp ()

  use Driver_interface,            ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Opacity_data,                ONLY : op_totalSpecies,              &
                                          op_totalElements,             &
                                          op_nEnergyGroups,             &
                                          op_element2AtomicNumber,      &
                                          op_writeOpacityInfo
                               

  use op_lowTempData,              ONLY : op_initializedLowTemp,    &
                                          op_ignoreKleinNishina,    &
                                          op_maxElements,           &
                                          op_maxNstepsLowTemp,      &
                                          op_maxJmax,               &
                                          op_A1group,               &
                                          op_A2group,               &
                                          op_A3group,               &
                                          op_A4group,               &
                                          op_intLimits,             &
                                          op_Aij4,                  &
                                          op_Jmax,                  &
                                          op_PEenergyRange,         &
                                          op_elementAij4,           &
                                          op_elementJmax,           &
                                          op_elementPEenergyRange,  &
                                          op_tableLowTemp,          &
                                          op_PlanckLowTempTables,   &
                                          op_RosselandLowTempTables
                                  
  use op_interface,                ONLY : op_setPEcoeffsAij4,           &
                                          op_setPEarrayJmax,            &
                                          op_setPEenergyRange,          &
                                          op_calculateLowTempOpacities, &
                                          op_writeCompletePEdata,       &
                                          op_writeElementsPEdata,       &
                                          op_writeLowTempTables

  implicit none

# include "Opacity.h"

  integer :: n
  integer :: status
  integer :: Z
!
!
!    ...Query handling of the Klein-Nishina opacities.
!
!
  call RuntimeParameters_get ("opacity_ignoreKleinNishina",  op_ignoreKleinNishina)
!
!
!   ...Allocate the big arrays:
!
!           1) j-index delimiter array
!                  load the j-index delimiter array -> determine overall maximum
!           2) photoelectric cross section A(i,j,4) coefficient array
!           3) photoelectronic energy range array
!
!
  allocate (op_Jmax (1:op_maxElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_Jmax () allocate failed')
  end if

  call op_setPEarrayJmax   ()

  op_maxJmax = maxval (op_Jmax)

  allocate (op_Aij4 (1:4,1:op_maxJmax,1:op_maxElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_Aij4 () allocate failed')
  end if

  allocate (op_PEenergyRange (LOW:HIGH,1:op_maxJmax,1:op_maxElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_PEenergyRange () allocate failed')
  end if
!
!
!   ...Load the rest of the data into the big arrays.
!
!
  call op_setPEcoeffsAij4  ()
  call op_setPEenergyRange ()
!
!
!   ...Write out the big arrays (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeCompletePEdata ()
  end if
!
!
!   ...Allocate the small arrays depending on the # of elements:
!
!           1) photoelectric cross section A(i,j,4) coefficient array
!           2) j-index delimiter array
!           3) photoelectronic energy range array
!           4) low temperature cutoff array
!
!
  if ((op_totalElements < 1) .or. (op_totalElements > op_maxElements)) then
       call Driver_abortFlash ('[op_initLowTemp] ERROR: # of atomic elements < 1 or > op_maxElements')
  end if

  allocate (op_elementAij4 (1:4,1:op_maxJmax,1:op_totalElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_elementAij4 () allocate failed')
  end if

  allocate (op_elementJmax (1:op_totalElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_elementJmax () allocate failed')
  end if

  allocate (op_elementPEenergyRange (LOW:HIGH,1:op_maxJmax,1:op_totalElements), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_elementPEenergyRange () allocate failed')
  end if
!
!
!   ...Load the data into the small arrays.
!
!
  do n = 1,op_totalElements

     Z = op_element2AtomicNumber (n)

     if (Z < 1 .or. Z > op_maxElements) then
         call Driver_abortFlash ('[op_initLowTemp] ERROR: Atomic number Z < 1 or > op_maxElements')
     end if

     op_elementJmax                                (n) = op_Jmax                                (Z)
     op_elementAij4               (1:4,1:op_maxJmax,n) = op_Aij4               (1:4,1:op_maxJmax,Z)
     op_elementPEenergyRange (LOW:HIGH,1:op_maxJmax,n) = op_PEenergyRange (LOW:HIGH,1:op_maxJmax,Z)

  end do
!
!
!   ...Deallocate the big arrays.
!
!
  deallocate (op_Aij4)
  deallocate (op_Jmax)
  deallocate (op_PEenergyRange)
!
!
!   ...Write out the small arrays (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeElementsPEdata ()
  end if
!
!
!   ...Allocate the arrays needed for integration. Note the extra +1 spot in these arrays
!      to treat the 'zero' region of the integrals for the photoelectron low energy region.
!      Below a certain radiation energy there is no photoelectron effect and this region
!      must be treated as having zero Opacity.
!
!
  allocate (op_A1group (1:op_maxJmax+1), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_A1group () allocate failed')
  end if

  allocate (op_A2group (1:op_maxJmax+1), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_A2group () allocate failed')
  end if

  allocate (op_A3group (1:op_maxJmax+1), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_A3group () allocate failed')
  end if

  allocate (op_A4group (1:op_maxJmax+1), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_A4group () allocate failed')
  end if

  allocate (op_intLimits (1:op_maxJmax+1+1), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_intLimits () allocate failed')
  end if
!
!
!   ...Allocate the low temperature grid array and the low temperature opacity tables.
!
!
  allocate (op_tableLowTemp (1:op_maxNstepsLowTemp), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_tableLowTemp () allocate failed')
  end if

  allocate (op_PlanckLowTempTables (1:op_maxNstepsLowTemp,1:op_totalSpecies,1:op_nEnergyGroups), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_PlanckLowTempTables () allocate failed')
  end if

  allocate (op_RosselandLowTempTables (1:op_maxNstepsLowTemp,1:op_totalSpecies,1:op_nEnergyGroups), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initLowTemp] ERROR: op_RosselandLowTempTables () allocate failed')
  end if
!
!
!   ...Calculate the low temperature opacity tables.
!
!
  call op_calculateLowTempOpacities ()
!
!
!   ...Write out the low temperature opacity tables (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeLowTempTables ()
  end if
!
!
!   ...Set initialization status.
!
!
  op_initializedLowTemp = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initLowTemp
