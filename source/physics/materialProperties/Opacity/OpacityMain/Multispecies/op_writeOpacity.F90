!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_writeOpacity
!!
!! NAME
!!
!!  op_writeOpacity
!!
!! SYNOPSIS
!!
!!  call op_writeOpacity ()
!!
!! DESCRIPTION
!!
!!  Prints out details of the main Opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeOpacity ()

  use Opacity_data,   ONLY : op_nEnergyGroups,            &
                             op_totalSpecies,             &
                             op_totalElements,            &
                             op_atomName,                 &
                             op_atomWeight,               &
                             op_element2AtomicNumber,     &
                             op_elementNumberofSpecies,   &
                             op_elements2Species,         &
                             op_speciesNumberofElements,  &
                             op_species2FractionElements, &
                             op_speciesLowTempCutoff,     &
                             op_speciesWeights,           &
                             op_species2elements,         &
                             op_energyGroupBoundaries

  use Timers_interface, ONLY : Timers_start, Timers_stop 
  implicit none

  integer :: group
  integer :: element
  integer :: fileUnit, ut_getFreeFileUnit
  integer :: n
  integer :: nElements
  integer :: nSpecies
  integer :: species
  integer :: Z

  real    :: cutoff
  real    :: Elower,Eupper
  real    :: fraction
  real    :: weight
!
!
!    ...Open the printout file.
!
!
  call Timers_start("op_writeOpacity")

  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout.txt", &
        form = 'formatted')
!
!
!   ...Print out the title.
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   OPACITY UNIT PRINTOUT"
  write (fileUnit,*)
!
!
!   ...Print out the elements. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) ' Elements / Atomic weights'
  write (fileUnit,*)

  do element = 1,op_totalElements
     Z = op_element2AtomicNumber (element)
     write (fileUnit,'(8X,A12,F12.4)') op_atomName (Z),op_atomWeight (Z)
  end do
!
!
!   ...Print out the species -> elements composition. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) ' Species -> Elements / Number fraction composition'
  write (fileUnit,*)

  do species = 1,op_totalSpecies

     nElements = op_speciesNumberofElements  (species)
     element   = op_species2elements         (species,1)
     Z         = op_element2AtomicNumber     (element)
     fraction  = op_species2FractionElements (species,element)

     write (fileUnit,'(1X,I3,A4,A12,F14.8)') species,' -> ',op_atomName (Z),fraction

     do n = 2,nElements

        element   = op_species2elements         (species,n)
        Z         = op_element2AtomicNumber     (element)
        fraction  = op_species2FractionElements (species,element)

        write (fileUnit,'(8X,A12,F14.8)') op_atomName (Z),fraction
     end do

     write (fileUnit,*)

  end do
!
!
!   ...Print out the species -> elements index map. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) ' Species -> Elements (index map)'
  write (fileUnit,*)

  do species = 1,op_totalSpecies
     nElements = op_speciesNumberofElements  (species)
     write (fileUnit,'(1X,I3,A4,80I3)') species,' -> ', &
                                        (op_species2elements (species,element),element=1,nElements)
  end do
!
!
!   ...Print out the elements -> species index map. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) ' Elements -> Species (index map)'
  write (fileUnit,*)

  do element = 1,op_totalElements
     nSpecies = op_elementNumberofSpecies (element)
     write (fileUnit,'(1X,I3,A4,80I3)') element,' -> ', &
                                        (op_elements2Species (element,species),species=1,nSpecies)
  end do
!
!
!   ...Print out the species -> weight / low temperature cutoff. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) ' Species -> Weight / Low Temperature Cutoff (in Kelvin)'
  write (fileUnit,*)

  do species = 1,op_totalSpecies

     weight = op_speciesWeights       (species)
     cutoff = op_speciesLowTempCutoff (species)

     write (fileUnit,'(1X,I3,A4,2ES14.6)') species,' -> ',weight,cutoff

  end do
!
!
!   ...Print out the energy group boundaries. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) ' Energy Group -> Lower Bound (eV) / Upper Bound (eV)'
  write (fileUnit,*)

  do group = 1,op_nEnergyGroups

     Elower = op_energyGroupBoundaries (group)
     Eupper = op_energyGroupBoundaries (group+1)

     write (fileUnit,'(1X,I6,A4,2ES20.12)') group,' -> ',Elower,Eupper

  end do
!
!
!    ...Close the printout file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!
  call Timers_stop("op_writeOpacity")
  return
end subroutine op_writeOpacity
