!!****ih* source/physics/utilities/PlasmaState/localAPI/pls_interface
!!
!! NAME
!!
!!  pls_interface
!!
!! SYNOPSIS
!!
!!   use pls_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the PlasmaState unit that defines its
!!  private interfaces.
!!
!!***

Module pls_interface

  interface
     subroutine pls_initComposition () 
     end subroutine pls_initComposition
  end interface

  interface
     subroutine pls_setSpeciesElementsData () 
     end subroutine pls_setSpeciesElementsData
  end interface

end Module pls_interface
