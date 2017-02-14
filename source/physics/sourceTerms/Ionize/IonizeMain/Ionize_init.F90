!!****if* source/physics/sourceTerms/Ionize/IonizeMain/Ionize_init
!!
!! NAME
!!  
!!  Ionize_init
!!
!!
!! SYNOPSIS
!! 
!!  call Ionize_init()
!!
!!  
!! DESCRIPTION
!!
!!  Perform various initializations (apart from the problem-dependent ones)
!!  for the Ionize unit.
!!
!!
!! ARGUMENTS
!!
!!   
!!
!!***

subroutine Ionize_init()
  use ion_interface, ONLY : ion_readTable
  use Ionize_data,ONLY : ion_tneimin, ion_tneimax, &
       ion_dneimin, ion_dneimax, ion_smallx, ion_useIonize, &
       ion_emass

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  
#include "constants.h"
  
  implicit none
   
  

  

  call RuntimeParameters_get("smallx", ion_smallx)  
  call RuntimeParameters_get("tneimin", ion_tneimin)
  call RuntimeParameters_get("tneimax", ion_tneimax)
  call RuntimeParameters_get("dneimin", ion_dneimin)
  call RuntimeParameters_get("dneimax", ion_dneimax)
  call RuntimeParameters_get("useIonize",ion_useIonize)

  call PhysicalConstants_get("electron mass",ion_emass)

  !Read ionize coefficients and related data into the module level arrays: 
  !ion_nion12, ion_tp, ion_cfinz, ion_cfric.
  call ion_readTable() 

end subroutine Ionize_init
