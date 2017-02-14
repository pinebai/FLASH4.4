!!****f* source/Driver/Driver_initFlash
!!
!! NAME
!!  Driver_initFlash
!!
!! SYNOPSIS
!!
!!   Driver_initFlash()
!!
!! DESCRIPTION
!!
!!  Performs Flash initializations, which includes:
!!
!!  Call all 'init' routines in units, in the right order.
!!  Order does matter, particularly when restarting from a 
!!  checkpoint file.
!!
!!  For the most part, Driver_initFlash calls another units' init
!!  routines directly, like call IO_init or call Grid_init.  This
!!  routine also makes calls to other Driver initialization routines
!!  like Driver_initMaterialProperties or Driver_initSourceTerms.
!!  These routines then call the unit-specific initialization 
!!  routines.  This level of abstraction was added to simplify
!!  the initialization calls.
!!
!!
!!
!!***


!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"


subroutine Driver_initFlash()

  
implicit none
  return
end subroutine Driver_initFlash
