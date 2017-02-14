!!****f* source/physics/Eos/Eos_init
!!
!! NAME
!!
!!  Eos_init
!!
!! 
!! SYNOPSIS
!!
!!  call Eos_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars needed
!!  by the EOS unit from the runtime parameters and physical
!!  constants facilities
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!  
!!   Particular implementations (Gamma,Helmholz,etc) of the unit
!!   define their own runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might overwrite these values with the 
!!   flash.par values for your specific run.  
!!
!!
!!***


subroutine Eos_init()
    
    implicit none
    ! stub for eos initialization.  A non-stub implementation of his routine
    ! will be supplied, if required,  by the main subunit of Eos, normally
    ! located under Eos/EosMain.

    return
end subroutine Eos_init


