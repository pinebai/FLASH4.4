!!****f* source/physics/sourceTerms/Stir/Stir_init
!!
!! NAME
!!  Stir_init
!!
!! SYNOPSIS
!!
!!  call Stir_init()
!!            logical(in) :: restart)
!!
!! DESCRIPTION
!!  Apply the isothermal cooling and stirring opperator 
!!  on the list of blocks provided as input
!!
!! ARGUMENTS
!!   
!!  
!!   restart -indicates if run is starting from scratch or restarting from 
!!            checkpoint
!!
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in the Stir unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!    st_decay [REAL]
!!        correlation time for driving
!!    st_energy [REAL]
!!        energy input/mode
!!    st_freq [INTEGER]
!!        frequency of stirring
!!    st_seed [INTEGER]
!!        random number generator seed
!!    st_stirmax [REAL]
!!        maximum stirring *wavenumber*
!!    st_stirmin [REAL]
!!        minimum stirring *wavenumber*
!!    useStir [BOOLEAN]
!!        runtime switch to turn Stir on/off
!!***

subroutine Stir_init(restart)
  implicit none
  
  
  logical, intent(in) :: restart

  return
end subroutine Stir_init
