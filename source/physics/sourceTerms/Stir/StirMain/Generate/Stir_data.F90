!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/Stir_data
!!
!! NAME
!!
!!  Stir_data
!!
!! SYNOPSIS
!!  Stir_data()
!!
!! DESCRIPTION
!!  Stores the local data for Source Term: Stir
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
!!    st_computeDt {BOOLEAN]
!!        whether to restrict timestep based on stirring
!!    useStir [BOOLEAN]
!!        Switch to turn stirring on or off at runtime.
!!
!!
!!***

Module Stir_data

  integer,parameter :: st_maxmodes = 1000

  integer, save :: st_meshMe, st_meshComm

  !OU variance corresponding to decay time and energy input rate
  real,save :: st_OUvar
 
  !Number of modes
  integer, save :: st_nmodes

  real,save, dimension(3,st_maxmodes) :: st_mode, st_aka, st_akb
  real,save, dimension(6*st_maxmodes) ::st_OUphases

  logical,save :: st_useStir, st_computeDt

  real,save :: st_decay
  real,save :: st_energy
  real,save :: st_stirmin, st_stirmax
  integer,save  :: st_seed  
  integer,save :: st_freq

  integer,save,dimension(:), allocatable :: st_randseed 
  integer, save :: st_seedLen

  integer, save :: st_eosMode
  logical, save :: st_reproducible, st_saveReproducible

  integer, save :: st_randomSaveUnit
  
end Module Stir_data
