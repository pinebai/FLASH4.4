!!****if* source/physics/Eos/EosMain/Eos_data
!!
!! NAME
!!
!!  Eos_data
!!
!! 
!! SYNOPSIS
!!
!! use Eos_data
!!
!! DESCRIPTION
!!
!!  This is the data module for the Gamma law Eos implementation.
!!  It stores all the runtime parameters, and all the unit scope
!!  data. Some of the unit scope data is fecthed by the wrapper layer
!!  from elsewhere in the code and some is local unit data common
!!  multiple functions in the unit 
!! 
!! PARAMETERS
!!  
!!   These are the runtime parameters used by the Gamma law implementation
!!   of the Eos unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!   gamma[Real]   --- The ideal gas gamma from runtime parameters
!!   smalle[Real]  --- the smallest value for energy (runtime Paramters)
!!
!!***


module Eos_data
#include "Flash.h"
#include "Eos.h"
#include "Eos_map.h"

  real, save :: eos_gasConstant 
  real, save :: eos_smalle
  real, save :: eos_gamma
  real, save :: eos_singleSpeciesA
  real, save :: eos_singleSpeciesZ
!  real, save :: eos_gammam1
  real, save :: eos_eintSwitch
  integer, save :: eos_type, eos_meshMe, eos_meshNumProcs, eos_globalMe

  ! maximum number of iterations for the Newton loop to find T from e
  integer, save :: eos_maxNewton

  ! how accurately to do the Newton iteration to get T from e
  real, save :: eos_tol

  ! force the iterative solver to leave inputs alone (always true in MODE_DENS_TEMP)
  logical, save :: eos_forceConstantInput

  real, save :: eos_smallt
  real, save :: eos_largeT = 1.0e10 ! used by some modified implementations of Newton-Raphson
  real, save :: eos_smallRho

  integer, save :: eos_logLevel = 700

  ! Some stuff that is only used by multiTemp implementations
  integer, save :: eos_combinedTempRule = -1

  integer, save :: eos_entrEleScaleChoice = -1

  real, save :: eos_smallEion=0.0, eos_smallEele=0.0, eos_smallErad=0.0

#ifdef FIXEDBLOCKSIZE
  real,save,dimension(NSPECIES*MAXCELLS) :: eos_massFr
  real,save,dimension(EOS_NUM*MAXCELLS) :: eos_inOut
#else
  real,save, allocatable :: eos_inOut(:),eos_massFr(:)
#endif

  real,parameter :: eos_pradScaleFactor = 1.0 !hardwired for now
#ifdef FLLM_VAR
  integer,parameter :: eos_pradScaleVar = FLLM_VAR
#else
  integer,parameter :: eos_pradScaleVar = -1
#endif

  integer, save, dimension(1:EOSMAP_NUM_ROLES, 1:2, 1:5) :: eos_mapLookup
  logical, save :: eos_threadWithinBlock = .false.
end module Eos_data
