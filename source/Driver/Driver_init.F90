!!****f* source/Driver/Driver_init
!!
!! NAME
!!  Driver_init
!!
!! SYNOPSIS
!!  call Driver_init()
!!
!! DESCRIPTION
!!
!!  Perform the Driver unit initializations.
!!  Gets runtime parameters from the flash.par, or from checkpoint file
!!  if run is a restart.
!!  
!!  Initializes the simulation time, initial dt, timestep, whether 
!!  from scratch or restart.
!!  Also, verifies that the initial dt passes CFL criteria.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used by the basic Driver unit.
!!   Your specific implementation may have more runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You may have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    dt_init [REAL]
!!        Initial timestep DEV: inconsistent naming
!!    dtmax [REAL]
!!        Maximum timestep
!!    dtmin [REAL]
!!        Minimum timestep
!!    nbegin [INTEGER]
!!        First timestep
!!    nend [INTEGER]
!!        Maximum number of timesteps to take
!!    restart [BOOLEAN]
!!        Is this a restart run?
!!    tinitial [REAL]
!!        Initial simulation time
!!    tmax [REAL]
!!        Maximum simulation time
!!    wall_clock_time_limit [REAL]
!!        Total wall clock time limit (seconds)
!!    tstep_change_factor [REAL]
!!        factor to allow multiplicative increase in dt until it
!!        hits the CFL condition. This lets initial dt be
!!        very conservative initially, but increase rapidly to find the
!!        the optimum value.
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


subroutine Driver_init()

  
  implicit none

  return
end subroutine Driver_init








