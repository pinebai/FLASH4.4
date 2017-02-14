!!****f* source/Driver/Driver_computeDt
!!
!! NAME
!!
!!  Driver_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDt(integer(IN) :: nbegin,
!!                  integer(IN) :: nstep, 
!!                  real(IN)    :: simTime,
!!                  real(IN)    :: dtOld,
!!                  real(OUT)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Determine the stability-limited time step.
!!  This timestep is determined using information from the included
!!  physics modules - many different timestep limiters are polled.
!!
!!  The global driver might use a different (hopefully smaller) time
!!  step, to match a file write time (tplot or trstr) or if the
!!  simulation end time has been reached; such possibilities are
!!  not considered here.
!!
!! ARGUMENTS
!!  nbegin - first step of the simulation (this is only used
!!              to determine if a label header should be written to
!!              the screen)
!!  nstep - current step of the simulation
!!  simTime - current simulation time of the run
!!  dtOld - the dt from the timestep that we just finished 
!!         (it's old because we be using dtOld to calculate 
!!          and return the dt for the next timestep (dtNew)
!!  dtNew - returned value of the dt calculated for the next timestep
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


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_computeDt(nbegin, nstep, &
                    simTime, dtOld, dtNew)


  implicit none


  integer, intent(IN) :: nbegin, nstep
  real,    intent(IN) :: simTime    !! current simulation time
  real, intent(IN) :: dtOld      !! last time step we used
  real, intent(OUT):: dtNew      !! the new timestep we get. to be returned.
 

  !dummy values for stubs

  dtNew = 0.0

  return
end subroutine Driver_computeDt
