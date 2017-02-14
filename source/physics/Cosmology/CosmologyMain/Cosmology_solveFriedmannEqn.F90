!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_solveFriedmannEqn
!!
!! NAME
!!
!!  csm_solveFriedmannEqn
!!
!!
!! SYNOPSIS
!!
!!  Cosmology_solveFriedmannEqn(real,intent(IN) :: timeEndAdvance,
!!                                   real, intent(IN) :: dt)
!!  
!!
!! DESCRIPTION
!! 
!! Numerically solve the Friedmann equation, deriving the scale
!! factor at time t+dt from the scale factor at time t.  This
!!  version assumes a matter-dominated universe.
!! 
!! ARGUMENTS
!! 
!!   timeEndAdvance - the time at the end of last timestep advance
!!   dt             - time step
!!
!!***

subroutine Cosmology_solveFriedmannEqn (timeEndAdvance, dt)


  use Cosmology_data, ONLY : csm_oldscaleFactor, csm_scaleFactor
  implicit none

  real,intent(IN)  :: timeEndAdvance, dt

  real            :: t1, t2, dt_alt
   

! Save the old scale factor and redshift.

  csm_oldScaleFactor = csm_scaleFactor

! Start and end times (to avoid roundoff)

  t1     = timeEndAdvance
  t2     = timeEndAdvance + dt
  dt_alt = t2 - t1
  t2     = t1 + dt_alt

  call csm_integrateFriedmann(csm_oldScaleFactor,timeEndAdvance,dt,csm_scaleFactor)
  
  return
end subroutine Cosmology_solveFriedmannEqn




