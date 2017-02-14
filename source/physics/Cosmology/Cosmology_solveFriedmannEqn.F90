!!****f* source/physics/Cosmology/Cosmology_solveFriedmannEqn
!!
!! NAME
!!
!!  Cosmology_solveFriedmannEqn
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

  implicit none

  real,intent(IN)  :: timeEndAdvance, dt
  
  return

end subroutine Cosmology_solveFriedmannEqn
