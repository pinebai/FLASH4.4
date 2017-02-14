!!****if* source/physics/Diffuse/DiffuseMain/Diffuse
!!
!!  NAME 
!!
!!  Diffuse
!!
!!  SYNOPSIS
!!
!!  call Diffuse( integer(IN) :: blockCount,
!!                integer(IN) :: blockList(blockCount),
!!                real(IN)    :: dt,
!!       OPTIONAL,integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This routine advances the solution for diffusive effects,
!!      such as the heat diffusion equation (i.e., heat conduction).
!!      Typically, implicit schemes will be used where available. 
!!
!!      Currently, only support for advancing the heat diffusion equation
!!      (i.e., heat conduction) is available in Diffuse implementations.
!!      As of FLASH 3.3, a flux-based implementation in connection with
!!      a Hydro solver must be used for other diffusive effects, in particular
!!      viscosity and magnetic resistivity.
!!
!!      Supported boundary conditions are ??? and
!!      periodic (1).  The same boundary conditions are currently applied
!!      in all directions.
!!
!!      DEV: Actually, boundary condition support depends on the solver implementation.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which the solution must be updated
!!   dt           : The time step
!!   pass         : Reverses the direction of solve in between half time steps
!!                  i.e X-Y-Z, Z-Y-X
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  updated temperature and some auxiliaries.  Invokes a solver (of the Heat
!!  diffusion equation). On return,
!!     TEMP_VAR:  contains updated temperature for the current simulation time.
!!     EINT_VAR, ENER_VAR, PRES_VAR, etc.: updated accordingly by EOS call
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Diffuse implementation.
!!  The following information is subject to change.
!!     COND_VAR:  contains conductivity that was passed to Grid_advanceDiffusion
!!
!!
!!***

!!REORDER(4): solnVec

subroutine Diffuse(blockCount,blockList,dt,pass)


  use diff_interface, ONLY : diff_advanceTherm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Diffuse_data, ONLY : useDiffuse

  implicit none

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in) :: dt
  integer, OPTIONAL, intent(IN):: pass

!=========================================================================
  if(.not.useDiffuse) return

  call Timers_start("Diffuse")

  call diff_advanceTherm(blockCount,blockList,dt,pass)

  ! DEV: Add more calls here as standalone implementations for viscosity, etc. become available. - KW

  call Timers_stop ("Diffuse")
  
  return
end subroutine Diffuse
