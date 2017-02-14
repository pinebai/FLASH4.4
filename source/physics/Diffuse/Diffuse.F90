!!****f* source/physics/Diffuse/Diffuse
!!
!!  NAME 
!!
!!  Diffuse
!!
!!  SYNOPSIS
!!
!!  call Diffuse( integer(IN)           :: blockCount,
!!                integer(IN)           :: blockList(blockCount),
!!                real(IN)              :: dt,
!!                integer(IN), OPTIONAL :: pass)
!!
!!  DESCRIPTION 
!!      This routine advances the heat diffusion equation (i.e., heat conduction).
!!      Typically, implicit schemes will be used where available. 
!!
!!      Currently, only support for advancing the heat diffusion equation
!!      (i.e., heat conduction) is available in Diffuse implementations.
!!      As of FLASH 3.3, a flux-based implementation in connection with
!!      a Hydro solver must be used for other diffusive effects, in particlular
!!      viscosity and magnetic resistivity.
!!
!!      Supported boundary conditions are ??? and
!!      periodic (1).  The same boundary conditions are currently applied
!!      in all directions.
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
!!     WTMP_VAR:  contains source term that was passed to Grid solver routines
!!     COND_VAR:  contains conductivity that was passed to Grid solver routines
!!     DTMP_VAR:  contains temperature increment for the current simulation time.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!
!!
!!***

subroutine Diffuse(blockCount,blockList,dt,pass)

  implicit none

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in) :: dt
  integer, OPTIONAL, intent(IN):: pass

end subroutine Diffuse
