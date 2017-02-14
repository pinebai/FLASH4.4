!!****if* source/physics/Diffuse/localAPI/diff_advanceTherm
!!
!!  NAME 
!!
!!  diff_advanceTherm
!!
!!  SYNOPSIS
!!
!!  call diff_advanceTherm(integer(IN) :: blockCount,
!!                            integer(IN) :: blockList(blockCount),
!!                            real(IN)    :: dt)
!!
!!  DESCRIPTION 
!!      This routine advances the heat diffusion equation (i.e., heat conduction).
!!      An implicit scheme is used. 
!!
!!      Supported boundary conditions are ??? and
!!      periodic (1).  The same boundary conditions are currently applied
!!      in all directions.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which the solution must be updated
!!  dt           - The time step
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
!! NOTES
!!
!!  This is a dummmy implementation; it will get invoked when the common
!!  implementation of Diffuse in the main subunit, DiffuseMain, gets called,
!!  but none of the subdirectories of DiffuseMain, which contain actual
!!  implementation of standalone diffusion solvers, got included in the
!!  simulation.  That configuation is not necessarily a mistake - it is
!!  the expected situation when the Diffuse unit gets included only for
!!  the DiffuseFluxBased routines.
!!
!!  The interface of this subroutine must be explicitly know to code that
!!  calls it.  The simplest way to make it so is to have something like
!!     use diff_interface,ONLY: diff_advanceTherm
!!  in the calling routine.
!!***

subroutine diff_advanceTherm(blockCount,blockList,dt,pass)

  implicit none

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in) :: dt
  integer, OPTIONAL, intent(IN):: pass

end subroutine diff_advanceTherm
