!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_setMinTimeStep
!!
!! NAME
!!
!!  hy_uhd_setMinTimeStep
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_setMinTimeStep(integer(in) :: blockID,
!!                        integer(in) :: i,
!!                        integer(in) :: j,
!!                        integer(in) :: k,
!!                        real(in)    :: delta,
!!                        real(in)    :: speed)
!!
!!
!! DESCRIPTION
!!
!!  Set minimum time step and its location internally.
!!
!!
!! ARGUMENTS
!!
!!  i - local index in x direction
!!  j - local index in y direction
!!  k - local index in z direction
!!  blockID - local blockID
!!  speed - local maximum wave speed
!!
!!***



subroutine hy_uhd_setMinTimeStep(blockID,i,j,k,delta,speed)

  use Hydro_data, ONLY : hy_dtmin, hy_dtminloc, hy_dtminValid, hy_dtminCfl, hy_cfl, &
                         hy_meshMe

  implicit none

  !!$ Argument list -------------------------------------
  integer, INTENT(in) :: blockID,i,j,k
  real, INTENT(in) :: delta,speed
  !!$ ---------------------------------------------------
  real :: dt_hydro

  dt_hydro = hy_cfl*delta/abs(speed)

  if( dt_hydro < hy_dtmin) then
    hy_dtmin = dt_hydro
    hy_dtminloc(1) = i
    hy_dtminloc(2) = j
    hy_dtminloc(3) = k
    hy_dtminloc(4) = blockID
    hy_dtminloc(5) = hy_meshMe
    hy_dtminCfl    = hy_cfl
    hy_dtminValid = .TRUE.
  end if

end subroutine hy_uhd_setMinTimeStep
