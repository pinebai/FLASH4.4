!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_setTimestep
!!
!! NAME
!!
!!  hy_8wv_setTimestep
!!
!!
!! SYNOPSIS
!!
!!  hy_8wv_setTimestep(integer(in) :: hy_meshMe,
!!                     real(inout) :: dtmin,
!!                     integer(in) :: i,
!!                     integer(in) :: j,
!!                     integer(in) :: k,
!!                     integer(in) :: blockID)
!!
!!
!! DESCRIPTION
!!
!!  Set minimum time step and its location internally.
!!
!!
!! ARGUMENTS
!!
!!  hy_meshMe - local processor number
!!  dtmin - local minimum time step
!!  i - local index in x direction
!!  j - local index in y direction
!!  k - local index in z direction
!!  blockID - local blockID
!!
!!***



subroutine hy_8wv_setTimestep(hy_meshMe,dtmin,i,j,k,blockID)

  use Hydro_data, only : hy_dtmin, hy_dtminloc

  implicit none

  !!$ Argument list -------------------------------------
  real, INTENT(in) :: dtmin
  integer, INTENT(in) :: i,j,k,blockID, hy_meshMe
  !!$ ---------------------------------------------------

  if( dtmin <= hy_dtmin ) then
    hy_dtmin = dtmin
    hy_dtminloc(1) = i
    hy_dtminloc(2) = j
    hy_dtminloc(3) = k
    hy_dtminloc(4) = blockID
    hy_dtminloc(5) = hy_meshMe
  end if

end subroutine hy_8wv_setTimestep
