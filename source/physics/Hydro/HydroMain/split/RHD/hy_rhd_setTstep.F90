!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_setTstep
!!
!! NAME
!!
!!  hy_rhd_setTstep
!!
!! SYNOPSIS
!!
!!  hy_rhd_setTstep( 
!!                  real(IN) :: dtime,
!!                  integer(IN)  :: i,
!!                  integer(IN)  :: j,
!!                  integer(IN)  :: k,
!!                  integer(IN)  :: blockID)
!!
!! DESCRIPTION
!!
!!  Minimum dt is actually computed in hy_rhd_sweep, where all required
!!  information is evaluated anyway. This routine exchanges results
!!  with hy_rhd_sweep through private data members. 
!!
!!
!! ARGUMENTS
!!
!!  dtime  -       timestep calculated in hy_rhd_sweep
!!  i,j,k  -       indices of the location of lowest tstep
!!  blockID -      local block ID
!!
!!
!!***

subroutine hy_rhd_setTstep( dtime,i,j,k,blockID)
  
  use Hydro_data, only : hy_dtmin, hy_dtminloc, hy_meshMe

  implicit none

  real, INTENT(in) :: dtime
  integer, INTENT(in) ::  i,j,k,blockID

  if( dtime <= hy_dtmin ) then
    hy_dtmin = dtime
    hy_dtminloc(1) = i
    hy_dtminloc(2) = j
    hy_dtminloc(3) = k
    hy_dtminloc(4) = blockID
    hy_dtminloc(5) = hy_meshMe
  end if

end subroutine hy_rhd_setTstep



