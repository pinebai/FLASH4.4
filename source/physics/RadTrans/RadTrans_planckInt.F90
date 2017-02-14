!!****f* source/physics/RadTrans/RadTrans_planckInt
!!
!!  NAME 
!!
!!  RadTrans_planckInt
!!
!!  SYNOPSIS
!!
!!  call RadTrans_planckInt(real(IN)  :: x,
!!                           real(OUT) :: p)
!!
!!  DESCRIPTION 
!!      This routine must evaluate the Planck integral
!!
!!      The Planck integral is:
!!      
!!             x/|        y^3
!!      P(x) =   |  dy ----------
!!             0 |/    exp(y) - 1
!!
!!
!! ARGUMENTS
!!
!!   x : The argument for the planck integral
!!   p : The value of the integral
!!
!!***
subroutine RadTrans_planckInt(x, p)

  implicit none

  ! Arguments:
  real, intent(in)  :: x
  real, intent(out) :: p

  ! Stub implementation
  p = 0.0
  return
end subroutine RadTrans_planckInt
