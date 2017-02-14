!!****if* source/physics/sourceTerms/Flame/FlameSpeed/Constant/Flame_laminarSpeed
!!
!! NAME
!!
!!  Flame_laminarSpeed
!!
!! SYNOPSIS
!!
!!  call Flame_laminarSpeed(real(in) :: dens,
!!                          real(out) :: s,
!!                          real, optional(out) :: ds,
!!                          real, dimension(:), optional(in) :: info)
!!
!! DESCRIPTION
!!
!!
!! Dean Townsley 008
!!
!! ARGUMENTS
!!
!!   dens : 
!!
!!   s : 
!!
!!   ds : 
!!
!!   info : 
!!
!!
!!***


subroutine Flame_laminarSpeed(dens, s, ds, info)

  use fl_fsData

  implicit none

  real, intent(in)   :: dens
  real, intent(out)  :: s
  real, optional, intent(out) :: ds
  real, dimension(:), optional, intent(in) :: info

  s = fl_fsConstFlameSpeed

  if (present(ds)) ds = fl_fsConstFlameWidth

end subroutine Flame_laminarSpeed
