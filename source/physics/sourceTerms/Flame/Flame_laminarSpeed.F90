!!****f* source/physics/sourceTerms/Flame/Flame_laminarSpeed
!!
!! NAME
!!
!!  Flame_laminarSpeed
!!
!! SYNOPSIS
!!
!!  call Flame_laminarSpeed (real(in)          :: dens,
!!                           real(out)         :: s ,
!!                          optional,real(in)  :: info(:))
!!
!! DESCRIPTION
!!
!!  Return the physical laminar flame speed s.
!!  Generally only used in initialization.
!!  Behavior is strongly dependent on what FlameSpeed subunit is included
!!  the info array is to provide imlementation-dependent information such
!!  as composition if applicable.
!!
!! ARGUMENTS
!!
!!   dens - density
!!   s    - returned speed
!!
!! SEE ALSO
!!
!!  See Flame_interface.F90 for possible updates
!!
!!***

! this is a stub for when the Flame Unit is not included
!
! Dean Townsley 2008
!
subroutine Flame_laminarSpeed(dens, s, ds, info)

  implicit none

  real, intent(in)   :: dens
  real, intent(out)  :: s
  real, optional, intent(out)  :: ds
  real, dimension(:), optional, intent(in) :: info

  s=0.0
  if (present(ds)) ds=0.0

end subroutine Flame_laminarSpeed
