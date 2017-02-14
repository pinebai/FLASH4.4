!!****if* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CONe/Flame_laminarSpeed
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
!! Dean Townsley 2008
!!
!!  In this implementation the info array contains the Carbon 12 and Neon 22 aboundances
!!  in the fuel as mass fractions.  The rest is assumed to be Oxygen 16.
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
!!
!!***


subroutine Flame_laminarSpeed(dens, s, ds, info)

  use Flame_data, ONLY: fl_useFlame

  implicit none

  real, intent(in)   :: dens
  real, intent(out)  :: s
  real, optional, intent(out) :: ds
  real, dimension(:), optional, intent(in) :: info


  real :: c12i, ne22i, flame_width

  if (present(info)) then
     c12i = info(1)
     ne22i = info(2)
  else
     ! default values
     c12i = 0.5
     ne22i = 0.0
  endif

  if (fl_useFlame) then
     call fl_fsConeInterp(c12i,ne22i,dens,s,flame_width)
  else
     s = 0.0
  end if

  if (present(ds)) ds = flame_width

end subroutine Flame_laminarSpeed
