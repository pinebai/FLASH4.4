!!****if* source/Grid/localAPI/gr_mpolePotentials
!!
!! NAME
!!
!!  gr_mpolePotentials
!!
!! SYNOPSIS
!!
!!  gr_mpolePotentials  (integer, intent(in) :: ipotvar,
!!                       real,    intent(in) :: Poisson_factor )
!!
!! DESCRIPTION
!!
!!  Computes the potential field using the mass moments already
!!  calculated. On output tha variable indexed by ipotvar contains
!!  the potential. The calculations are entirely local to each
!!  processor, since each processor has a local copy of the moments.
!!
!!  This routine calls the appropriate subroutines according to
!!  the geometry specified.
!!
!! ARGUMENTS
!!
!!  ipotvar        : index to variable containing the potential
!!  Poisson_factor : the name says it all 
!!
!!***

subroutine gr_mpolePotentials (ipotvar,Poisson_factor)

  implicit none
    
  integer, intent (in) :: ipotvar
  real,    intent (in) :: Poisson_factor

  return
end subroutine gr_mpolePotentials
