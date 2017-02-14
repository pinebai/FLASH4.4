!!****f* source/Particles/Particles_sinkSumAttributes
!!
!! NAME
!!
!!  Particles_sinkSumAttributes
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSumAttributes(  real(OUT) :: sums(:),
!!                                   integer(in) :: attribs(:),
!!                          OPTIONAL,integer(in) :: factor)
!!
!! DESCRIPTION
!!
!!  Compute global sums over all sink particles for the attributes in
!!  the list ATTRIBS, optionally multiplied by the value of another
!!  particle attribute FACTOR.
!!
!!
!! ARGUMENTS
!!
!!   sums - the sums are returned here
!!
!!   attribs - list of sink particle properties
!!
!!   factor  - a sink particle property to use as an optional  factor
!!
!! NOTES
!!
!!   written by Klaus Weide, 2014
!!
!!***

subroutine Particles_sinkSumAttributes(sums, attribs, factor)

  implicit none

  real,intent(OUT)   :: sums(:)
  integer,intent(in) :: attribs(:)
  integer,intent(in),OPTIONAL :: factor

  sums(:) = 0.0

  return

end subroutine Particles_sinkSumAttributes
