!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_minConfinement1DRec
!!
!! NAME
!!
!!  ed_minConfinement1DRec
!!
!! SYNOPSIS
!!
!!  ed_minConfinement1DRec (integer, intent (in) :: nc,
!!                          real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the minimum confinement rules for 1D rectangular ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 1 and <= size (y), not checked here!)
!!  y  : dependent variable array (rx,vx,P)
!!
!!***

function ed_minConfinement1DRec (nc,y)

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: ed_minConfinement1DRec (1:nc)    ! (compatible to Runge Kutta interfaces)

  ed_minConfinement1DRec (1:nc) = 0.0

  return
end function ed_minConfinement1DRec
