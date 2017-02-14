!!****f* source/physics/Diffuse/Diffuse_computeFluxLimiter
!!
!!  NAME 
!!
!!  Diffuse_computeFluxLimiter
!!
!!  SYNOPSIS
!!
!! 
!!  call Diffuse_computeFluxLimiter(integer(in) :: idcoef, 
!!                           integer(in) :: ifunc,
!!                           integer(in) :: ifl,
!!                           integer(in) :: iflOut,
!!                           integer(in) :: ieddi3,
!!                           integer(in) :: mode,
!!                           real(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:),
!!                           integer(in),value :: lbIU,lbUJ,lbUK,
!!                           integer(IN) :: blockID,
!!                  OPTIONAL,integer(IN) :: gcLayers )
!!
!!  DESCRIPTION 
!!      This routine computes flux limiter factors and returns them
!!      in UNK variables. It does not apply the limiter by modifying
!!      diffusion coefficients or in any other way.
!!
!! ARGUMENTS
!!
!!   idcoef  : index into solution vector, indicating a variable that holds
!!             the coefficient to which the limiter is applied.
!!             This variable is used only for input.
!!   ifunc   : index into solution vector giving the quantity whose flux
!!             is to be limited. This gives either the quantity directly
!!             if iDenForIfunc is < 1, or must be multiplied by the variable
!!             given by iDenForIfunc.
!!             This variable is used only for input.
!!   iDenForIfunc : A density to multiply ifunc by, or -1.
!!             This variable is used only for input.
!!   ifl     : index into solution vector giving the limiting value that
!!             would be applied to the variable given by idcoef such that
!!             the latter is less than or equal to the former.
!!             This variable is used only for input.
!!   iflOut  : index into solution vector to put flux limiter factor.
!!             Note that this may be one of the variables given by the
!!             previous four dummy arguments that are used as inputs.
!!   mode     : Flux limiter mode.
!!   solnData : The block data on which to operate
!!   lbUI,lbUJ,lbUK: lower bounds of assumed-shape array solnData
!!   blockID  : Identifies the block on which to operate
!!   gcLayers : In how many layers of guard cells, in additionm to the
!!              interior cells, the flux limiter coefficients are desired.
!!              Note that one more layer is required in the variable given
!!              by ifl for input, for computing the gradient.
!!
!! SIDE EFFECTS
!!
!!  Modifies flux limiter variables (as requested by argument iflOut) in the solution
!!  storage UNK.  On output it will be overwritten with a value
!! three times the flux limiter factor lambda, which should be
!! familiar from the literature on flux limiters. The value
!! of lambda should be in the range [ 0, 1/3 ].
!!
!!***

#include "FortranLangFeatures.fh"

subroutine Diffuse_computeFluxLimiter(idcoef, ifunc, iDenForIfunc, ifl, iflOut, ieddi3, mode, &
     solnData, lbUI,lbUJ,lbUK, &
     blockID, gcLayers)
  implicit none
  
  ! Arguments:
  integer, intent(in) :: idcoef
  integer, intent(in) :: ifunc
  integer, intent(in) :: iDenForIfunc
  integer, intent(in) :: ifl
  integer, intent(in) :: iflOut
  integer, intent(in) :: ieddi3
  integer, intent(in) :: mode
  integer, VALUE_INTENT(IN) :: lbUI,lbUJ,lbUK
  real,    intent(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:)
  integer, intent(IN) :: blockID
  integer, intent(IN),OPTIONAL :: gcLayers

end subroutine Diffuse_computeFluxLimiter
