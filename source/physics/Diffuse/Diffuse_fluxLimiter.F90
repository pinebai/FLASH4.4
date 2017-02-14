!!****f* source/physics/Diffuse/Diffuse_fluxLimiter
!!
!!  NAME 
!!
!!  Diffuse_fluxLimiter
!!
!!  SYNOPSIS
!!
!! 
!!  call Diffuse_fluxLimiter(integer(in) :: idcoef, 
!!                           integer(in) :: ifunc,
!!                           integer(in) :: ifl,
!!                           integer(in) :: mode,
!!                           integer(IN) :: blkcnt,
!!                           integer(IN) :: blklst(blkcnt))
!!
!!  DESCRIPTION 
!!      This routine modifies the diffusion coefficient (idcoef) 
!!      and applies limiting.
!!
!! ARGUMENTS
!!
!!   idcoef  : index into solution vector, indicating a variable that holds
!!             the coefficient to which the limiter is applied.
!!             This variable is used for input and output.
!!   ifunc   : index into solution vector giving the quantity whose flux
!!             is to be limited.
!!             This variable is used only for input.
!!   ifl     : index into solution vector giving flux limiter variable.
!!             This variable is used for input and output.
!!   mode    : Flux limiter mode.
!!   blkcnt  : The number of blocks in the list
!!   blklst  : The list of blocks on which the solution must be updated.
!!
!! SIDE EFFECTS
!!
!! Modifies the diffusion coefficient variable, indicated by the
!! argument idcoef.
!!
!! Modifies the flux limiter variable, indicated by the argument ifl.
!! On output it will be overwritten with a value three times the flux
!! limiter factor lambda, which should be familiar from the literature
!! on flux limiters.  The value of a flux limiter factor lambda should
!! be in the range ( 0, 1/3 ].  Thus the value returned in variable
!! ifl of the solution vector will be in the range [ 0.0, 1.0 ].
!!
!! NOTES
!!
!!  Results for idcoef and ifl are returned for interior cells
!!  and additionally one layer of guard cells.
!!
!!  Valid input data in ifunc is required in one layers of guard
!!  cells.  The caller is responsible for ensuring that this
!!  requirement is met, for example by calling Grid_fillGuardCells
!!  before calling this routine.
!!***
subroutine Diffuse_fluxLimiter(idcoef, ifunc, ifl, mode, blkcnt, blklst)
  implicit none
  
  ! Arguments:
  integer, intent(in) :: idcoef
  integer, intent(in) :: ifunc
  integer, intent(in) :: ifl
  integer, intent(in) :: mode
  integer, intent(IN) :: blkcnt
  integer, intent(IN) :: blklst(blkcnt)

end subroutine Diffuse_fluxLimiter
