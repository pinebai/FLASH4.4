!!****f* source/physics/Eos/Eos_nucOneZone
!!
!! NAME
!!
!!  Eos_nucOneZone
!!
!! SYNOPSIS
!!
!!  call Eos_nucOneZone(real(INOUT)  :: xdens,
!!                      real(INOUT)  :: xtemp,
!!                      real(IN)  :: xye,
!!                      real(INOUT)  :: xener,
!!                      real(INOUT)  :: xpres,
!!                      real(INOUT)  :: xentr,
!!                      real(OUT)  :: xVar,
!!                      integer(IN)  :: mode)
!!
!! DESCRIPTION
!!
!!  Routine for cranking the EOS on a single zone, in any old
!!  mode you like.  Calls kernel routine nuc_eos_full.  Can also
!!  be used to recover values from the table that are not provided
!!  by nuc_eos_short without doing a full EOS solve.
!!
!! ARGUMENTS
!!
!!   xdens : 
!!
!!   xtemp : 
!!
!!   xye : 
!!
!!   xener : 
!!
!!   xpres : 
!!
!!   xentr : 
!!
!!   xVar : 
!!
!!   mode : Eos mode
!!
!!
!!
!!***

subroutine Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xdedt,xCs2,xXp,xXn,xXa,xXh, xVar,varID,mode)

  implicit none 

  real, intent(INOUT) :: xDens
  real, intent(IN)    :: xYe
  real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
  real, intent(OUT)   :: xVar
  integer, intent(IN) :: mode, varID
  real, intent(OUT) :: xXp, xXn, xXa,xXh,xdedt,xCs2

  ! STUB

  xXp = 0.0;  xXn = 0.0;  xXa = 0.0; xXh = 0.0; xdedt = 0.0; xCs2 = 0.0
  xVar = 0.0
  return
end subroutine Eos_nucOneZone
