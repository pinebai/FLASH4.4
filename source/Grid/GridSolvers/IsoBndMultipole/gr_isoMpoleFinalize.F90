!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoMpoleFinalize
!!
!! NAME
!!
!!  gr_isoMpoleFinalize
!!
!! 
!! SYNOPSIS
!!
!!  gr_isoMpoleFinalize()
!!
!!
!! DESCRIPTION
!!
!!  Finalize the isolated boundary conditions using multipole
!!  solver.
!!
!!***

subroutine gr_isoMpoleFinalize

  use gr_isoMpoleData, ONLY : Moment, Momtmp, costable,sintable,&
                               rpower, Legk1,rprinv,Legk2,Leg_fact
  implicit none

  deallocate ( Moment)
  deallocate( Momtmp)
  deallocate( costable)
  deallocate( sintable)
  deallocate( rpower)
  deallocate( Legk1)
  deallocate( rprinv)
  deallocate( Legk2)
  deallocate( Leg_fact)

  return
end subroutine gr_isoMpoleFinalize
