!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleFinalize
!!
!! NAME
!!
!!  gr_mpoleFinalize
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleFinalize()
!!
!!
!! DESCRIPTION
!!
!!  Finalize the multipole Poisson solver.  Deallocate all storage
!!
!!***

subroutine gr_mpoleFinalize()

  use gr_mpoleData, ONLY : G_2DSPHERICAL, Moment, Momtmp, mpole_geometry,&
                         costable, sintable, rpower, &
                         Legk1, rprinv, Legk2, Leg_fact,&
                         pleg, pint, yzn, xzn, r2, gpot, &
                         MomtmpMatrix, mpole_useMatrixMPI

  implicit none


  deallocate ( Moment)
  
  if (.not. mpole_useMatrixMPI) then
     deallocate( Momtmp)
  else
     deallocate( MomtmpMatrix)
  end if

  if(mpole_geometry /= G_2DSPHERICAL) then
     deallocate( costable)
     deallocate( sintable)
     deallocate( rpower) 
     deallocate( Legk1)
     deallocate( rprinv)
     deallocate( Legk2)
     deallocate( Leg_fact)

  else
     
     deallocate( pleg)
     deallocate( pint)
     deallocate( yzn)
     deallocate( xzn)
     deallocate( r2)
     deallocate( gpot)

  endif


  return
end subroutine gr_mpoleFinalize
