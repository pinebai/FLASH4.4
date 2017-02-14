!!****if* source/Grid/GridSolvers/Pfft/gr_pfftFinalize
!!
!! NAME
!!
!!  gr_pfftFinalize
!!
!! 
!! SYNOPSIS
!!
!!  gr_pfftFinalize()
!!
!!
!! DESCRIPTION
!!
!!  Finalize the multipole Poisson solver.  Deallocate all storage
!!
!!***
#ifdef DEBUG_ALL
#define DEBUG_PFFT
#endif

#define DEBUG_PFFT

subroutine gr_pfftFinalize()

  use Grid_interface, ONLY : Grid_pfftFinalize
  use gr_pfftData, ONLY : pfft_setupOnce
  use Grid_data, ONLY : gr_meshMe

  implicit none

#ifdef DEBUG_PFFT
  if(gr_meshMe==0) print*,'gr_pfftFinalize: pfft_setupOnce is', pfft_setupOnce
#endif
  if(pfft_setupOnce)call Grid_pfftFinalize()
  
  return
end subroutine gr_pfftFinalize
