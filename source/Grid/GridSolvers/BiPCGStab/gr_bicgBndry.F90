!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bicgBndry
!!
!! NAME
!!
!!  gr_bicgBndry
!!
!! SYNOPSIS
!!
!!  call gr_bicgBndry(integer(in) :: ivar,
!!                    integer(in) :: nlayers,
!!                    logical(in) :: gcellflg)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ivar : 
!!
!!   nlayers : 
!!
!!   gcellflg : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     bicg_bndry()

!  Description: Update boundary zones on a given level for a given variable.

!  Parameters:  
!               ivar        Index of variable to update.
!               nlayers     Number of Layers to fill guardcells.
!               gcelflg     If 0 just copy ivar to work, if one copy to work and fill guardcells.

subroutine gr_bicgBndry(ivar, nlayers, gcellflg)

!===============================================================================
#include "Flash.h"

  use gr_bicgData, ONLY: ili, jli, kli, iui, jui, kui
  
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_fillGuardCells

  use tree, only : nodetype,lrefine
  use workspace, ONLY: work
  use paramesh_dimensions, only : nguard_work
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs

  use Driver_interface, ONLY : Driver_abortFlash

implicit none
#include "constants.h"


integer, intent(in) :: ivar, nlayers
logical, intent(in) :: gcellflg

integer :: i, j, k, lb, lnblocks

real,    pointer, dimension(:,:,:,:), save :: solnData

!!$logical, save :: first_call=.true.
!!$integer, parameter :: ndim = NDIM
!!$integer, parameter :: k2d = K2D
!!$integer, parameter :: k3d = K3D


!===============================================================================


call Grid_getLocalNumBlks(lnblocks)

#if NDIM == 2

do lb = 1, lnblocks
    if (nodetype(lb) .eq. 1) then
 
       ! Point to blocks center vars:
       call Grid_getBlkPtr(lb,solnData,CENTER)
            
       ! Copy data to work adding one layer of guardcells.
       do k = kli, kui
          do j = jli-1, jui+1
             do i = ili-1, iui+1
                work(i,j,k,lb,1) = solnData(ivar,i,j,k)
             enddo
          enddo
       enddo

       ! Release pointers:
       call Grid_releaseBlkPtr(lb,solnData,CENTER)
    endif
enddo

#elif NDIM == 3

do lb = 1, lnblocks
    if (nodetype(lb) .eq. 1) then
 
       ! Point to blocks center vars:
       call Grid_getBlkPtr(lb,solnData,CENTER)
            
       ! Copy data to work adding one layer of guardcells.
       do k = kli-1, kui+1
          do j = jli-1, jui+1
             do i = ili-1, iui+1
                work(i,j,k,lb,1) = solnData(ivar,i,j,k)
             enddo
          enddo
       enddo

       ! Release pointers:
       call Grid_releaseBlkPtr(lb,solnData,CENTER)
    endif
enddo

#endif

if (gcellflg) then
  

   if (nlayers > nguard_work) &
   call Driver_abortFlash("bicg_Bndry: you must set nlayers <= nguard_work")

   ! Fills only one layer of guardcells from Work Data structure, enough 
   ! to complete 2nd order difference centered stencil on block boundaries. 
   call Grid_fillGuardCells( WORK, ALLDIR, minLayers=1)

endif

!===============================================================================

return
end subroutine gr_bicgBndry
