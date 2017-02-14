!!****if* source/Grid/GridMain/paramesh/Paramesh2/tot_bnd
!! NAME
!!  tot_bnd
!!
!! SYNOPSIS
!!
!!  tot_bnd(integer(IN) :: idiag
!!          integer(IN) :: idir)
!!  
!! DESCRIPTION 
!!
!!  This routine is called by amr_guardcell in Paramesh2 to 
!!  apply boundary condition to "unk", which is represented by 
!!  the constant "CENTER" in FLASH. Instead of calculating the boundary
!!  conditions directly, this routine makes the necessary checks
!!  and prepares the blocks for another routine "gr_bcApplyToAllBlks", one 
!!  block at a time. It then passed control to that routine, which
!!  applies boundary conditions appropriately
!! 
!! ARGUMENTS
!!  
!!    idiag -   This parameters is redundant in FLASH. It is kept
!!              only for Paramesh2 interface compatibility. In other
!!              codes it indicates whether to fill corner guardcells
!!              In FLASH it is always assumed to be true.
!!   idir   -   Direction to update:  either ALLDIR (for
!!              all directions) or #AXIS, where #=I,J or K
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_BOUNDARY
#endif

subroutine tot_bnd (idiag, idir)
  use tree, ONLY : lnblocks,neigh
  use Grid_data,ONLY : gr_blkList,gr_allPeriodic,gr_meshMe,gr_bndOrder,gr_domainBC
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks
  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: idiag, idir
  integer lb,i,j,n

  integer :: varCount, loc_idir
  logical :: isWork = .false.
  
  if(gr_allPeriodic) then
     return
  end if

  if (idir == 0) then
     do i = 0,NDIM-1
        loc_idir = gr_bndOrder(NDIM-i)
        call gr_bcApplyToAllBlks(loc_idir,isWork)
     end do
  else
     loc_idir=idir
     call gr_bcApplyToAllBlks(loc_idir,isWork)
  end if

  return
end subroutine tot_bnd


!!****if* source/Grid/GridMain/paramesh/Paramesh2/tot_bnd_work
!! NAME
!!  tot_bnd_work
!!
!! SYNOPSIS
!!
!!  tot_bnd_work(integer(IN) :: idiag
!!          integer(IN) :: idir)
!!  
!! DESCRIPTION 
!!
!!  This routine is called by amr_guardcell in Paramesh2 to 
!!  apply boundary condition to "work", which is represented by 
!!  the constant "WORK" in FLASH. Instead of calculating the boundary
!!  conditions directly, this routine makes the necessary checks
!!  and prepares the blocks for another routine "gr_bcApplyToAllBlks", one 
!!  block at a time. It then passed control to that routine, which
!!  applies boundary conditions appropriately
!! 
!! ARGUMENTS
!!  
!!    idiag -   This parameters is redundant in FLASH. It is kept
!!              only for Paramesh2 interface compatibility. In other
!!              codes it indicates whether to fill corner guardcells
!!              In FLASH it is always assumed to be true.
!!    idir   -   Direction to update:  either ALLDIR (for
!!              all directions) or #AXIS, where #=I,J or K
!!
!!***



subroutine tot_bnd_work (idiag, idir)

  use tree, ONLY : lnblocks,neigh
  use Grid_data,ONLY : gr_allPeriodic,gr_bndOrder
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks
  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: idiag, idir
  integer lb,i,j,n

  integer :: varCount, loc_idir
  logical :: isWork=.true.

  if(gr_allPeriodic) then
     return
  end if
  
  if (idir == 0) then
     do i = 0,NDIM-1
        loc_idir = gr_bndOrder(NDIM-i)
        call gr_bcApplyToAllBlks(loc_idir,isWork)
     end do
  else
     loc_idir=idir
     call gr_bcApplyToAllBlks(loc_idir,isWork)
  end if
  return
end subroutine tot_bnd_work




