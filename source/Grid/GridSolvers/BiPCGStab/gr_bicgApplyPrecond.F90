!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bicgApplyPrecond
!!
!! NAME
!!
!!  gr_bicgApplyPrecond
!!
!! SYNOPSIS
!!
!!  call gr_bicgApplyPrecond(integer(in) :: isoln,
!!                           integer(in) :: isrc,
!!                           integer, dimension(6)(in) :: bc_types,
!!                           real, dimension(2,6)(in) :: bc_values)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isoln : 
!!
!!   isrc : 
!!
!!   bc_types : 
!!
!!   bc_values : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!
!
! gr_bicgApplyPrecond:
!
! Subroutine that calls preconditioner solves for BiPCGStab.
! In this default case no preconditioner is used, therefore
! The preconditioning matrix M is the identity.
!
!

  subroutine gr_bicgApplyPrecond(isoln, isrc, bc_types, bc_values)

#include "Flash.h"

  use Grid_interface,    ONLY : Grid_getListOfBlocks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr

  use gr_bicgData, ONLY : gcellflg, ili, iui, jli, jui, kli, kui

  implicit none
#include "constants.h"

  integer, intent(in) :: isoln, isrc
  integer, dimension(6), intent(in) :: bc_types
  real, dimension(2,6), intent(in)  :: bc_values

  ! Local vars:
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer :: i, j, k, lb, blockID

  real, pointer, dimension(:,:,:,:) :: unk

  ! Get list of processors leaf blocks:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do lb = 1, blockCount
     blockID = blockList(lb)
     ! Point to blocks center vars:
     call Grid_getBlkPtr(blockID,unk,CENTER)
     do k=kli,kui
        do j=jli,jui
           do i=ili,iui
              unk(isoln,i,j,k)=unk(isrc,i,j,k)
           enddo
        enddo
     enddo
     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,unk,CENTER)
  enddo

  ! We need to fill guardcells before A*x operation in BiPCGStab 
  gcellflg = .true.

  return
  end subroutine gr_bicgApplyPrecond
