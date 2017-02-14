!!****if* source/physics/RadTrans/RadTransMain/RadTrans_sumEnergy
!!
!!  NAME 
!!
!!  RadTrans_sumEnergy
!!
!!  SYNOPSIS
!!
!!  call RadTrans_sumEnergy( integer(IN) :: ivar,
!!                           integer(IN) :: nblk,
!!                           integer(IN) :: blklst(nblk),
!!
!!  DESCRIPTION 
!!
!!    This subroutine is useful when mesh replication is active
!!    (meshCopyCount > 1). It takes an unk variable and adds it across
!!    all the meshes. For example, if each mesh had computed ERAD
!!    separately using its own energy groups, this subroutine could be
!!    used to add all of the ERADs to compute the total radiation
!!    energy.
!!
!! ARGUMENTS
!!
!!   ivar   : the unk variable to add
!!   blklst : the list of blocks to cover
!!   nblk   : the number of blocks
!!
!!***
subroutine RadTrans_sumEnergy(ivar, nblk, blklst)
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits
  use RadTrans_data, ONLY: rt_meshCopyCount, rt_acrossComm
  use Timers_interface, ONLY : Timers_start, Timers_stop 
    implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer, intent(in) :: ivar
  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)

  ! Local Variables:
  integer :: lb, i, j, k, is, js, ks, ni, nj, nk
  integer :: blkLimitsGC(LOW:HIGH,MDIM), blkLimits(LOW:HIGH,MDIM)    
  integer :: ierr
  real, pointer :: blkPtr(:,:,:,:)
  real, allocatable :: local_var(:,:,:)
  real, allocatable :: global_var(:,:,:)

  if(rt_meshCopyCount == 1) return

  call Timers_start("RadTrans_sumEnergy")
  
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     ni = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
     nj = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
     nk = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
     is = blkLimits(LOW, IAXIS)
     js = blkLimits(LOW, JAXIS)
     ks = blkLimits(LOW, KAXIS)

     allocate(local_var(ni, nj, nk))
     allocate(global_var(ni, nj, nk))

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)              
              local_var(i-is+1, j-js+1, k-ks+1) = blkPtr(ivar,i,j,k)
           end do
        end do
     end do

     call mpi_allreduce(local_var, global_var, ni*nj*nk, FLASH_REAL, & 
          MPI_SUM, rt_acrossComm, ierr )

     if(ierr /= MPI_SUCCESS) then
        call Driver_abortFlash("MPI_ALLREDUCE Error in Simulation initBlock")
     end if

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)              
              blkPtr(ivar,i,j,k) = global_var(i-is+1, j-js+1, k-ks+1)
           end do
        end do
     end do
     
     deallocate(local_var)
     deallocate(global_var)
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

  call Timers_stop("RadTrans_sumEnergy")

end subroutine RadTrans_sumEnergy
