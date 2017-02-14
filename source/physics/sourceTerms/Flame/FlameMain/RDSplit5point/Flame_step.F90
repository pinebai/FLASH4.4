!!****if* source/physics/sourceTerms/Flame/FlameMain/RDSplit5point/Flame_step
!!
!! NAME
!!
!!  Flame_step
!!
!! SYNOPSIS
!!
!!  call Flame_step(integer(in) :: num_blocks,
!!                  integer(in) :: blocklist(num_blocks),
!!                  real(in)    :: dt)
!!
!! DESCRIPTION
!!
!!  see Flame_interface.F90 for function description
!!
!! ARGUMENTS
!!
!!   num_blocks : 
!!
!!   blocklist : 
!!
!!   dt : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! see Flame_interface.F90 at top level for function description
!
! Dean Townsley 2008
!

! Implementation details
!
! This implementation is operator split betwee reaction and diffusion
! but dimensionally unsplit.  Diffusion is treated by directly computing
! the Laplacian instead of differencing fluxes, so that this is not a
! conservative diffusion operator.
!
! Computation of the Laplacian itself is done in a subroutine, as its
! form depends heavily on the mesh geometry

!!REORDER(4): solnData

#define DEBUG_GRID_GCMASK

#include "Flash.h"
#include "constants.h"
subroutine Flame_step( num_blocks, blockList, dt  )    

  use Flame_data

  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkIndexLimits, Grid_fillGuardCells
  use fl_fsInterface, only : fl_flameSpeed
  use fl_effInterface, only: fl_effects
  use fl_interface, only : fl_laplacian
  use Driver_interface, only : Driver_abortFlash
  use Timers_interface, only : Timers_start, Timers_stop
  use Logfile_interface, only : Logfile_stamp, Logfile_stampVarMask
       
  implicit none
  integer, INTENT(in)                        :: num_blocks
  integer, INTENT(in), DIMENSION(num_blocks) :: blockList
  real,    INTENT(in)                        :: dt

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, fspeedLimits
  integer :: istat
  integer :: n, bid

  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:,:,:) :: flam, flamdot, flamespeed, lapl

  real :: f, inv_dt
  integer :: i,j,k
  integer :: sizeI, sizeJ, sizeK

  if( .not. fl_useFlame ) return

  call Timers_start("flame")

  inv_dt = 1.0/dt

  if (fl_gcDoLogMask) then
     call Logfile_stamp('calling guardcell fill with mask logging on','[Flame_step]')
#ifdef DEBUG_GRID_GCMASK
     call Logfile_stampVarMask(fl_gcMask, fl_gcDoEos, '[Flame_step]', 'gcMask')
#endif
  end if
  call Grid_fillGuardCells(CENTER, ALLDIR, eosMode=MODE_DENS_EI, &
    doEos=fl_gcDoEos, maskSize=fl_gcMaskSize, mask=fl_gcMask, &
    makeMaskConsistent=.false., doLogMask=fl_gcDoLogMask)
  fl_gcDoLogMask=.false.

  do n = 1,num_blocks
     bid = blockList(n)

     call Grid_getBlkPtr(bid,solnData)

     call Grid_getBlkIndexLimits(bid, blkLimits, blkLimitsGC)

     sizeI=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     sizeJ=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     sizeK=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

     allocate(flam(sizeI,sizeJ,sizeK), STAT=istat)
     if (istat /= 0) call Driver_abortFlash("Cannot allocate flam in Flame_step")
     allocate(flamdot(sizeI,sizeJ,sizeK), STAT=istat)
     if (istat /= 0) call Driver_abortFlash("Cannot allocate flamdot in Flame_step")
     allocate(flamespeed(sizeI,sizeJ,sizeK), STAT=istat)
     if (istat /= 0) call Driver_abortFlash("Cannot allocate flamespeed in Flame_step")
     allocate(lapl(sizeI,sizeJ,sizeK), STAT=istat)
     if (istat /= 0) call Driver_abortFlash("Cannot allocate lapl in Flame_step")

     ! extract flam variable, should make cache work better
     ! need two layers in GCs becausee of RD splitting
     flam( blkLimits(LOW,IAXIS)-2 : blkLimits(HIGH,IAXIS)+2 , &
           blkLimits(LOW,JAXIS)-2*K2D : blkLimits(HIGH,JAXIS)+2*K2D , &
           blkLimits(LOW,KAXIS)-2*K3D : blkLimits(HIGH,KAXIS)+2*K3D) &
         = solnData(FLAM_MSCALAR, blkLimits(LOW,IAXIS)-2 : blkLimits(HIGH,IAXIS)+2 , &
                 blkLimits(LOW,JAXIS)-2*K2D : blkLimits(HIGH,JAXIS)+2*K2D , &
                 blkLimits(LOW,KAXIS)-2*K3D : blkLimits(HIGH,KAXIS)+2*K3D)


     call fl_flameSpeed(solnData, flamespeed, bid, 2)

     do k = blkLimits(LOW,KAXIS)-2*K3D, blkLimits(HIGH,KAXIS)+2*K3D
        do j = blkLimits(LOW,JAXIS)-2*K2D, blkLimits(HIGH,JAXIS)+2*K2D
           do i = blkLimits(LOW,IAXIS)-2, blkLimits(HIGH,IAXIS)+2
              f = flam(i,j,k)
              flam(i,j,k) = f + dt*fl_R_over_s*flamespeed(i,j,k)*(f-fl_epsilon_0)*(1.0+fl_epsilon_1-f)
           enddo
        enddo
     enddo

     ! 1 specifies the step size should be 1 grid cell
     ! cannot be any larger because flam is filled with only 2 guard cell layers
     call fl_laplacian(lapl, flam, 1, bid)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              flam(i,j,k) = max(0.0, min(1.0, flam(i,j,k) + dt*fl_kappa_over_s*flamespeed(i,j,k)*lapl(i,j,k) ) )
              flamdot(i,j,k) = (flam(i,j,k) - solnData(FLAM_MSCALAR,i,j,k))*inv_dt
              solnData(FLAM_MSCALAR, i,j,k) = flam(i,j,k)
           enddo
        enddo
     enddo

     deallocate(lapl)
     deallocate(flamespeed)
     deallocate(flam)

     call fl_effects( solnData, flamdot, dt, bid)

     deallocate(flamdot)

     call Grid_releaseBlkPtr(bid, solnData)

  enddo

  call Timers_stop("flame")

  return
end subroutine Flame_step
