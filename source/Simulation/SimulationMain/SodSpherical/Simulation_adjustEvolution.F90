!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

#include "constants.h"
#include "Flash.h"

  use Simulation_data

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr
  use Grid_interface, ONLY: Grid_getBlkBC

  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real,    intent(in) :: dt
  real,    intent(in) :: stime

  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer, dimension(2,MDIM) :: bcs
  integer :: axis(MDIM)
  integer :: i
  integer :: j
  integer :: k
  integer :: lb
  real    :: newCfl

  real, pointer :: blkPtr(:,:,:,:)

#ifdef CFL_VAR

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     if (NDIM > 1) call Grid_getBlkBC(blklst(lb),bcs)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              newCfl = 1.0
#if NDIM > 1
              if (bcs(LOW,IAXIS) == REFLECTING .OR. bcs(LOW,IAXIS) == AXISYMMETRIC) then
                 if (i==blkLimits(LOW,IAXIS)) then
                    newCfl = .25
                    if (bcs(LOW,JAXIS) == REFLECTING) then
                       if (j .LE. blkLimits(LOW,JAXIS)+3) then
                          newCfl = .2333 ! 0.275 / 1.2
                       end if
                    end if
                 end if
              end if
#endif

              blkPtr(CFL_VAR,i,j,k) = min(newCfl,blkPtr(CFL_VAR,i,j,k))
           enddo
        end do
     end do

     call Grid_releaseBlkPtr(blklst(lb), blkPtr)

  end do

#endif
  
end subroutine Simulation_adjustEvolution
