!!****if* source/Simulation/SimulationMain/MGDInfinite/Simulation_adjustEvolution
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

  use Simulation_data, ONLY: sim_fileUnit
  use Simulation_data, ONLY: sim_globalMe

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr

  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real,    intent(in) :: dt
  real,    intent(in) :: stime

  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: i
  integer :: j
  integer :: k
  integer :: lb

  real, pointer :: blkPtr(:,:,:,:)

  if(sim_globalMe /= MASTER_PE) return

  ! Get the temperatures in any cell:
  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! Write temperatures (and whatever else) to a file. If
              ! you add additional quantities here, you might want to
              ! also change the file header in Simulation_init.F90:
              write(sim_fileUnit, '(i10,5(1pe15.6))') &
                   nstep, stime, dt, &
                   blkPtr(TION_VAR,i,j,k), &
                   blkPtr(TELE_VAR,i,j,k), &
                   blkPtr(TRAD_VAR,i,j,k)
              exit
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     exit
  end do
  
end subroutine Simulation_adjustEvolution
